"""A pile contains subpiles which contain tracesfiles which contain traces."""

import cPickle as pickle
import copy
import logging
import math
import operator
import os
import re
import sys
import time
import weakref

import config
import util

import io
import trace


def sl(s):
    return map(str, sorted(list(s)))


class Counter(dict):
    def __missing__(self, k):
        return 0

    def update(self, other):
        for k, v in other.iteritems():
            self[k] += v

    def subtract(self, other):
        for k, v in other.iteritems():
            self[k] -= v
            if self[k] <= 0:
                del self[k]

    def subtract1(self, k):
        self[k] -= 1
        if self[k] <= 0:
            del self[k]


import avl

pjoin = os.path.join
logger = logging.getLogger("pyrocko.pile")

from trace import degapper


def avl_remove_exact(avltree, element):
    ilo, ihi = avltree.span(element)
    for i in xrange(ilo, ihi):
        if avltree[i] is element:
            avltree.remove_at(i)
            return

    raise ValueError("avl_remove_exact(avltree, element): element not in avltree")


def cmpfunc(key):
    if isinstance(key, str):
        # special cases; these run about 50% faster than the generic one on
        # Python 2.5
        if key == "tmin":
            return lambda a, b: cmp(a.tmin, b.tmin)
        if key == "tmax":
            return lambda a, b: cmp(a.tmax, b.tmax)

        key = operator.attrgetter(key)

    return lambda a, b: cmp(key(a), key(b))


class Dummy:
    pass


class Sorted(object):
    def __init__(self, values=[], key=None):
        self._set_key(key)
        self._avl = avl.new(values, self._cmp)

    def _set_key(self, key):
        self._key = key
        self._cmp = cmpfunc(key)
        if isinstance(key, str):

            class Dummy:
                def __init__(self, k):
                    setattr(self, key, k)

            self._dummy = Dummy

    def __getstate__(self):
        state = list(self._avl.iter()), self._key
        return state

    def __setstate__(self, state):
        l, key = state
        self._set_key(key)
        self._avl = avl.from_iter(iter(l), len(l))

    def insert(self, value):
        self._avl.insert(value)

    def remove(self, value):
        avl_remove_exact(self._avl, value)

    def insert_many(self, values):
        for value in values:
            self._avl.insert(value)

    def remove_many(self, values):
        for value in values:
            avl_remove_exact(self._avl, value)

    def __iter__(self):
        return iter(self._avl)

    def with_key_in(self, kmin, kmax):
        omin, omax = self._dummy(kmin), self._dummy(kmax)
        ilo, ihi = self._avl.span(omin, omax)
        return self._avl[ilo:ihi]

    def min(self):
        return self._avl.min()

    def max(self):
        return self._avl.max()

    def __len__(self):
        return len(self._avl)


class TracesFileCache(object):
    """Manages trace metainformation cache.

    For each directory with files containing traces, one cache file is
    maintained to hold the trace metainformation of all files which are
    contained in the directory.
    """

    caches = {}

    def __init__(self, cachedir):
        """Create new cache.

        :param cachedir: directory to hold the cache files.

        """

        self.cachedir = cachedir
        self.dircaches = {}
        self.modified = set()
        util.ensuredir(self.cachedir)

    def get(self, abspath):
        """Try to get an item from the cache.

        :param abspath: absolute path of the object to retrieve

        :returns: a stored object is returned or None if nothing could be found.

        """

        dircache = self._get_dircache_for(abspath)
        if abspath in dircache:
            return dircache[abspath]
        return None

    def put(self, abspath, tfile):
        """Put an item into the cache.

        :param abspath: absolute path of the object to be stored
        :param tfile: object to be stored
        """

        cachepath = self._dircachepath(abspath)
        # get lock on cachepath here
        dircache = self._get_dircache(cachepath)
        dircache[abspath] = tfile
        self.modified.add(cachepath)

    def dump_modified(self):
        """Save any modifications to disk."""

        for cachepath in self.modified:
            self._dump_dircache(self.dircaches[cachepath], cachepath)
            # unlock

        self.modified = set()

    def clean(self):
        """Weed out missing files from the disk caches."""

        self.dump_modified()

        for fn in os.listdir(self.cachedir):
            try:
                i = int(fn)  # valid filenames are integers
                cache = self._load_dircache(pjoin(self.cachedir, fn))
                self._dump_dircache(cache, pjoin(self.cachedir, fn))

            except ValueError:
                pass

    def _get_dircache_for(self, abspath):
        return self._get_dircache(self._dircachepath(abspath))

    def _get_dircache(self, cachepath):
        if cachepath not in self.dircaches:
            if os.path.isfile(cachepath):
                self.dircaches[cachepath] = self._load_dircache(cachepath)
            else:
                self.dircaches[cachepath] = {}

        return self.dircaches[cachepath]

    def _dircachepath(self, abspath):
        cachefn = "%i" % abs(hash(os.path.dirname(abspath)))
        return pjoin(self.cachedir, cachefn)

    def _load_dircache(self, cachefilename):

        f = open(cachefilename, "r")
        cache = pickle.load(f)
        f.close()

        # weed out files which no longer exist
        for fn in cache.keys():
            if not os.path.isfile(fn):
                del cache[fn]

        for v in cache.values():
            v.trees_from_content(v.traces)
            for tr in v.traces:
                tr.file = v

            v.data_use_count = 0
            v.data_loaded = False
        return cache

    def _dump_dircache(self, cache, cachefilename):

        if not cache:
            if os.path.exists(cachefilename):
                os.remove(cachefilename)
            return

        # make a copy without the parents and the binsearch trees
        cache_copy = {}
        for fn in cache.keys():
            trf = copy.copy(cache[fn])
            trf.parent = None
            trf.by_tmin = None
            trf.by_tmax = None
            trf.by_tlen = None
            trf.by_mtime = None
            trf.data_use_count = 0
            trf.data_loaded = False
            trf.traces = [tr.copy(data=False) for tr in trf.traces]
            for tr in trf.traces:
                tr.file = trf
            cache_copy[fn] = trf

        tmpfn = cachefilename + ".%i.tmp" % os.getpid()
        f = open(tmpfn, "w")

        pickle.dump(cache_copy, f)
        f.close()
        os.rename(tmpfn, cachefilename)


def get_cache(cachedir):
    """Get global TracesFileCache object for given directory."""
    if cachedir not in TracesFileCache.caches:
        TracesFileCache.caches[cachedir] = TracesFileCache(cachedir)

    return TracesFileCache.caches[cachedir]


def loader(
    filenames,
    fileformat,
    cache,
    filename_attributes,
    show_progress=True,
    update_progress=None,
):
    class Progress:
        def __init__(self, label, n):
            self._label = label
            self._n = n
            self._bar = None
            if show_progress:
                self._bar = util.progressbar(label, self._n)

            if update_progress:
                update_progress(label, 0, self._n)

        def update(self, i):
            if self._bar:
                if i < self._n - 1:
                    self._bar.update(i)
                else:
                    self._bar.finish()
                    self._bar = None

            if update_progress:
                update_progress(self._label, i, self._n)

    if not filenames:
        logger.warn("No files to load from")
        return

    regex = None
    if filename_attributes:
        regex = re.compile(filename_attributes)

    progress = Progress("Looking at files", len(filenames))

    failures = []
    to_load = []
    for i, filename in enumerate(filenames):
        try:
            abspath = os.path.abspath(filename)

            substitutions = None
            if regex:
                m = regex.search(filename)
                if not m:
                    raise FilenameAttributeError(
                        "Cannot get attributes with pattern '%s' from path '%s'"
                        % (filename_attributes, filename)
                    )
                substitutions = {}
                for k in m.groupdict():
                    if k in ("network", "station", "location", "channel"):
                        substitutions[k] = m.groupdict()[k]

            mtime = os.stat(filename)[8]
            tfile = None
            if cache:
                tfile = cache.get(abspath)
            to_load.append(
                (
                    (not tfile or tfile.mtime != mtime or substitutions),
                    mtime,
                    abspath,
                    substitutions,
                    tfile,
                )
            )

        except (OSError, FilenameAttributeError) as xerror:
            failures.append(abspath)
            logger.warn(xerror)
        progress.update(i + 1)

    progress.update(len(filenames))

    to_load.sort()

    nload = len([1 for x in to_load if x[0]])
    iload = 0

    count_all = False
    if nload < 0.01 * len(to_load):
        nload = len(to_load)
        count_all = True

    if to_load:
        progress = Progress("Scanning files", nload)

        for (mustload, mtime, abspath, substitutions, tfile) in to_load:
            try:
                if mustload:
                    tfile = TracesFile(
                        None,
                        abspath,
                        fileformat,
                        substitutions=substitutions,
                        mtime=mtime,
                    )
                    if cache and not substitutions:
                        cache.put(abspath, tfile)

                    if not count_all:
                        iload += 1

                if count_all:
                    iload += 1

            except (io.FileLoadError, OSError) as xerror:
                failures.append(abspath)
                logger.warn(xerror)
            else:
                yield tfile

            progress.update(iload + 1)

        progress.update(nload)

    if failures:
        logger.warn(
            "The following file%s caused problems and will be ignored:\n"
            % util.plural_s(len(failures))
            + "\n".join(failures)
        )

    if cache:
        cache.dump_modified()


def tlen(x):
    return x.tmax - x.tmin


class TracesGroup(object):

    """Trace container base class.

    Base class for Pile, SubPile, and TracesFile, i.e. anything containing
    a collection of several traces. A TracesGroup object maintains lookup sets
    of some of the traces meta-information, as well as a combined time-range
    of its contents.
    """

    def __init__(self, parent):
        self.parent = parent
        self.empty()
        self.nupdates = 0
        self.abspath = None

    def set_parent(self, parent):
        self.parent = parent

    def get_parent(self):
        return self.parent

    def empty(self):
        (
            self.networks,
            self.stations,
            self.locations,
            self.channels,
            self.nslc_ids,
            self.deltats,
        ) = [Counter() for x in range(6)]
        self.by_tmin = Sorted([], "tmin")
        self.by_tmax = Sorted([], "tmax")
        self.by_tlen = Sorted([], tlen)
        self.by_mtime = Sorted([], "mtime")
        self.tmin, self.tmax = None, None
        self.deltatmin, self.deltatmax = None, None

    def trees_from_content(self, content):
        self.by_tmin = Sorted(content, "tmin")
        self.by_tmax = Sorted(content, "tmax")
        self.by_tlen = Sorted(content, tlen)
        self.by_mtime = Sorted(content, "mtime")
        self.adjust_minmax()

    def add(self, content):

        if isinstance(content, trace.Trace) or isinstance(content, TracesGroup):
            content = [content]

        for c in content:

            if isinstance(c, TracesGroup):
                self.networks.update(c.networks)
                self.stations.update(c.stations)
                self.locations.update(c.locations)
                self.channels.update(c.channels)
                self.nslc_ids.update(c.nslc_ids)
                self.deltats.update(c.deltats)

                self.by_tmin.insert_many(c.by_tmin)
                self.by_tmax.insert_many(c.by_tmax)
                self.by_tlen.insert_many(c.by_tlen)
                self.by_mtime.insert_many(c.by_mtime)

            elif isinstance(c, trace.Trace):
                self.networks[c.network] += 1
                self.stations[c.station] += 1
                self.locations[c.location] += 1
                self.channels[c.channel] += 1
                self.nslc_ids[c.nslc_id] += 1
                self.deltats[c.deltat] += 1

                self.by_tmin.insert(c)
                self.by_tmax.insert(c)
                self.by_tlen.insert(c)
                self.by_mtime.insert(c)

        self.adjust_minmax()

        self.nupdates += 1
        self.notify_listeners("add")

        if self.parent is not None:
            self.parent.add(content)

    def remove(self, content):

        if isinstance(content, trace.Trace) or isinstance(content, TracesGroup):
            content = [content]

        for c in content:

            if isinstance(c, TracesGroup):
                self.networks.subtract(c.networks)
                self.stations.subtract(c.stations)
                self.locations.subtract(c.locations)
                self.channels.subtract(c.channels)
                self.nslc_ids.subtract(c.nslc_ids)
                self.deltats.subtract(c.deltats)

                self.by_tmin.remove_many(c.by_tmin)
                self.by_tmax.remove_many(c.by_tmax)
                self.by_tlen.remove_many(c.by_tlen)
                self.by_mtime.remove_many(c.by_mtime)

            elif isinstance(c, trace.Trace):
                self.networks.subtract1(c.network)
                self.stations.subtract1(c.station)
                self.locations.subtract1(c.location)
                self.channels.subtract1(c.channel)
                self.nslc_ids.subtract1(c.nslc_id)
                self.deltats.subtract1(c.deltat)

                self.by_tmin.remove(c)
                self.by_tmax.remove(c)
                self.by_tlen.remove(c)
                self.by_mtime.remove(c)

        self.adjust_minmax()

        self.nupdates += 1
        self.notify_listeners("remove")

        if self.parent is not None:
            self.parent.remove(content)

    def relevant(self, tmin, tmax, group_selector=None, trace_selector=None):

        if not self.by_tmin or not self.is_relevant(tmin, tmax, group_selector):
            return []

        return [
            tr
            for tr in self.by_tmin.with_key_in(tmin - self.tlenmax, tmax)
            if tr.is_relevant(tmin, tmax, trace_selector)
        ]

    def adjust_minmax(self):
        if self.by_tmin:
            self.tmin = self.by_tmin.min().tmin
            self.tmax = self.by_tmax.max().tmax
            t = self.by_tlen.max()
            self.tlenmax = t.tmax - t.tmin
            self.mtime = self.by_mtime.max().mtime
            deltats = self.deltats.keys()
            self.deltatmin = min(deltats)
            self.deltatmax = max(deltats)
        else:
            self.tmin = None
            self.tmax = None
            self.tlenmax = None
            self.mtime = None
            self.deltatmin = None
            self.deltatmax = None

    def notify_listeners(self, what):
        pass

    def get_update_count(self):
        return self.nupdates

    def overlaps(self, tmin, tmax):
        # return not (tmax < self.tmin or self.tmax < tmin)
        return tmax >= self.tmin and self.tmax >= tmin

    def is_relevant(self, tmin, tmax, group_selector=None):
        # return  not (tmax <= self.tmin or self.tmax < tmin) and (selector is
        # None or selector(self))
        return (
            tmax >= self.tmin
            and self.tmax >= tmin
            and (group_selector is None or group_selector(self))
        )


class MemTracesFile(TracesGroup):

    """This is needed to make traces without an actual disc file to be inserted
    into a Pile."""

    def __init__(self, parent, traces):
        TracesGroup.__init__(self, parent)
        self.add(traces)
        self.mtime = time.time()

    def add(self, traces):
        if isinstance(traces, trace.Trace):
            traces = [traces]

        for tr in traces:
            tr.file = self

        TracesGroup.add(self, traces)

    def load_headers(self, mtime=None):
        pass

    def load_data(self):
        pass

    def use_data(self):
        pass

    def drop_data(self):
        pass

    def reload_if_modified(self):
        return False

    def iter_traces(self):
        for trace in self.by_tmin:
            yield trace

    def get_traces(self):
        return list(self.by_tmin)

    def gather_keys(self, gather, selector=None):
        keys = set()
        for trace in self.by_tmin:
            if selector is None or selector(trace):
                keys.add(gather(trace))

        return keys

    def __str__(self):

        s = "MemTracesFile\n"
        s += "file mtime: %s\n" % util.time_to_str(self.mtime)
        s += "number of traces: %i\n" % len(self.by_tmin)
        s += "timerange: %s - %s\n" % (
            util.time_to_str(self.tmin),
            util.time_to_str(self.tmax),
        )
        s += "networks: %s\n" % ", ".join(sl(self.networks.keys()))
        s += "stations: %s\n" % ", ".join(sl(self.stations.keys()))
        s += "locations: %s\n" % ", ".join(sl(self.locations.keys()))
        s += "channels: %s\n" % ", ".join(sl(self.channels.keys()))
        s += "deltats: %s\n" % ", ".join(sl(self.deltats.keys()))
        return s


class TracesFile(TracesGroup):
    def __init__(self, parent, abspath, format, substitutions=None, mtime=None):
        TracesGroup.__init__(self, parent)
        self.abspath = abspath
        self.format = format
        self.traces = []
        self.data_loaded = False
        self.data_use_count = 0
        self.substitutions = substitutions
        self.load_headers(mtime=mtime)
        self.mtime = mtime

    def load_headers(self, mtime=None):
        logger.debug("loading headers from file: %s" % self.abspath)
        if mtime is None:
            self.mtime = os.stat(self.abspath)[8]

        self.remove(self.traces)
        for tr in io.load(
            self.abspath,
            format=self.format,
            getdata=False,
            substitutions=self.substitutions,
        ):
            self.traces.append(tr)
            tr.file = self

        self.add(self.traces)

        self.data_loaded = False
        self.data_use_count = 0

    def load_data(self, force=False):
        file_changed = False
        if not self.data_loaded or force:
            logger.debug("loading data from file: %s" % self.abspath)

            for itr, tr in enumerate(
                io.load(
                    self.abspath,
                    format=self.format,
                    getdata=True,
                    substitutions=self.substitutions,
                )
            ):
                if itr < len(self.traces):
                    xtr = self.traces[itr]
                    if (
                        xtr.mtime != tr.mtime
                        or xtr.tmin != tr.tmin
                        or xtr.tmax != tr.tmax
                    ):
                        logger.warn(
                            "file may have changed since last access (trace number %i has changed): %s"
                            % (itr, self.abspath)
                        )
                        self.remove(xtr)
                        self.traces.remove(xtr)
                        xtr.file = None
                        self.traces.append(tr)
                        self.add(tr)
                        tr.file = self
                        file_changed = True
                    else:
                        xtr.ydata = tr.ydata

                else:
                    self.traces.add(tr)
                    self.add(tr)
                    logger.warn(
                        "file may have changed since last access (new trace found): %s"
                        % self.abspath
                    )
                    file_changed = True
            self.data_loaded = True
        return file_changed

    def use_data(self):
        if not self.data_loaded:
            raise Exception("Data not loaded")
        self.data_use_count += 1

    def drop_data(self):
        if self.data_loaded:
            if self.data_use_count == 1:
                logger.debug("forgetting data of file: %s" % self.abspath)
                for tr in self.traces:
                    tr.drop_data()

                self.data_loaded = False

            self.data_use_count -= 1
        else:
            self.data_use_count = 0

    def reload_if_modified(self):
        mtime = os.stat(self.abspath)[8]
        if mtime != self.mtime:
            logger.debug("mtime=%i, reloading file: %s" % (mtime, self.abspath))
            self.mtime = mtime
            if self.data_loaded:
                self.load_data(force=True)
            else:
                self.load_headers()

            return True

        return False

    def iter_traces(self):
        for trace in self.traces:
            yield trace

    def gather_keys(self, gather, selector=None):
        keys = set()
        for trace in self.by_tmin:
            if selector is None or selector(trace):
                keys.add(gather(trace))

        return keys

    def __str__(self):
        s = "TracesFile\n"
        s += "abspath: %s\n" % self.abspath
        s += "file mtime: %s\n" % util.time_to_str(self.mtime)
        s += "number of traces: %i\n" % len(self.traces)
        s += "timerange: %s - %s\n" % (
            util.time_to_str(self.tmin),
            util.time_to_str(self.tmax),
        )
        s += "networks: %s\n" % ", ".join(sl(self.networks.keys()))
        s += "stations: %s\n" % ", ".join(sl(self.stations.keys()))
        s += "locations: %s\n" % ", ".join(sl(self.locations.keys()))
        s += "channels: %s\n" % ", ".join(sl(self.channels.keys()))
        s += "deltats: %s\n" % ", ".join(sl(self.deltats.keys()))
        return s


class FilenameAttributeError(Exception):
    pass


class SubPile(TracesGroup):
    def __init__(self, parent):
        TracesGroup.__init__(self, parent)
        self.files = []
        self.empty()

    def add_file(self, file):
        self.files.append(file)
        file.set_parent(self)
        self.add(file)

    def remove_file(self, file):
        self.files.remove(file)
        file.set_parent(None)
        self.remove(file)

    def remove_files(self, files):
        for file in files:
            self.files.remove(file)
            file.set_parent(None)
        self.remove(files)

    def gather_keys(self, gather, selector=None):
        keys = set()
        for file in self.files:
            keys |= file.gather_keys(gather, selector)

        return keys

    def iter_traces(
        self,
        load_data=False,
        return_abspath=False,
        group_selector=None,
        trace_selector=None,
    ):
        for file in self.files:

            if group_selector and not group_selector(file):
                continue

            must_drop = False
            if load_data:
                file.load_data()
                file.use_data()
                must_drop = True

            for trace in file.iter_traces():
                if trace_selector and not trace_selector(trace):
                    continue

                if return_abspath:
                    yield file.abspath, trace
                else:
                    yield trace

            if must_drop:
                file.drop_data()

    def iter_files(self):
        for file in self.files:
            yield file

    def reload_modified(self):
        modified = False
        for file in self.files:
            modified |= file.reload_if_modified()

        return modified

    def __str__(self):
        s = "SubPile\n"
        s += "number of files: %i\n" % len(self.files)
        s += "timerange: %s - %s\n" % (
            util.time_to_str(self.tmin),
            util.time_to_str(self.tmax),
        )
        s += "networks: %s\n" % ", ".join(sl(self.networks.keys()))
        s += "stations: %s\n" % ", ".join(sl(self.stations.keys()))
        s += "locations: %s\n" % ", ".join(sl(self.locations.keys()))
        s += "channels: %s\n" % ", ".join(sl(self.channels.keys()))
        s += "deltats: %s\n" % ", ".join(sl(self.deltats.keys()))
        return s


class Pile(TracesGroup):
    def __init__(self):
        TracesGroup.__init__(self, None)
        self.subpiles = {}
        self.open_files = {}
        self.listeners = []
        self.abspaths = set()

    def add_listener(self, obj):
        self.listeners.append(weakref.ref(obj))

    def notify_listeners(self, what):
        for ref in self.listeners:
            obj = ref()
            if obj:
                obj.pile_changed(what)

    def load_files(
        self,
        filenames,
        filename_attributes=None,
        fileformat="mseed",
        cache=None,
        show_progress=True,
        update_progress=None,
    ):
        l = loader(
            filenames,
            fileformat,
            cache,
            filename_attributes,
            show_progress=show_progress,
            update_progress=update_progress,
        )
        self.add_files(l)

    def add_files(self, files):
        for file in files:
            self.add_file(file)

    def add_file(self, file):
        if file.abspath is not None and file.abspath in self.abspaths:
            logger.warn("File already in pile: %s" % file.abspath)
            return

        subpile = self.dispatch(file)
        subpile.add_file(file)
        if file.abspath is not None:
            self.abspaths.add(file.abspath)

    def remove_file(self, file):
        subpile = file.get_parent()
        subpile.remove_file(file)
        if file.abspath is not None:
            self.abspaths.remove(file.abspath)

    def remove_files(self, files):
        subpile_files = {}
        for file in files:
            subpile = file.get_parent()
            if subpile not in subpile_files:
                subpile_files[subpile] = []

            subpile_files[subpile].append(file)

        for subpile, files in subpile_files.iteritems():
            subpile.remove_files(files)
            for file in files:
                if file.abspath is not None:
                    self.abspaths.remove(file.abspath)

    def dispatch_key(self, file):
        dt = int(math.floor(math.log(file.deltatmin)))
        return dt

    def dispatch(self, file):
        k = self.dispatch_key(file)
        if k not in self.subpiles:
            self.subpiles[k] = SubPile(self)

        return self.subpiles[k]

    def get_deltats(self):
        return self.deltats.keys()

    def chop(
        self,
        tmin,
        tmax,
        group_selector=None,
        trace_selector=None,
        snap=(round, round),
        include_last=False,
        load_data=True,
    ):
        chopped = []
        used_files = set()

        traces = self.relevant(tmin, tmax, group_selector, trace_selector)
        if load_data:
            files_changed = False
            for tr in traces:
                if tr.file not in used_files:
                    if tr.file.load_data():
                        files_changed = True

                    used_files.add(tr.file)

            if files_changed:
                traces = self.relevant(tmin, tmax, group_selector, trace_selector)

        for tr in traces:
            try:
                chopped.append(
                    tr.chop(
                        tmin, tmax, inplace=False, snap=snap, include_last=include_last
                    )
                )
            except trace.NoData:
                pass

        return chopped, used_files

    def _process_chopped(
        self, chopped, degap, maxgap, maxlap, want_incomplete, wmax, wmin, tpad
    ):
        chopped.sort(lambda a, b: cmp(a.full_id, b.full_id))
        if degap:
            chopped = degapper(chopped, maxgap=maxgap, maxlap=maxlap)

        if not want_incomplete:
            wlen = (wmax + tpad) - (wmin - tpad)
            chopped_weeded = []
            for tr in chopped:
                emin = tr.tmin - (wmin - tpad)
                emax = tr.tmax + tr.deltat - (wmax + tpad)
                if abs(emin) <= 0.5 * tr.deltat and abs(emax) <= 0.5 * tr.deltat:
                    chopped_weeded.append(tr)
                elif degap:
                    if 0.0 < emin <= 5.0 * tr.deltat and -5.0 * tr.deltat <= emax < 0.0:
                        tr.extend(
                            wmin - tpad, wmax + tpad - tr.deltat, fillmethod="repeat"
                        )
                        chopped_weeded.append(tr)

            chopped = chopped_weeded

        for tr in chopped:
            tr.wmin = wmin
            tr.wmax = wmax

        return chopped

    def chopper(
        self,
        tmin=None,
        tmax=None,
        tinc=None,
        tpad=0.0,
        group_selector=None,
        trace_selector=None,
        want_incomplete=True,
        degap=True,
        maxgap=5,
        maxlap=None,
        keep_current_files_open=False,
        accessor_id=None,
        snap=(round, round),
        include_last=False,
        load_data=True,
    ):

        if tmin is None:
            tmin = self.tmin + tpad

        if tmax is None:
            tmax = self.tmax - tpad

        if tinc is None:
            tinc = tmax - tmin

        if not self.is_relevant(tmin - tpad, tmax + tpad, group_selector):
            return

        if accessor_id not in self.open_files:
            self.open_files[accessor_id] = set()

        open_files = self.open_files[accessor_id]

        iwin = 0
        while True:
            chopped = []
            wmin, wmax = tmin + iwin * tinc, min(tmin + (iwin + 1) * tinc, tmax)
            eps = tinc * 1e-6
            if wmin >= tmax - eps:
                break
            chopped, used_files = self.chop(
                wmin - tpad,
                wmax + tpad,
                group_selector,
                trace_selector,
                snap,
                include_last,
                load_data,
            )
            for file in used_files - open_files:
                # increment datause counter on newly opened files
                file.use_data()

            open_files.update(used_files)

            processed = self._process_chopped(
                chopped, degap, maxgap, maxlap, want_incomplete, wmax, wmin, tpad
            )
            yield processed

            unused_files = open_files - used_files

            while unused_files:
                file = unused_files.pop()
                file.drop_data()
                open_files.remove(file)

            iwin += 1

        if not keep_current_files_open:
            while open_files:
                file = open_files.pop()
                file.drop_data()

    def all(self, *args, **kwargs):
        alltraces = []
        for traces in self.chopper(*args, **kwargs):
            alltraces.extend(traces)

        return alltraces

    def iter_all(self, *args, **kwargs):
        for traces in self.chopper(*args, **kwargs):
            for trace in traces:
                yield trace

    def chopper_grouped(self, gather, progress=None, *args, **kwargs):
        keys = self.gather_keys(gather)
        if len(keys) == 0:
            return
        outer_group_selector = None
        if "group_selector" in kwargs:
            outer_group_selector = kwargs["group_selector"]

        outer_trace_selector = None
        if "trace_selector" in kwargs:
            outer_trace_selector = kwargs["trace_selector"]

        # the use of this gather-cache makes it impossible to modify the pile
        # during chopping
        gather_cache = {}
        pbar = None
        if progress is not None:
            pbar = util.progressbar(progress, len(keys))

        for ikey, key in enumerate(keys):

            def tsel(tr):
                return gather(tr) == key and (
                    outer_trace_selector is None or outer_trace_selector(tr)
                )

            def gsel(gr):
                if gr not in gather_cache:
                    gather_cache[gr] = gr.gather_keys(gather)

                return key in gather_cache[gr] and (
                    outer_group_selector is None or outer_group_selector(gr)
                )

            kwargs["trace_selector"] = tsel
            kwargs["group_selector"] = gsel

            for traces in self.chopper(*args, **kwargs):
                yield traces

            if pbar:
                pbar.update(ikey + 1)

        if pbar:
            pbar.finish()

    def gather_keys(self, gather, selector=None):
        keys = set()
        for subpile in self.subpiles.values():
            keys |= subpile.gather_keys(gather, selector)

        return sorted(keys)

    def iter_traces(
        self,
        load_data=False,
        return_abspath=False,
        group_selector=None,
        trace_selector=None,
    ):
        for subpile in self.subpiles.values():
            if not group_selector or group_selector(subpile):
                for tr in subpile.iter_traces(
                    load_data, return_abspath, group_selector, trace_selector
                ):
                    yield tr

    def iter_files(self):
        for subpile in self.subpiles.values():
            for file in subpile.iter_files():
                yield file

    def reload_modified(self):
        modified = False
        for subpile in self.subpiles.values():
            modified |= subpile.reload_modified()

        return modified

    def get_tmin(self):
        return self.tmin

    def get_tmax(self):
        return self.tmax

    def get_deltatmin(self):
        return self.deltatmin

    def get_deltatmax(self):
        return self.deltatmax

    def __str__(self):
        s = "Pile\n"
        s += "number of subpiles: %i\n" % len(self.subpiles)
        s += "timerange: %s - %s\n" % (
            util.time_to_str(self.tmin),
            util.time_to_str(self.tmax),
        )
        s += "networks: %s\n" % ", ".join(sl(self.networks.keys()))
        s += "stations: %s\n" % ", ".join(sl(self.stations.keys()))
        s += "locations: %s\n" % ", ".join(sl(self.locations.keys()))
        s += "channels: %s\n" % ", ".join(sl(self.channels.keys()))
        s += "deltats: %s\n" % ", ".join(sl(self.deltats.keys()))
        return s

    def snuffle(self, **kwargs):
        """Visualize it.

        :param stations: list of `pyrocko.model.Station` objects or ``None``
        :param events: list of `pyrocko.model.Event` objects or ``None``
        :param markers: list of `pyrocko.gui_util.Marker` objects or ``None``
        :param ntracks: float, number of tracks to be shown initially (default: 12)
        :param follow: time interval (in seconds) for real time follow mode or ``None``
        :param controls: bool, whether to show the main controls (default: ``True``)
        :param opengl: bool, whether to use opengl (default: ``False``)
        """

        from pyrocko.snuffler import snuffle

        snuffle(self, **kwargs)


def make_pile(
    paths=None,
    selector=None,
    regex=None,
    fileformat="mseed",
    cachedirname=config.cache_dir,
    show_progress=True,
):
    """Create pile from given file and directory names.

    :param paths: filenames and/or directories to look for traces. If paths is
        ``None`` ``sys.argv[1:]`` is used.
    :param selector: lambda expression taking group dict of regex match object as
        a single argument and which returns true or false to keep or reject
        a file
    :param regex: regular expression which filenames have to match
    :param fileformat: format of the files ('mseed', 'sac', 'kan',
        'from_extension', 'try')
    :param cachedirname: loader cache is stored under this directory. It is
        created as neccessary.
    :param show_progress: show progress bar and other progress information
    """
    if isinstance(paths, str):
        paths = [paths]

    if paths is None:
        paths = sys.argv[1:]

    fns = util.select_files(paths, selector, regex, show_progress=show_progress)

    cache = get_cache(cachedirname)
    p = Pile()
    p.load_files(
        sorted(fns), cache=cache, fileformat=fileformat, show_progress=show_progress
    )
    return p


class Injector(trace.States):
    def __init__(
        self,
        pile,
        fixation_length=None,
        path=None,
        format="from_extension",
        forget_fixed=False,
    ):
        trace.States.__init__(self)
        self._pile = pile
        self._fixation_length = fixation_length
        self._format = format
        self._path = path
        self._forget_fixed = forget_fixed

    def set_fixation_length(self, l):
        """Set length after which the fixation method is called on buffer traces.

        The length should be given in seconds. Give None to disable.
        """
        self.fixate_all()
        self._fixation_length = l  # in seconds

    def set_save_path(
        self,
        path="dump_%(network)s.%(station)s.%(location)s.%(channel)s_%(tmin)s_%(tmax)s.mseed",
    ):
        self.fixate_all()
        self._path = path

    def inject(self, trace):
        logger.debug("Received a trace: %s" % trace)

        buf = self.get(trace)
        if buf is None:
            trbuf = trace.copy()
            buf = MemTracesFile(None, [trbuf])
            self._pile.add_file(buf)
            self.set(trace, buf)

        else:
            self._pile.remove_file(buf)
            trbuf = buf.get_traces()[0]
            buf.remove(trbuf)
            trbuf.append(trace.ydata)
            buf.add(trbuf)
            self._pile.add_file(buf)
            self.set(trace, buf)

        trbuf = buf.get_traces()[0]
        if self._fixation_length is not None:
            if trbuf.tmax - trbuf.tmin > self._fixation_length:
                self._fixate(buf, complete=False)

    def fixate_all(self):
        for state in self._states.values():
            self._fixate(state[-1])

        self._states = {}

    def free(self, buf):
        self._fixate(buf)

    def _fixate(self, buf, complete=True):
        trbuf = buf.get_traces()[0]
        del_state = True
        if self._path:
            if self._fixation_length is not None:
                ttmin = trbuf.tmin
                ytmin = util.year_start(ttmin)
                n = int(math.floor((ttmin - ytmin) / self._fixation_length))
                tmin = ytmin + n * self._fixation_length
                traces = []
                t = tmin
                while t <= trbuf.tmax:
                    try:
                        traces.append(
                            trbuf.chop(
                                t,
                                t + self._fixation_length,
                                inplace=False,
                                snap=(math.ceil, math.ceil),
                            )
                        )
                    except trace.NoData:
                        pass
                    t += self._fixation_length

                if (
                    abs(traces[-1].tmax - (t - trbuf.deltat)) < trbuf.deltat / 100.0
                    or complete
                ):
                    self._pile.remove_file(buf)

                else:  # reinsert incomplete last part
                    new_trbuf = traces.pop()
                    self._pile.remove_file(buf)
                    buf.remove(trbuf)
                    buf.add(new_trbuf)
                    self._pile.add_file(buf)
                    del_state = False

            else:
                traces = [trbuf]
                self._pile.remove_file(buf)

            fns = io.save(traces, self._path, format=self._format)

            if not self._forget_fixed:
                self._pile.load_files(fns, show_progress=False, fileformat=self._format)

        if del_state:
            del self._states[trbuf.nslc_id]

    def __del__(self):
        self.fixate_all()
