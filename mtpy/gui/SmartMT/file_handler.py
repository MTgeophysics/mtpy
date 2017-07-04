"""
    Description:
        todo to be written

    Usage:
        todo tobe written

    Author: YingzhiGou
    Date: 20/06/2017
"""

import os

import mtpy.core.mt as mt

from mtpy.utils.mtpylog import MtPyLog

DEFAULT_GROUP_PREFIX = 'Group'
DEFAULT_GROUP = "Default Group"


class FileHandler:
    """
        Description:
            container that holds all file references and MT object created from the files
    """

    def __init__(self):
        self._station_dict = dict()
        self._logger = MtPyLog().get_mtpy_logger(__name__)
        self._file_dict = dict()
        self._group_dict = dict()
        self._file_to_groups = dict()

    def add_file(self, file_name, group_id=None):
        """

        :param file_name:
        :param group_id:
        :type file_name: str
        :type group_id: str
        :return:
        """
        file_ref = mt_obj = None
        if isinstance(file_name, str):
            if os.path.isfile(file_name):
                if file_name in self._file_dict and self._file_dict[file_name] is not None:
                    self._logger.warning("File %s already loaded." % file_name)
                    file_ref = file_name
                    mt_obj = self.get_MT_obj(file_name)
                else:
                    file_ref = file_name
                    self._logger.info("loading %s" % file_name)
                    mt_obj = mt.MT(file_name)
        elif isinstance(file_name, mt.MT):
            mt_obj = file_name
            file_ref = mt_obj.fn
        else:
            raise FileHandlingException("Unsupported input type %s" % type(file_name))

        # add file in to container
        self._logger.info("referencing %s to %s" % (file_ref, mt_obj.station))
        if file_ref not in self._file_dict:
            self._file_dict[file_ref] = mt_obj
            self._station_dict[mt_obj.station] = file_ref
            self._file_to_groups[file_ref] = set()
        # add file to group
        return self.add_to_group(group_id, file_ref)

    def add_files(self, file_list, group_id=None):
        for file_name in file_list:
            self.add_file(file_name, group_id)
        return True

    def station2ref(self, station):
        return self._station_dict[station]

    def add_to_group(self, group_ids, file_ref):
        """
        add a file ref to a group
        :param group_id:
        :type group_id str
        :param file_ref:
        :type file_ref str
        :return: True
        """
        if not group_ids:
            group_ids = [DEFAULT_GROUP]
        elif isinstance(group_ids, str):
            group_ids = [group_ids]
        for group_id in group_ids:
            if isinstance(group_id, str):
                if group_id not in self._group_dict:
                    self._group_dict[group_id] = set()
                if file_ref in self._file_dict:
                    self._logger.info("adding %s to group \"%s\"" % (self._file_dict[file_ref].station, group_id))
                    self._group_dict[group_id].add(file_ref)
                    self._file_to_groups[file_ref].add(group_id)
                    self._logger.info("%s now in group %s" % (self._file_dict[file_ref].station, ", ".join(self._file_to_groups[file_ref])))
                else:
                    self._logger.error("File %s has not yet been loaded." % file_ref)
                    return False
            else:
                self._logger.warning("Unsupported group ID \"%s\", add file %s to \"%s\"" % (
                    type(group_id), file_ref, DEFAULT_GROUP))
                return self.add_to_group(DEFAULT_GROUP, file_ref)
        return True

    def remove_file_from_group(self, group_id, file_ref):
        """

        :param group_id:
        :type group_id str
        :param file_ref:
        :type file_ref str
        :return:
        """
        self._logger.info("Remove %s from group %s" % (self.get_MT_obj(file_ref).station, group_id))
        if group_id in self._group_dict:
            try:
                self._group_dict[group_id].remove(file_ref)
                self._file_to_groups[file_ref].remove(group_id)
                if not self._file_to_groups[file_ref]:
                    self.add_to_group(DEFAULT_GROUP, file_ref)
                return True
            except KeyError:
                return False
        return False

    def unload(self, file_ref):
        self._logger.info("unload %s" % file_ref)
        for group in self._file_to_groups[file_ref]:
            self._group_dict[group].remove(file_ref)
        station = self.get_MT_obj(file_ref).station
        del self._station_dict[station]
        del self._file_to_groups[file_ref]
        del self._file_dict[file_ref]


    def remove_group(self, group_id):
        self._logger.info("Remove group %s" % group_id)
        members = self.get_group_members(group_id)
        if members is not None:
            for ref in members:
                self._file_to_groups[ref].remove(group_id)
                if not self._file_to_groups[ref]:
                    self.add_to_group(DEFAULT_GROUP, ref)
            del self._group_dict[group_id]

    def get_groups(self):
        return self._group_dict.keys()
    # properties

    def get_group_members(self, group):
        if group in self._group_dict:
            return self._group_dict[group]
        else:
            self._logger.error("Group \"%s\" does not exist." % group)
            return None

    def get_MT_obj(self, ref):
        if ref in self._file_dict:
            return self._file_dict[ref]
        else:
            self._logger.error("File \"%s\" is not loaded" % ref)
            return None

    def get_file_refs(self):
        return self._file_dict.keys()

    def create_group(self, group_id):
        if group_id not in self._group_dict:
            self._group_dict[group_id] = set()
            return True
        else:
            self._logger.warning("Group %s exists!" % group_id)
            return False


class FileHandlingException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
