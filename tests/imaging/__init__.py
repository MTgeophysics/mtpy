from __future__ import print_function

import functools
import inspect
import os
import sys
import threading
from unittest import TestCase

import matplotlib

from mtpy.utils.mtpylog import MtPyLog
from tests import TEST_DIR, make_temp_dir, TEST_TEMP_DIR

if os.name == "posix" and "DISPLAY" not in os.environ:
    print(
        "MATPLOTLIB: No Display found, using non-interactive svg backend",
        file=sys.stderr,
    )
    matplotlib.use("svg")
    import matplotlib.pyplot as plt

    MTPY_TEST_HAS_DISPLAY = False
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt

    MTPY_TEST_HAS_DISPLAY = True
    plt.ion()

MtPyLog.get_mtpy_logger(__name__).info(
    "Testing using matplotlib backend {}".format(matplotlib.rcParams["backend"])
)


def reset_matplotlib():
    # save some important params
    interactive = matplotlib.rcParams["interactive"]
    backend = matplotlib.rcParams["backend"]
    # reset
    matplotlib.rcdefaults()  # reset the rcparams to default
    # recover
    matplotlib.rcParams["backend"] = backend
    matplotlib.rcParams["interactive"] = interactive
    logger = MtPyLog().get_mtpy_logger(__name__)
    logger.info(
        "Testing using matplotlib backend {}".format(matplotlib.rcParams["backend"])
    )


class ImageCompare(object):
    """
    - to enable the image comparison tests in the code, add @ImageCompare(**kwargs)
    to the function that generates image with matplotlib, by default, the new image
    will be saved and check if the image is an empty image. to
    - enable the image comparison tests against the baseline image stored, set environment variable
    MTPY_TEST_COMPARE_IMAGE = True
    """

    def __init__(self, *args, **kwargs):
        self.baseline_dir = kwargs.pop(
            "baseline_dir",
            # 'tests/baseline_images/matplotlib_{ver}'.format(ver=matplotlib.__version__).replace('.', '_')
            os.path.normpath(os.path.join(TEST_DIR, "baseline_images")),
        )
        self.result_dir = kwargs.pop(
            "result_dir",
            # 'tests/result_images/matplotlib_{ver}'.format(ver=matplotlib.__version__).replace('.', '_')
            os.path.normpath(os.path.join(TEST_TEMP_DIR, "image_compare_tests")),
        )
        self.filename = kwargs.pop("filename", None)
        self.extensions = kwargs.pop("extensions", ["png"])
        self.savefig_kwargs = kwargs.pop("savefig_kwargs", {"dpi": 80})
        self.tolerance = kwargs.pop("tolerance", 2)
        self.fig_size = kwargs.pop("fig_size", None)
        self.is_compare_image = self.to_bool(
            os.getenv("MTPY_TEST_COMPARE_IMAGE", False)
        )
        self._logger = MtPyLog().get_mtpy_logger(__name__)
        self._logger.info(
            "Image Comparison Test: {stat}".format(
                stat="ENABLED" if self.is_compare_image else "DISABLED"
            )
        )
        self.on_fail = kwargs.pop("on_fail", None)
        self.on_compare_fail = kwargs.pop("on_compare_fail", None)
        self.on_empty_image = kwargs.pop("on_empty_image", None)

        ImageCompare._thread_lock.acquire()
        if not os.path.exists(self.result_dir):
            os.mkdir(self.result_dir)
        ImageCompare._thread_lock.release()

    def __call__(self, original):
        import matplotlib.pyplot as plt
        from matplotlib.testing.compare import compare_images

        if self.filename is None:
            filename = original.__name__
        else:
            filename = self.filename

        filename = filename.replace("[", "_").replace("]", "_").replace("/", "_")
        filename = filename.strip(" _")

        test_suite_name = original.__module__.split(".")[-1]

        @functools.wraps(original)
        def new_test_func(*args, **kwargs):
            if inspect.ismethod(original):
                result = original.__func__(*args, **kwargs)
            else:
                result = original(*args, **kwargs)
            for (
                baseline_image,
                test_image,
                baseline_rcparams,
                test_rcparams,
            ) in self._get_baseline_result_pairs(
                test_suite_name, filename, self.extensions
            ):
                # save image
                fig = plt.gcf()
                if fig is not None:
                    if self.fig_size is not None:
                        fig.set_size_inches(self.fig_size)
                        fig.set_tight_layout(True)
                    fig.savefig(test_image, **self.savefig_kwargs)
                    import pytest

                    if self.is_compare_image and os.path.exists(baseline_image):
                        msg = compare_images(
                            baseline_image, test_image, tol=self.tolerance
                        )
                        if msg is not None:
                            msg += "\n"
                            msg += self.compare_rcParam(
                                baseline_rcparams, test_rcparams
                            )
                            # print image in base64
                            # print("====================")
                            # print("Expected Image:")
                            # self._print_image_base64(baseline_image)
                            # print("Actual Image:")
                            # self._print_image_base64(test_image)
                            # print("====================")
                            self.print_image_testing_note(file=sys.stderr)
                            if self.on_compare_fail is not None:
                                self.on_compare_fail()
                            if self.on_fail is not None:
                                self.on_fail()
                            pytest.fail(msg, pytrace=False)
                        else:
                            # clearup the image as they are the same with the baseline
                            os.remove(test_image)
                            os.remove(test_rcparams)
                            if not os.listdir(os.path.dirname(test_image)):
                                os.rmdir(os.path.dirname(test_image))
                    else:
                        # checking if the created image is empty
                        # verify(test_image) # if issues with test_image, nonexistent? raise exception.
                        actual_image = matplotlib.pyplot.imread(test_image)
                        actual_image = actual_image[
                            :, :, :3
                        ]  # remove the alpha channel (if exists)
                        import numpy as np

                        if np.any(actual_image):
                            self.print_image_testing_note(file=sys.stderr)
                            if self.is_compare_image:
                                pytest.skip(
                                    "Image file not found for comparison test "
                                    "(This is expected for new tests.)\nGenerated Image: "
                                    "\n\t{test}".format(test=test_image)
                                )
                            else:
                                self._logger.info(
                                    "\nGenerated Image: {test}".format(test=test_image)
                                )
                        else:
                            # empty image created
                            if self.on_empty_image is not None:
                                self.on_empty_image()
                            if self.on_fail is not None:
                                self.on_fail()
                            pytest.fail(
                                "Image file not found for comparison test "
                                "(This is expected for new tests.),"
                                " but the new image created is empty."
                            )
            return result

        return new_test_func

    def _get_baseline_result_pairs(self, test_suite_name, fname, extensions):
        if test_suite_name is None:
            baseline = self.baseline_dir
            result = self.result_dir
        else:
            baseline = os.path.join(self.baseline_dir, test_suite_name)
            result = os.path.join(self.result_dir, test_suite_name)
            if not os.path.exists(baseline):
                os.makedirs(baseline)
            if not os.path.exists(result):
                os.makedirs(result)

        for ext in extensions:
            name = "{fname}.{ext}".format(fname=fname, ext=ext)
            rc_name = "{fname}_rcParams.txt".format(fname=fname)
            yield os.path.normpath(os.path.join(baseline, name)), os.path.normpath(
                os.path.join(result, name)
            ), os.path.normpath(os.path.join(baseline, rc_name)), os.path.normpath(
                os.path.join(result, rc_name)
            )

    @staticmethod
    def _print_image_base64(image_file_name):
        with open(image_file_name, "rb") as image_file:
            image_data = image_file.read()
            print(
                '<img src="data:image/{};base64,{}" style="display:block; max-width:800px; width: auto; height: '
                'auto;" />'.format(
                    os.path.splitext(image_file_name)[1].strip(" ."),
                    image_data.encode("base64"),
                )
            )

    @staticmethod
    def print_image_testing_note(file=sys.stdout):
        print("====================", file=file)
        print("matplotlib Version: " + matplotlib.__version__, file=file)
        print(
            "NOTE: The test result may be different in different versions of matplotlib.",
            file=file,
        )
        print("====================", file=file)

    @staticmethod
    def compare_rcParam(baseline_rcparams, test_rcparams):
        # check if the rcParams are different
        msg = "Comparing matplotlib.rcParam:\n"
        if os.path.isfile(baseline_rcparams):
            with open(baseline_rcparams, "r") as fbaserc:
                with open(test_rcparams, "r") as ftestrc:
                    import difflib

                    lines = [
                        line
                        for line in difflib.unified_diff(
                            fbaserc.readlines(),
                            ftestrc.readlines(),
                            fromfile="baseline",
                            tofile="test",
                            n=0,
                        )
                    ]
                    if lines:
                        msg += "  Found differences:\n    " + "    ".join(lines)
                    else:
                        msg += " NO differences found."
        else:
            msg += "  Baseline rcParams file not found."
        return msg

    def to_bool(self, param):
        if isinstance(param, str) and param:
            param = param.lower()
            if param in ("true", "1", "t"):
                return True
            elif param in ("false", "f", "0"):
                return False
        else:
            return bool(param)


ImageCompare._thread_lock = threading.Lock()
ImageCompare.print_image_testing_note(file=sys.stderr)

_thread_lock = threading.Lock()


class ImageTestCase(TestCase):
    @classmethod
    def setUpClass(cls):
        _thread_lock.acquire()
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__.split(".")[-1])

    @classmethod
    def tearDownClass(cls):
        plt_close("all")
        _thread_lock.release()

    def setUp(self):
        if plt.get_fignums():
            plt.clf()
        reset_matplotlib()
        if plt.isinteractive():
            plt.show(block=False)  # show an empty window first for drawing

    def tearDown(self):
        if MTPY_TEST_HAS_DISPLAY:
            plt_wait(1)
        plt_close()


def plt_wait(seconds):
    if plt.isinteractive() and plt.get_fignums():
        try:
            plt.pause(seconds)
        except:
            pass


def plt_close(to_close=None):
    if plt.get_fignums():
        if to_close is None:
            plt.close()
        else:
            plt.close(to_close)
