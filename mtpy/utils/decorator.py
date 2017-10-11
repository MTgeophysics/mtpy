"""
    Description:
        used to mark functions, methods and classes deprecated, and prints warning message when it called
        decorators based on https://stackoverflow.com/a/40301488

    Usage:
        todo: write usage

    Author: YingzhiGou
    Date: 20/06/2017
"""
from __future__ import print_function
import functools
import inspect
import os
from matplotlib import _png

import matplotlib
import sys

from matplotlib.testing.compare import verify


class deprecated(object):
    def __init__(self, reason):
        if inspect.isclass(reason) or inspect.isfunction(reason):
            raise TypeError("Reason for deprecation must be supplied")
        self.reason = reason

    def __call__(self, cls_or_func):
        if inspect.isfunction(cls_or_func):
            if hasattr(cls_or_func, 'func_code'):
                _code = cls_or_func.func_code
            else:
                _code = cls_or_func.__code__
            fmt = "Call to deprecated function or method {name} ({reason})."
            filename = _code.co_filename
            lineno = _code.co_firstlineno + 1

        elif inspect.isclass(cls_or_func):
            fmt = "Call to deprecated class {name} ({reason})."
            filename = cls_or_func.__module__
            lineno = 1

        else:
            raise TypeError(type(cls_or_func))

        msg = fmt.format(name=cls_or_func.__name__, reason=self.reason)

        @functools.wraps(cls_or_func)
        def new_func(*args, **kwargs):
            import warnings
            warnings.simplefilter('always', DeprecationWarning)  # turn off filter
            warnings.warn_explicit(msg, category=DeprecationWarning, filename=filename, lineno=lineno)
            warnings.simplefilter('default', DeprecationWarning)  # reset filter
            return cls_or_func(*args, **kwargs)

        return new_func


class ImageCompare(object):
    def __init__(self, *args, **kwargs):
        self.baseline_dir = kwargs.pop(
            'baseline_dir',
            'tests/baseline_images/matplotlib_{ver}'.format(ver=matplotlib.__version__).replace('.', '_'))
        self.result_dir = kwargs.pop(
            'result_dir',
            'tests/result_images/matplotlib_{ver}'.format(ver=matplotlib.__version__).replace('.', '_'))
        self.filename = kwargs.pop('filename', None)
        self.extensions = kwargs.pop('extensions', ['png'])
        self.savefig_kwargs = kwargs.pop('savefig_kwargs', {'dpi': 80})
        self.tolerance = kwargs.pop('tolerance', 2)
        self.fig_size = kwargs.pop('fig_size', None)

        if self.result_dir and not os.path.exists(self.result_dir):
            os.makedirs(self.result_dir)
        if self.baseline_dir and not os.path.exists(self.baseline_dir):
            os.makedirs(self.baseline_dir)

    def __call__(self, original):
        import matplotlib.pyplot as plt
        from matplotlib.testing.compare import compare_images

        if self.filename is None:
            filename = original.__name__
        else:
            filename = self.filename

        filename = filename.replace('[', '_').replace(']', '_').replace('/', '_')
        filename = filename.strip(' _')

        test_suite_name = original.__module__.split('.')[-1]

        @functools.wraps(original)
        def new_test_func(*args, **kwargs):
            if inspect.ismethod(original):
                result = original.__func__(*args, **kwargs)
            else:
                result = original(*args, **kwargs)
            for baseline_image, test_image in self._get_baseline_result_pairs(
                    test_suite_name,
                    filename,
                    self.extensions
            ):
                # save image
                fig = plt.gcf()
                if fig is not None:
                    if self.fig_size is not None:
                        fig.set_size_inches(self.fig_size)
                        fig.set_tight_layout(True)
                    fig.savefig(test_image, **self.savefig_kwargs)
                    import pytest
                    if os.path.exists(baseline_image):
                        msg = compare_images(baseline_image, test_image, tol=self.tolerance)
                        if msg is not None:
                            # print image in base64
                            # print("====================")
                            # print("Expected Image:")
                            # self._print_image_base64(baseline_image)
                            # print("Actual Image:")
                            # self._print_image_base64(test_image)
                            # print("====================")
                            self.print_image_testing_note(file=sys.stderr)
                            pytest.fail(msg, pytrace=False)
                        else:
                            # clearup the image as they are the same with the baseline
                            os.remove(test_image)
                            if not os.listdir(os.path.dirname(test_image)):
                                os.rmdir(os.path.dirname(test_image))
                    else:
                        # checking if the created image is empty
                        verify(test_image)
                        actual_image = _png.read_png_int(test_image)
                        actual_image = actual_image[:, :, :3]  # remove the alpha channel (if exists)
                        import numpy as np
                        if np.any(actual_image):
                            self.print_image_testing_note(file=sys.stderr)
                            pytest.xfail("Image file not found for comparison test "
                                        "(This is expected for new tests.)\nGenerated Image: "
                                        "\n\t{test}".format(test=test_image))
                        else:
                            # empty image created
                            pytest.fail("Image file not found for comparison test "
                                        "(This is expected for new tests.),"
                                        " but the new image created is empty.")
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
            name = '{fname}.{ext}'.format(fname=fname, ext=ext)
            yield os.path.normpath(
                os.path.join(
                    baseline,
                    name
                )
            ), os.path.normpath(
                os.path.join(
                    result,
                    name
                )
            )

    def _print_image_base64(self, image_file_name):
        with open(image_file_name, "rb") as image_file:
            image_data = image_file.read()
            print(
                "<img src=\"data:image/{};base64,{}\" style=\"display:block; max-width:800px; width: auto; height: auto;\" />".format(
                    os.path.splitext(image_file_name)[1].strip(" ."),
                    image_data.encode("base64")))

    @staticmethod
    def print_image_testing_note(file=sys.stdout):
        print("====================", file=file)
        print("matplotlib Version: " + matplotlib.__version__, file=file)
        print("NOTE: The test result may be different in different versions of matplotlib.", file=file)
        print("====================", file=file)
