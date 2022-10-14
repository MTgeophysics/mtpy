"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

from mtpy.utils import exceptions as mtex

# =============================================================================


class ControlInv(object):
    """
    read and write control file for how the inversion starts and how it is run

    """

    def __init__(self, **kwargs):

        self.output_fn = "MODULAR_NLCG"
        self.lambda_initial = 10
        self.lambda_step = 10
        self.model_search_step = 1
        self.rms_reset_search = 2.0e-3
        self.rms_target = 1.05
        self.lambda_exit = 1.0e-4
        self.max_iterations = 100
        self.save_path = Path().cwd()
        self.fn_basename = "control.inv"

        self._control_keys = [
            "Model and data output file name",
            "Initial damping factor lambda",
            "To update lambda divide by",
            "Initial search step in model units",
            "Restart when rms diff is less than",
            "Exit search when rms is less than",
            "Exit when lambda is less than",
            "Maximum number of iterations",
        ]

        self._control_dict = dict(
            [
                (key, value)
                for key, value in zip(
                    self._control_keys,
                    [
                        self.output_fn,
                        self.lambda_initial,
                        self.lambda_step,
                        self.model_search_step,
                        self.rms_reset_search,
                        self.rms_target,
                        self.lambda_exit,
                        self.max_iterations,
                    ],
                )
            ]
        )
        self._string_fmt_dict = dict(
            [
                (key, value)
                for key, value in zip(
                    self._control_keys,
                    [
                        "<",
                        "<.1f",
                        "<.1f",
                        "<.1f",
                        "<.1e",
                        "<.2f",
                        "<.1e",
                        "<.0f",
                    ],
                )
            ]
        )

        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def control_fn(self):
        return self.save_path.joinpath(self.fn_basename)

    @control_fn.setter
    def control_fn(self, value):
        if value is not None:
            value = Path(value)
            self.save_path = value.parent
            self.fn_basename = value.name

    def write_control_file(
        self, control_fn=None, save_path=None, fn_basename=None
    ):
        """
        write control file

        Arguments:
        ------------
            **control_fn** : string
                             full path to save control file to
                             *default* is save_path/fn_basename

            **save_path** : string
                            directory path to save control file to
                            *default* is cwd

            **fn_basename** : string
                              basename of control file
                              *default* is control.inv

        """

        if control_fn is not None:
            self.control_fn = control_fn

        if save_path is not None:
            self.save_path = Path(save_path)

        if fn_basename is not None:
            self.fn_basename = fn_basename

        self._control_dict = dict(
            [
                (key, value)
                for key, value in zip(
                    self._control_keys,
                    [
                        self.output_fn,
                        self.lambda_initial,
                        self.lambda_step,
                        self.model_search_step,
                        self.rms_reset_search,
                        self.rms_target,
                        self.lambda_exit,
                        self.max_iterations,
                    ],
                )
            ]
        )

        clines = []
        for key in self._control_keys:
            value = self._control_dict[key]
            str_fmt = self._string_fmt_dict[key]
            clines.append("{0:<35}: {1:{2}}\n".format(key, value, str_fmt))

        with open(self.control_fn, "w") as cfid:
            cfid.writelines(clines)

        print("Wrote ModEM control file to {0}".format(self.control_fn))

    def read_control_file(self, control_fn=None):
        """
        read in a control file
        """

        if control_fn is not None:
            self.control_fn = control_fn

        if self.control_fn is None:
            raise mtex.MTpyError_file_handling(
                "control_fn is None, input " "control file"
            )

        if not self.control_fn.is_file():
            raise mtex.MTpyError_file_handling(
                "Could not find {0}".format(self.control_fn)
            )

        with open(self.control_fn, "r") as cfid:
            clines = cfid.readlines()
            for cline in clines:
                clist = cline.strip().split(":")
                if len(clist) == 2:

                    try:
                        self._control_dict[clist[0].strip()] = float(clist[1])
                    except ValueError:
                        self._control_dict[clist[0].strip()] = clist[1]

        # set attributes
        attr_list = [
            "output_fn",
            "lambda_initial",
            "lambda_step",
            "model_search_step",
            "rms_reset_search",
            "rms_target",
            "lambda_exit",
            "max_iterations",
        ]
        for key, kattr in zip(self._control_keys, attr_list):
            setattr(self, kattr, self._control_dict[key])
