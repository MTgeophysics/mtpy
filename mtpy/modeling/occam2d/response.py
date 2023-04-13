# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 19:09:35 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path
import numpy as np

# =============================================================================
class Response(object):
    """
    Reads .resp file output by Occam.  Similar structure to Data.data.

    If resp_fn is given in the initialization of Response, read_response_file
    is called.

    Arguments:
    ------------
        **resp_fn** : string
                      full path to .resp file

    Attributes:
    -------------
        **resp** : is a list of dictioinaries containing the data for each
                   station.  keys include:

                   * 'te_res' -- TE resisitivity in linear scale
                   * 'tm_res' -- TM resistivity in linear scale
                   * 'te_phase' -- TE phase in degrees
                   * 'tm_phase' --  TM phase in degrees in first quadrant
                   * 're_tip' -- real part of tipper along profile
                   * 'im_tip' -- imaginary part of tipper along profile

               each key is a np.ndarray(2, num_freq)
               index 0 is for model response
               index 1 is for normalized misfit

    :Example: ::
        >>> resp_obj = occam2d.Response(r"/home/occam/line1/inv1/test_01.resp")



    """

    def __init__(self, resp_fn=None, **kwargs):
        self.resp_fn = resp_fn

        self.resp = None
        self.occam_dict = {
            "1": "log_te_res",
            "2": "te_phase",
            "3": "re_tip",
            "4": "im_tip",
            "5": "log_tm_res",
            "6": "tm_phase",
            "9": "te_res",
            "10": "tm_res",
        }

        if resp_fn is not None:
            self.read_response_file()

    def read_response_file(self, resp_fn=None):
        """
        read in response file and put into a list of dictionaries similar
        to Data
        """

        if resp_fn is not None:
            self.resp_fn = resp_fn

        if self.resp_fn is None:
            raise ValueError("resp_fn is None, please input response file")

        if not self.resp_fn.is_file():
            raise ValueError(f"Could not find {self.resp_fn}")

        r_arr_tmp = np.loadtxt(self.resp_fn)
        r_arr = np.zeros(
            len(r_arr_tmp),
            dtype=[
                ("station", np.int),
                ("freq", np.int),
                ("comp", np.int),
                ("z", np.int),
                ("data", np.float),
                ("resp", np.float),
                ("err", np.float),
            ],
        )
        r_arr["station"] = r_arr_tmp[:, 0]
        r_arr["freq"] = r_arr_tmp[:, 1]
        r_arr["comp"] = r_arr_tmp[:, 2]
        r_arr["z"] = r_arr_tmp[:, 3]
        r_arr["data"] = r_arr_tmp[:, 4]
        r_arr["resp"] = r_arr_tmp[:, 5]
        r_arr["err"] = r_arr_tmp[:, 6]

        num_stat = r_arr["station"].max()
        num_freq = r_arr["freq"].max()

        # set zero array size the first row will be the data and second the
        # error
        asize = (2, num_freq)

        # make a list of dictionaries for each station.
        self.resp = [
            {
                "te_phase": np.zeros(asize),
                "tm_phase": np.zeros(asize),
                "re_tip": np.zeros(asize),
                "im_tip": np.zeros(asize),
                "te_res": np.zeros(asize),
                "tm_res": np.zeros(asize),
            }
            for ss in range(num_stat)
        ]

        for line in r_arr:
            # station index -1 cause python starts at 0
            ss = line["station"] - 1

            # frequency index -1 cause python starts at 0
            ff = line["freq"] - 1
            # data key
            key = self.occam_dict[str(line["comp"])]
            # put into array
            if line["comp"] == 1 or line["comp"] == 5:
                self.resp[ss][key[4:]][0, ff] = 10 ** line["resp"]
                # error
                self.resp[ss][key[4:]][1, ff] = line["err"] * np.log(10)
            else:
                self.resp[ss][key][0, ff] = line["resp"]
                # error
                self.resp[ss][key][1, ff] = line["err"]
