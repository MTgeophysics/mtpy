"""
Description:
    To compute anc encaptutlate the properties of EDI file collection

Author: fei.zhang@ga.gov.au

Date: 2017-04-20
"""

from __future__ import print_function
import sys
import os, glob
import logging

import numpy as np

import mtpy.core.mt as mt
from mtpy.utils.mtpylog import MtPyLog

logger = MtPyLog().get_mtpy_logger(__name__)
logger.setLevel(logging.DEBUG)


class EdiCollection(object):
    """
    A super class to encapsulate the properties pertinent to a set of EDI files
    """

    def __init__(self, edilist, ptol=0.05):
        """
        constructor
        :param edilist: a list of edifiles with full path, for read-only
        :param ptol: period tolerance considered as equal, default 0.05 means 5%
        this param controls what freqs/periods are grouped together:
        10% may result more double counting of freq/period data than 5%.
        eg: E:\Data\MT_Datasets\WenPingJiang_EDI 18528 rows vs 14654 rows
        """
        self.edifiles = edilist
        logger.info("number of edi files in this collection: %s",
                    len(self.edifiles))
        assert len(self.edifiles) > 0

        self.ptol = ptol

        if self.edifiles is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edifiles]
        else:
            logger.error("None Edi file set")

        # get all frequencies from all edi files
        self.all_frequencies = None
        self.all_periods = self._get_all_periods()

        return

        self.all_periods = self._get_all_periods()
        self.bounding_box = self._get_bounding_box()

    def _get_all_periods(self):
        """
        from the list of edi files get a list of all unique frequencies.
        """
        if self.all_frequencies is not None:  # already initialized
            return

        # get all frequencies from all edi files
        all_freqs = []
        for mt_obj in self.mt_obj_list:
            all_freqs.extend(list(mt_obj.Z.freq))

        # sort all frequencies so that they are in ascending order,
        # use set to remove repeats and make an array
        self.all_frequencies = sorted(list(set(all_freqs)))

        logger.info("Number of MT Frequencies: %s", len(self.all_frequencies))

        all_periods = 1.0 / np.array(sorted(self.all_frequencies, reverse=True))

        logger.debug("Type of all_periods %s", type(all_periods))
        logger.info("Number of MT Periods: %s", len(all_periods))
        logger.debug(all_periods)

        return all_periods

    def _get_bounding_box(self, epsgcode=4326):
        """

        :return: bounding box in given proj coord system
        """
        bdict={"MinLat" : -30,"MaxLat" : -20, "MinLon": 140.1, "MaxLon" : 145.5}

        return bdict

    def show_prop(self):
        """
        show all properties
        :return:
        """
        print (len(self.all_periods), 'unique periods (s)', self.all_periods)
        print (len(self.all_frequencies), 'unique frequencies (Hz)', self.all_frequencies)

        return

if __name__ == "__main__":

    if len(sys.argv)<2:
        print ("USAGE: %s edi_dir OR edi_list " % sys.argv[0])
        sys.exit(1)
    else:
        argv1 = sys.argv[1]
        if os.path.isdir(argv1):
            edilist = glob.glob(argv1+'/*.edi')
            obj=EdiCollection(edilist)
        elif os.path.isfile(argv1) and argv1.endswith('.edi'):
            obj=EdiCollection(sys.argv[1:])  # assume input is a list of EDI files
        else:
            pass

        obj.show_prop()
