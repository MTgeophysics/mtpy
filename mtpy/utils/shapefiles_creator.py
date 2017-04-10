"""
create shape files for MT datasets.
Phase Tensor, Tipper Real/Imag, MT-site locations,etc

fei.zhang@ga.gov.au
2017-03-06
"""
import os
import sys
import glob

import numpy as np
import geopandas

import mtpy.analysis.pt as mtpt
import mtpy.core.mt as mt
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps

class ShapeFilesCreator:

    def __init__(self, edifile_list, outdir):
        """
        loop through a list of edi files, create required shapefiles
        :param edifile_list: [path2edi,...]
        :param outdir: path2output dir, where the shpe file weill be written.
        """

        self.edifiles = edifile_list
        print ("number of edi files to be processed:", len(self.edifiles))
        assert len(self.edifiles) > 0

        self.outputdir = outdir

        if self.edifiles is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edifiles]


        # get all frequencies from all edi files
        all_freqs=[]
        for mt_obj in self.mt_obj_list:
            all_freqs.extend(list(mt_obj.Z.freq))

        # sort all frequencies so that they are in descending order,
        # use set to remove repeats and make an array
        self.all_frequencies = sorted(list(set(all_freqs)))

        self.all_periods=1. / np.array(sorted(list(set(all_freqs)),
                                              reverse=True))
        print self.all_periods

        print ("Number of MT Periods:", len(self.all_periods))

        return


    def create_csv(self):
        """
        generate a csv db for MT sites, phase tensor attributes, etc, from edi files
        :return:
        """

        for per in self.all_frequencies:
            ptm=PlotPhaseTensorMaps(fn_list=self.edifiles, plot_freq=per)

            ptm.export_params_to_file(save_path=self.outputdir)

        pass


    def create_phase_tensor_shp(self):

        pass

    def create_tipper_shp(self):

        pass

    def create_MTsites_shp(self):

        pass


# ==================================================================
if __name__ == "__main__":

    edidir = sys.argv[1]

    edifiles = glob.glob(os.path.join(edidir, "*.edi"))

    path2out = sys.argv[2]

    shp_maker = ShapeFilesCreator(edifiles, path2out)

    shp_maker.create_csv()

    #shp_maker.create_MTsites_shp()
