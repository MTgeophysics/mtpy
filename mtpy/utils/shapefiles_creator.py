"""
create shape files for MT datasets:
Phase Tensor, Tipper Real/Imag, MT-site locations,etc

fei.zhang@ga.gov.au
2017-03-06
"""
import os, sys
import glob
import geopandas

class ShapeFilesCreator:

    def __init__(self, edifile_list, outdir):
        """
        loop through a list of edi files, create required shapefiles
        :param edifile_list: [path2edi,...]
        :param outdir: path2output dir, where the shpe file weill be written.
        """

        self.edifiles= edifile_list
        self.outputdir=outdir

        print ("number of edi files to be processed:", len(self.edifiles))

        assert len(self.edifiles) > 0


    def create_phase_tensor_shp(self):

        pass

    def create_tipper_shp(self):

        pass


    def create_MTsites_shp(self):

        pass


# ==================================================================
if __name__=="__main__":

    edidir= sys.argv[1]

    edifiles=glob.glob(os.path.join(edidir,"*.edi"))

    path2out=sys.argv[2]

    shp_maker= ShapeFilesCreator(edifiles, path2out)

    shp_maker.create_MTsites_shp()