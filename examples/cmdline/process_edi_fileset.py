import os, sys
import glob
import click
from mtpy.core.edi_collection import EdiCollection

##################################################################
#
# python mtpy\core\edi_collection.py --input=examples/data/edi_files
# or
# python mtpy\core\edi_collection.py --input=
# "examples/data/edi_files/pb23c.edi examples/data/edi_files/pb23c.edi"
#
##################################################################


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-i",
    "--input",
    type=str,
    default="examples/data/edi_files",
    help="input directory to edi files or string of edi files separated by space"
    + "\n\n"
    + "python mtpy/core/edi_collection.py --input=examples/data/edi_files"
    + "\n"
    + "-or-"
    + "\n"
    + "python mtpy/core/edi_collection.py --input="
    + "\n"
    + '"examples/data/edi_files/pb23c.edi examples/data/edi_files/pb25c.edi"'
    + "\n",
)
def process_edi_files(input):
    print("Directory for edi files or single file  ---------> {}".format(input))
    edis = []
    if not (" " in input):
        if not os.path.isdir(input):
            print("Invalid Ditectory Input")
            sys.exit()
        if os.path.isdir(input):
            edis = glob.glob(input + "/*.edi")

            print(edis)

            if len(edis) == 0:
                print("Directory edi files {} empty".format(input))
                sys.exit()
            obj = EdiCollection(edilist=edis)
    else:
        edis = input.split(" ")
        for fl in edis:
            if not os.path.isfile(fl):
                print("Invalid Input File {}".format(fl))
                sys.exit()
        obj = EdiCollection(edilist=edis)
    # Compute distances
    min_dist, max_dist = obj.get_min_max_distance(obj)
    #     mt_distances = obj.get_stations_distances_stats()
    #     min_dist = mt_distances.get("MIN_DIST")
    #     max_dist = mt_distances.get("MAX_DIST")
    print("Min Distance = {}".format(min_dist))
    print("Max Distance = {}".format(max_dist))


if __name__ == "__main__":
    process_edi_files()
