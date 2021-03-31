#! /usr/bin/env python
"""
Description:
    Example template python script structure.
    .......
    .......
    
References: 
    https://gajira.atlassian.net/browse/ALAMP-49

CreationDate:   23/03/2018
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     23/03/2018   FZ
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

# import section
import os
import click
import mtpy.imaging.penetration_depth1d as pd1d


# Section to define functions or class
def fun1():
    """
    define my function1
    :return:
    """
    print("begin fun1")

    return


def main():
    """
    define my main function
    :return:
    """
    print("Template main()")

    return


###############################################################################
# Following is code for click making inputs to the plot depth
###############################################################################


@click.command()
@click.option(
    "-i",
    "--input",
    type=str,
    default="examples/data/edi_files",
    help="directory or edsi data files",
)
@click.option(
    "-o", "--output_file", type=str, default="temp", help="save jpg image file"
)
def plot_penetration_depth(input, output_file):
    if os.path.isfile(input):
        pd1d.plot_edi_file(input, savefile=output_file)
    elif os.path.isdir(input):
        pd1d.plot_edi_dir(input, rholist=["det"])
    else:
        pass


# =============================================
# Section for quick test of this script
# ---------------------------------------------

if __name__ == "__main__":
    plot_penetration_depth()
