
#from shutil import which
import os
import shutil


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which

    return shutil.which(name) is not None


def check_output_dir(
    dir,
    overwrite = False):

    #check if output folder exists, and whether to overwrite
    if os.path.exists(dir) and overwrite == False:
        raise ValueError("{dir} exists and overwrite not set".format(dir = dir))
    elif os.path.exists(dir) and overwrite == True:
        shutil.rmtree(dir)
        os.makedirs(dir)
    else: # make output directory
        os.makedirs(dir)
