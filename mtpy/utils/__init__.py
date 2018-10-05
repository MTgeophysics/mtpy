# Check for gdal availability at module level so we don't have to
# do this every time a function in gis_tools is being called.
from .decorator import gdal_data_check
import os, re

HAS_GDAL = gdal_data_check(None)._gdal_data_found

if (not HAS_GDAL):
    try:
        import pyproj
    except ImportError:
        raise RuntimeError("Either GDAL or PyProj must be installed")
        # end try
# end if

EPSG_DICT = {}
try:
    import pyproj

    epsgfn = os.path.join(pyproj.pyproj_datadir, 'epsg')

    f = open(epsgfn, 'r')
    lines = f.readlines()

    for line in lines:
        if ('#' in line): continue

        epsg_code_val = re.compile('<(\d+)>').findall(line)

        # print( "epsg_code_val", epsg_code_val)

        if epsg_code_val is not None and len(epsg_code_val) > 0 and epsg_code_val[0].isdigit():
            epsg_code = int(epsg_code_val[0])
            epsg_string = re.compile('>(.*)<').findall(line)[0].strip()

            EPSG_DICT[epsg_code] = epsg_string
        else:
            pass  #print("epsg_code_val NOT found for this line ", line, epsg_code_val)

except ImportError:
    pass
# end try
