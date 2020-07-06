"""
Utilities for converting EDI data to MARE2DEM compatible format for use
with MARE2DEM software.

Note in this module comments referring to "original script" mean the
EDI2Mare2DEM_withOccam2D_new.py script.

02-07-2020
brenainn.moushall@.ga.gov.au
"""
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy

import mtpy.modeling.occam2d as o2d
from mtpy.utils import mesh_tools, gis_tools, filehandling


def station_list(edi_dir):
    """
    Get list of station names from a directory containing EDI files.
    """
    return [os.path.splitext(os.path.basename(fn))[0]
            for fn in os.listdir(edi_dir) if fn.endswith('.edi')]


def line_length(x0, y0, x1, y1):
    return np.sqrt(abs(y1 - y0) ** 2 + abs(x1 - x0) ** 2)


def points_o2d_to_m2d(eastings, northings, profile_length=None):
    """
    Converts Occam2D points to Mare2D system. This is assuming
    O2D profile origin is start of line and Mare2D profile origin is
    middle of line.
    """
    converted_points = []
    if profile_length is None:
        profile_length = line_length(eastings[0], northings[0], eastings[-1], northings[-1])
    mid = profile_length / 2
    for e, n in zip(eastings, northings):
        point = line_length(e, n, eastings[0], northings[0])
        if point < mid:
            converted_points.append(-1 * (mid - point))
        elif point == mid:
            converted_points.append(0)
        elif point > mid:
            converted_points.append(point - mid)
    return np.array(converted_points)


def plot(m2d_profile, profile_elevation, site_locations, site_elevations, figsize=None):
    """
    Generate line plot of Mare2DEM profile and site locations against
    elevation.
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(m2d_profile, profile_elevation, 'b', label='M2D Profile')
    ax.plot(site_locations, site_elevations, color='r', marker='*', linestyle='None', label='Sites')
    ax.invert_yaxis()
    ax.set_xlabel("Position (Mare2D coordinates)")
    ax.set_ylabel("Elevation (meters)")
    ax.legend()
    fig.suptitle("Mare2D profile with elevation", fontsize=12)
    return fig


def occam2d_to_mare2dem(o2d_data, surface_file, elevation_sample_n=300,
                        flip_elevation=True):
    """
    Converts Occam2D profile to Mare2D format, giving station locations
    and elevations in Mare2D format.

    Parameters
    ----------
    o2d_data : mtpy.modeling.occam2d.Data
        Occam2D Data object that has been initialised with EDI data.
    surface_file : str or bytes
        Full path to ASCII grid file containing topography data.
    elevation_sample_n : int, optional
        Number of points to sample from the Occam2D profile for
        for loading into Mare2D. According to original script this
        should be changed to something sensible based on length of the
        profile line.
    flip_elevation : bool, optional
        If True, elevation is multiplied by -1 after interpolation.

    Returns
    -------
    tuple of float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray
        A tuple containing (mare origin x, mare origin y, utm_zone,
        site locations, site elevations, Mare2DEM profile,
        profile elevation). Mare origin X and Y are the middle of the
        profile line in UTM coordinates. utm_zone is the UTM string of
        the Mare2D profile origin. Site locations are the locations of
        the site in Mare2D coordinates. Site elevations are elevation
        of sites in metres (interpolated from surface_file). Mare2DEM
        profile is the full profile line in Mare2D coordinates.
        Elevation is the elevation of the Mare2D profile line.
    """
    # Get site location eastings and northings from Occam2D profile
    site_easts = []
    site_norths = []
    for i in range(len(o2d_data.station_list)):
        site_easts.append(o2d_data.edi_list[i].projected_east)
        site_norths.append(o2d_data.edi_list[i].projected_north)
    site_easts = np.array(site_easts)
    site_norths = np.array(site_norths)

    # Extract information from Occam2D profile
    m, c1 = o2d_data.profile_line
    # Find start of profile
    # Assume that the start of the profile always has minimum eastings
    x0 = site_easts.min()
    y0 = m * x0 + c1
    # Find end of the profile
    # Assume y1 as the maximum northing
    x1 = site_easts.max()
    y1 = m * x1 + c1
    # Mare2D origin (middle of profile line)
    mare_origin_x = (site_easts.min() + site_easts.max()) / 2
    mare_origin_y = (site_norths.min() + site_norths.max()) / 2
    profile_length = line_length(x0, y0, x1, y1)
    # UTM zone of Mare2D profile
    epsg = o2d_data.model_epsg
    projected_origin = gis_tools.project_point_utm2ll(mare_origin_x, mare_origin_y,
                                                      None, epsg=epsg)
    utm_zone = gis_tools.get_utm_zone(projected_origin[0], projected_origin[1])[2]

    # Select samples from Occam2D profile for loading into mare2dem
    # For whatever reason, original script ignores first and last coordinates (x0, y0), (x1, y1)
    o2d_easts = np.delete(np.linspace(x0, x1, elevation_sample_n, endpoint=False), 0)
    o2d_norths = np.delete(np.linspace(y0, y1, elevation_sample_n, endpoint=False), 0)
    # Add exact site locations
    o2d_easts = np.concatenate((o2d_easts, site_easts))
    o2d_norths = np.concatenate((o2d_norths, site_norths))
    # Make sure station indices align between east and north arrays
    sort_inds = np.argsort(o2d_easts)
    o2d_easts = o2d_easts[sort_inds]
    o2d_norths = o2d_norths[sort_inds]

    # Interpolate elevation across profile points as a grid
    elevation = mesh_tools.interpolate_elevation_to_grid(
        o2d_easts, o2d_norths, epsg=epsg, surfacefile=surface_file,
        method='cubic')
    elevation = elevation * -1 if flip_elevation else elevation

    # Get profile elevation (which is a straight line through the elevation grid)
    profile_elevation = []
    for i in range(len(o2d_easts)):
        profile_elevation.append(elevation[i, i])
    profile_elevation = np.array(profile_elevation)

    # Convert profile to Mare2D system
    m2d_profile = points_o2d_to_m2d(o2d_easts, o2d_norths, profile_length)

    # Get elevation at sites
    site_elevations = []
    for i, site_east in enumerate(site_easts):
        j = np.where(o2d_easts == site_east)
        if o2d_norths[j] == site_norths[i]:
            site_elevations.append(elevation[j, j])
    site_elevations = np.squeeze(np.array(site_elevations, dtype=np.float64))

    # Convert site locations to Mare2D system (0 is profile middle)
    site_locations = points_o2d_to_m2d(site_easts, site_norths, profile_length)

    return (mare_origin_x, mare_origin_y, utm_zone, site_locations, site_elevations,
            m2d_profile, profile_elevation)


def write_elevation_file(m2d_profile, profile_elevation, savepath=None):
    """
    Write an elevation file that can be imported into Mamba2D.
    Currently broken due to shape mismatch.

    m2d_profile : np.ndarray
        1D array of the m2d profile in m2d coordinates.
    profile_elevation : np.ndarray
        Z values for each point in the profile.
    savepath : str or bytes, optional
        Full path including file of where to save elevation file.
    """
    elevation_model = np.stack((m2d_profile, profile_elevation), axis=1).astype(np.float64)
    if savepath is None:
        savepath = os.path.join(os.getcwd(), 'elevation.txt')
    np.savetxt(savepath, elevation_model)


def write_mare2dem_data(o2d_filepath, site_locations, site_elevations,
                        mare_origin, gstrike, solve_statics=False, savepath=None):
    """
    Uses an Occam2D data file and site locations + elevations to
    generate a MARE2DEM data file.

    Parameters
    ----------
    o2d_filepath : bytes or str
        Full path to the Occam2D data file created from EDI data.
    site_locations : np.ndarray
        Array of shape (nstations). According to the original comments
        in the `EDI2Mare2DEM_withOccam2D_new` script, this is the
        site elevation but in Mare2D coordinates, however it gets set
        as the 'y' component of receiver locations in the Mare2DEM
        data file.
    site_elevations : np.ndarray
        Array of shape (nstations). I think this is the site elevation
        in UTM coordinates. Gets set as the 'z' component of receiver
        locations in the Mare2DEM data file.
    mare_origin : tuple
        Tuple of float (x, y, utm_zone). The Mare2D origin in UTM c
        oordinates. Note according to the original script, Mare2D origin
        is the middle of the profile line. utm_zone is the UTM string,
        e.g. '54S'.
    gstrike : int
        The 2D strike, same as
        `geoelectric_strike` in Occam2D model. This information is used
        for the UTM origin line in the data file.
    solve_statics : bool or list, optional
        If boolean, sets whether to solve statics for all stations.
        A list of site names can be passed. Sites in the list will
        have solve_statics set to True, any sites not included
        are set to False. This writes a '1' or a '0' in
        the Receiver section of the Mare2D data file in the SovleStatics
        column.
    savepath : bytes or str, optional
        Full path of where to save the Mare2D data file. If not
        provided, will be saved as 'Mare2D_data.txt' in working
        directory.
    """
    # Prepare O2D data for data block
    o2d_sites = []
    o2d_freqs = []
    o2d_types = []
    o2d_datums = []
    o2d_errors = []
    with open(o2d_filepath, 'r') as f:
        read_data = f.readlines()
        reading_data = False
        for line in read_data:
            if line.startswith('SITE '):
                reading_data = True
                continue
            elif reading_data:
                parts = line.split()
                o2d_sites.append(parts[0])
                o2d_freqs.append(parts[1])
                o2d_types.append(parts[2])
                o2d_datums.append(parts[3])
                o2d_errors.append(parts[4])

    sites = np.array(o2d_sites, dtype=np.int8)
    freqs = np.array(o2d_freqs, dtype=np.int8)
    types = np.array(o2d_types, dtype=np.int8)
    datums = np.array(o2d_datums, dtype=np.float64)
    errors = np.array(o2d_errors, dtype=np.float64)
    # Convert occam2d types to mare2d types
    # The below is: for each element in types array, return corresponding element in conversion
    # dict, if not found in dict return original element
    type_conversion = {1: 123, 2: 104, 3: 133, 4: 134, 5: 125, 6: 106, 9: 103, 10: 105}
    types = np.vectorize(lambda x: type_conversion.get(x, x))(types)
    # Put into dataframe for easier stringifying
    # Note: TX# == RX# == site ID for MT stations
    data_df = pd.DataFrame((types, freqs, sites, sites, datums, errors), dtype=np.object).T
    # Bit of a hack: add the '!' to the data frame header because the 'type' integer is small
    # enough that the 'Type' header will have no left whitespace padding, so we can't prepend
    # it with '!' without throwing off the alignment.
    data_df.columns = ['! Type', 'Freq #', 'Tx #', 'Rx #', 'Data', 'StdErr']
    data_str = data_df.to_string(index=False, float_format=lambda x: '%.4f' % x)

    # Prepare data for the Reciever block
    # Zeros of shape (n_sites) for X (as float), Theta, Alpha, Beta and Length (ints) columns
    x_col = np.zeros(site_locations.shape, dtype=np.float64)
    zero_ints = np.zeros(site_locations.shape, dtype=np.int8)
    t_col, a_col, b_col, l_col = zero_ints, zero_ints, zero_ints, zero_ints
    # add 0.1 m (shift the sites 10 cm beneath subsurface as recommended)
    site_elevations += 0.1
    # According to original script, need to reread the Occam2D file to get stations in the correct
    # order
    reread_o2d = o2d.Data()
    reread_o2d.read_data_file(o2d_filepath)
    site_names = np.array([sn.split('_')[0] for sn in reread_o2d.station_list])
    if isinstance(solve_statics, bool):
        statics = np.ones(site_locations.shape, dtype=np.int8) if solve_statics else zero_ints
    else:
        statics = np.zeros(site_locations.shape)
        for sn in solve_statics:
            statics[np.where(site_names == sn)] = 1
    # Put into dataframe for easier stringifying
    recv_df = pd.DataFrame((x_col, site_locations, site_elevations, t_col, a_col, b_col, l_col,
                            statics, site_names)).T
    recv_df.columns = ['X', 'Y', 'Z', 'Theta', 'Alpha', 'Beta', 'Length', 'SolveStatic', 'Name']
    recv_str = list(recv_df.to_string(index=False, float_format=lambda x: '%.6f' % x))
    # Replace the first char of header with Mare2DEM comment symbol '!'
    # This way the header is correct but Pandas handles the alignment and spacing
    recv_str[0] = '!'
    recv_str = "".join(recv_str)

    if savepath is None:
        savepath = os.path.join(os.getcwd(), 'Mare2D_data.txt')
    with open(savepath, 'w') as output:
        # 1. header
        fstring = 'Format:  EMData_2.2\n'
        fstring += 'UTM of x,y origin (UTM zone, N, E, 2D strike):'
        # TODO: fix hardocded UTM zone
        fstring += ' {:s}{:>13.1f}{:>13.1f}\t{:d}\n'.format(
            mare_origin[2], mare_origin[0], mare_origin[1], gstrike)

        # 2. frequencies
        fstring += '# MT Frequencies:    {}\n'.format(len(reread_o2d.freq))
        fstring += '\n'.join([str(round(f, 8)) for f in reread_o2d.freq])

        # 3. receiver info
        fstring += '\n# MT Receivers:      {}\n'.format(len(site_names))
        fstring += recv_str
        fstring += '\n'

        # 4. data
        fstring += '# Data:       {}\n'.format(len(datums))
        fstring += data_str

        output.write(fstring)
