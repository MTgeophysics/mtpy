##==============================================================================
## Add in elevation to the model
##==============================================================================
#
##--> read in ascii dem file
# def read_dem_ascii(ascii_fn, cell_size=500, model_center=(0, 0), rot_90=0):
#    """
#    read in dem which is ascii format
#
#    The ascii format is assumed to be:
#    ncols         3601
#    nrows         3601
#    xllcorner     -119.00013888889
#    yllcorner     36.999861111111
#    cellsize      0.00027777777777778
#    NODATA_value  -9999
#    elevation data W --> E
#    N
#    |
#    V
#    S
#    """
#    dfid = file(ascii_fn, 'r')
#    d_dict = {}
#    for ii in range(6):
#        dline = dfid.readline()
#        dline = dline.strip().split()
#        key = dline[0].strip().lower()
#        value = float(dline[1].strip())
#        d_dict[key] = value
#
#    x0 = d_dict['xllcorner']
#    y0 = d_dict['yllcorner']
#    nx = int(d_dict['ncols'])
#    ny = int(d_dict['nrows'])
#    cs = d_dict['cellsize']
#
#    # read in the elevation data
#    elevation = np.zeros((nx, ny))
#
#    for ii in range(1, int(ny)+2):
#        dline = dfid.readline()
#        if len(str(dline)) > 1:
#            #needs to be backwards because first line is the furthest north row.
#            elevation[:, -ii] = np.array(dline.strip().split(' '), dtype='float')
#        else:
#            break
#
#    dfid.close()
#
#    # create lat and lon arrays from the dem fle
#    lon = np.arange(x0, x0+cs*(nx), cs)
#    lat = np.arange(y0, y0+cs*(ny), cs)
#
#    # calculate the lower left and uper right corners of the grid in meters
#    ll_en = gis_tools.project_point_ll2utm(lat[0], lon[0])
#    ur_en = gis_tools.project_point_ll2utm(lat[-1], lon[-1])
#
#    # estimate cell sizes for each dem measurement
#    d_east = abs(ll_en[0]-ur_en[0])/nx
#    d_north = abs(ll_en[1]-ur_en[1])/ny
#
#    # calculate the number of new cells according to the given cell size
#    # if the given cell size and cs are similar int could make the value 0,
#    # hence the need to make it one if it is 0.
#    num_cells = max([1, int(cell_size/np.mean([d_east, d_north]))])
#
#    # make easting and northing arrays in meters corresponding to lat and lon
#    east = np.arange(ll_en[0], ur_en[0], d_east)
#    north = np.arange(ll_en[1], ur_en[1], d_north)
#
#    #resample the data accordingly
#    new_east = east[np.arange(0, east.size, num_cells)]
#    new_north = north[np.arange(0, north.size, num_cells)]
#    new_x, new_y = np.meshgrid(np.arange(0, east.size, num_cells),
#                               np.arange(0, north.size, num_cells),
#                               indexing='ij')
#    elevation = elevation[new_x, new_y]
#    # make any null values set to minimum elevation, could be dangerous
#    elevation[np.where(elevation == -9999.0)] = elevation[np.where(elevation != -9999.0)].min()
#
#    # estimate the shift of the DEM to relative model coordinates
#    mid_east = np.where(new_east >= model_center[0])[0][0]
#    mid_north = np.where(new_north >= model_center[1])[0][0]
#
#    new_east -= new_east[mid_east]
#    new_north -= new_north[mid_north]
#
#    # need to rotate cause I think I wrote the dem backwards
#    if rot_90 == 1 or rot_90 == 3:
#        elevation = np.rot90(elevation, rot_90)
#        return new_north, new_east, elevation
#    else:
#        elevation = np.rot90(elevation, rot_90)
#
#        return new_east, new_north, elevation
#
# def interpolate_elevation(elev_east, elev_north, elevation, model_east,
#                          model_north, pad=3):
#    """
#    interpolate the elevation onto the model grid.
#
#    Arguments:
#    ---------------
#
#        *elev_east* : np.ndarray(num_east_nodes)
#                      easting grid for elevation model
#
#        *elev_north* : np.ndarray(num_north_nodes)
#                      northing grid for elevation model
#
#        *elevation* : np.ndarray(num_east_nodes, num_north_nodes)
#                     elevation model assumes x is east, y is north
#                     Units are meters
#
#        *model_east* : np.ndarray(num_east_nodes_model)
#                     relative easting grid of resistivity model
#
#        *model_north* : np.ndarray(num_north_nodes_model)
#                     relative northin grid of resistivity model
#
#        *pad* : int
#                number of cells to repeat elevation model by.  So for pad=3,
#                then the interpolated elevation model onto the resistivity
#                model grid will have the outer 3 cells will be repeats of
#                the adjacent cell.  This is to extend the elevation model
#                to the resistivity model cause most elevation models will
#                not cover the entire area.
#
#    Returns:
#    --------------
#
#        *interp_elev* : np.ndarray(num_north_nodes_model, num_east_nodes_model)
#                        the elevation model interpolated onto the resistivity
#                        model grid.
#
#    """
#    # need to line up the elevation with the model
#    grid_east, grid_north = np.broadcast_arrays(elev_east[:, None],
#                                                elev_north[None, :])
#    # interpolate onto the model grid
#    interp_elev = spi.griddata((grid_east.ravel(), grid_north.ravel()),
#                               elevation.ravel(),
#                               (model_east[:, None],
#                                model_north[None, :]),
#                                method='linear',
#                                fill_value=elevation.mean())
#
#    interp_elev[0:pad, pad:-pad] = interp_elev[pad, pad:-pad]
#    interp_elev[-pad:, pad:-pad] = interp_elev[-pad-1, pad:-pad]
#    interp_elev[:, 0:pad] = interp_elev[:, pad].repeat(pad).reshape(
#                                                interp_elev[:, 0:pad].shape)
#    interp_elev[:, -pad:] = interp_elev[:, -pad-1].repeat(pad).reshape(
#                                                interp_elev[:, -pad:].shape)
#
#    # transpose the modeled elevation to align with x=N, y=E
#    interp_elev = interp_elev.T
#
#    return interp_elev
#
# def make_elevation_model(interp_elev, model_nodes_z, elevation_cell=30,
#                         pad=3, res_air=1e12, fill_res=100, res_sea=0.3):
#    """
#    Take the elevation data of the interpolated elevation model and map that
#    onto the resistivity model by adding elevation cells to the existing model.
#
#    ..Note: that if there are large elevation gains, the elevation cell size
#            might need to be increased.
#
#    Arguments:
#    -------------
#        *interp_elev* : np.ndarray(num_nodes_north, num_nodes_east)
#                        elevation model that has been interpolated onto the
#                        resistivity model grid. Units are in meters.
#
#        *model_nodes_z* : np.ndarray(num_z_nodes_of_model)
#                          vertical nodes of the resistivity model without
#                          topography.  Note these are the nodes given in
#                          relative thickness, not the grid, which is total
#                          depth.  Units are meters.
#
#        *elevation_cell* : float
#                           height of elevation cells to be added on.  These
#                           are assumed to be the same at all elevations.
#                           Units are in meters
#
#        *pad* : int
#                number of cells to look for maximum and minimum elevation.
#                So if you only want elevations within the survey area,
#                set pad equal to the number of padding cells of the
#                resistivity model grid.
#
#        *res_air* : float
#                    resistivity of air.  Default is 1E12 Ohm-m
#
#        *fill_res* : float
#                     resistivity value of subsurface in Ohm-m.
#
#    Returns:
#    -------------
#        *elevation_model* : np.ndarray(num_north_nodes, num_east_nodes,
#                                       num_elev_nodes+num_z_nodes)
#                         Model grid with elevation mapped onto it.
#                         Where anything above the surface will be given the
#                         value of res_air, everything else will be fill_res
#
#        *new_nodes_z* : np.ndarray(num_z_nodes+num_elev_nodes)
#                        a new array of vertical nodes, where any nodes smaller
#                        than elevation_cell will be set to elevation_cell.
#                        This can be input into a modem.Model object to
#                        rewrite the model file.
#
#    """
#
#    # calculate the max elevation within survey area
#    elev_max = interp_elev[pad:-pad, pad:-pad].max()
#
#    # need to set sea level to 0 elevation
#    elev_min = max([0, interp_elev[pad:-pad, pad:-pad].min()])
#
#    # scale the interpolated elevations to fit within elev_max, elev_min
#    interp_elev[np.where(interp_elev > elev_max)] = elev_max
#    #interp_elev[np.where(interp_elev < elev_min)] = elev_min
#
#    # calculate the number of elevation cells needed
#    num_elev_cells = int((elev_max-elev_min)/elevation_cell)
#    print 'Number of elevation cells: {0}'.format(num_elev_cells)
#
#    # find sea level if it is there
#    if elev_min < 0:
#        sea_level_index = num_elev_cells-abs(int((elev_min)/elevation_cell))-1
#    else:
#        sea_level_index = num_elev_cells-1
#
#    print 'Sea level index is {0}'.format(sea_level_index)
#
#
#    # make an array of just the elevation for the model
#    # north is first index, east is second, vertical is third
#    elevation_model = np.ones((interp_elev.shape[0],
#                               interp_elev.shape[1],
#                               num_elev_cells+model_nodes_z.shape[0]))
#
#    elevation_model[:, :, :] = fill_res
#
#
#
#    # fill in elevation model with air values.  Remeber Z is positive down, so
#    # the top of the model is the highest point and index 0 is highest
#    # elevation
#    for nn in range(interp_elev.shape[0]):
#        for ee in range(interp_elev.shape[1]):
#            # need to test for ocean
#            if interp_elev[nn, ee] < 0:
#                # fill in from bottom to sea level, then rest with air
#                elevation_model[nn, ee, 0:sea_level_index] = res_air
#                dz = sea_level_index+abs(int((interp_elev[nn, ee])/elevation_cell))+1
#                elevation_model[nn, ee, sea_level_index:dz] = res_sea
#            else:
#                dz = int((elev_max-interp_elev[nn, ee])/elevation_cell)
#                elevation_model[nn, ee, 0:dz] = res_air
#
#    # make new z nodes array
#    new_nodes_z = np.append(np.repeat(elevation_cell, num_elev_cells),
#                            model_nodes_z)
#
#    new_nodes_z[np.where(new_nodes_z < elevation_cell)] = elevation_cell
#
#    return elevation_model, new_nodes_z
#
# def add_topography_to_model(dem_ascii_fn, model_fn, model_center=(0,0),
#                            rot_90=0, cell_size=500, elev_cell=30, pad=1):
#    """
#    Add topography to an existing model from a dem in ascii format.
#
#    The ascii format is assumed to be:
#    ncols         3601
#    nrows         3601
#    xllcorner     -119.00013888889
#    yllcorner     36.999861111111
#    cellsize      0.00027777777777778
#    NODATA_value  -9999
#    elevation data W --> E
#    N
#    |
#    V
#    S
#
#    Arguments:
#    -------------
#        *dem_ascii_fn* : string
#                         full path to ascii dem file
#
#        *model_fn* : string
#                     full path to existing ModEM model file
#
#        *model_center* : (east, north) in meters
#                         Sometimes the center of the DEM and the center of the
#                         model don't line up.  Use this parameter to line
#                         everything up properly.
#
#        *rot_90* : [ 0 | 1 | 2 | 3 ]
#                   rotate the elevation model by rot_90*90 degrees.  Sometimes
#                   the elevation model is flipped depending on your coordinate
#                   system.
#
#        *cell_size* : float (meters)
#                      horizontal cell size of grid to interpolate elevation
#                      onto.  This should be smaller or equal to the input
#                      model cell size to be sure there is not spatial aliasing
#
#        *elev_cell* : float (meters)
#                      vertical size of each elevation cell.  This value should
#                      be about 1/10th the smalles skin depth.
#
#    Returns:
#    ---------------
#        *new_model_fn* : string
#                         full path to model file that contains topography
#
#    """
#     ### 1.) read in the dem and center it onto the resistivity model
#    e_east, e_north, elevation = read_dem_ascii(dem_ascii_fn,
#                                                cell_size=cell_size,
#                                                model_center=model_center,
#                                                rot_90=rot_90)
#    m_obj = Model()
#    m_obj.read_model_file(model_fn)
#    ### 2.) interpolate the elevation model onto the model grid
#    m_elev = interpolate_elevation(e_east, e_north, elevation,
#                                   m_obj.grid_east, m_obj.grid_north, pad=pad)
#
#    m_elev[np.where(m_elev == -9999.0)] = m_elev[np.where(m_elev != -9999.0)].min()
#    ### 3.) make a resistivity model that incoorporates topography
#    mod_elev, elev_nodes_z = make_elevation_model(m_elev, m_obj.nodes_z,
#                                                  elevation_cell=elev_cell)
#
#    ### 4.) write new model file
#    m_obj.nodes_z = elev_nodes_z
#    m_obj.res_model = mod_elev
#    m_obj.model_fn = None
#    m_obj.save_path = os.path.dirname(model_fn)
#    m_obj.write_model_file(model_fn_basename='{0}_topo.rho'.format(
#                           os.path.basename(model_fn)[0:-4]))
#
#    return m_obj.model_fn
#
# def change_data_elevation(data_fn, model_fn, new_data_fn=None, res_air=1e12):
#    """
#    At each station in the data file rewrite the elevation, so the station is
#    on the surface, not floating in air.
#
#    Arguments:
#    ------------------
#        *data_fn* : string
#                    full path to a ModEM data file
#
#        *model_fn* : string
#                    full path to ModEM model file that has elevation
#                    incoorporated.
#
#        *new_data_fn* : string
#                        full path to new data file name.  If None, then
#                        new file name will add _elev.dat to input filename
#
#        *res_air* : float
#                    resistivity of air.  Default is 1E12 Ohm-m
#    Returns:
#    -------------
#        *new_data_fn* : string
#                        full path to new data file.
#    """
#
#    d_obj = Data()
#    d_obj.read_data_file(data_fn)
#
#    m_obj = Model()
#    m_obj.read_model_file(model_fn)
#
#    s_locations = d_obj.station_locations.station_locations.copy()
#
#    # need to subtract one because we are finding the cell next to it
#    for s_arr in s_locations:
#        e_index = np.where(m_obj.grid_east >= s_arr['rel_east'])[0][0]-1
#        n_index = np.where(m_obj.grid_north >= s_arr['rel_north'])[0][0]-1
#        z_index = np.where(m_obj.res_model[n_index, e_index, :] < res_air*.9)[0][0]
#        s_index = np.where(d_obj.data_array['station']==s_arr['station'])[0][0]
#        d_obj.data_array[s_index]['elev'] = m_obj.grid_z[z_index]
#
#        print s_arr['station'], s_arr['elev'], n_index, e_index, z_index
#        print s_arr['rel_north'], s_arr['rel_east']
#        print m_obj.grid_north[n_index], m_obj.grid_east[e_index]
#        print '-'*20
##
##    for key in d_obj.mt_dict.keys():
##        mt_obj = d_obj.mt_dict[key]
##        e_index = np.where(m_obj.grid_east > mt_obj.grid_east)[0][0]
##        n_index = np.where(m_obj.grid_north > mt_obj.grid_north)[0][0]
##        z_index = np.where(m_obj.res_model[n_index, e_index, :] < res_air*.9)[0][0]
##        s_index = np.where(d_obj.data_array['station']==key)[0][0]
##        d_obj.data_array[s_index]['elev'] = m_obj.grid_z[z_index]
##
##        mt_obj.grid_elev = m_obj.grid_z[z_index]
##
#    if new_data_fn is None:
#        new_dfn = '{0}{1}'.format(data_fn[:-4], '_elev.dat')
#    else:
#        new_dfn=new_data_fn
#
#    d_obj.write_data_file(save_path=os.path.dirname(new_dfn),
#                          fn_basename=os.path.basename(new_dfn),
#                          compute_error=False,
#                          fill=False,
#                          elevation=True)
#
#    return new_dfn
#
# def center_stations(data_fn, model_fn, new_data_fn=None):
#    """
#    center station locations to the middle of cells, might be useful for
#    topography.
#    """
#
#    d_obj = Data()
#    d_obj.read_data_file(data_fn)
#
#    m_obj = Model()
#    m_obj.read_model_file(model_fn)
#
#    for s_arr in d_obj.station_locations.station_locations:
#        e_index = np.where(m_obj.grid_east >= s_arr['rel_east'])[0][0]-1
#        n_index = np.where(m_obj.grid_north >= s_arr['rel_north'])[0][0]-1
#
#        mid_east = m_obj.grid_east[e_index:e_index+2].mean()
#        mid_north = m_obj.grid_north[n_index:n_index+2].mean()
#
#        s_index = np.where(d_obj.data_array['station']==s_arr['station'])[0][0]
#
#        d_obj.data_array[s_index]['rel_east'] = mid_east
#        d_obj.data_array[s_index]['rel_north'] = mid_north
#
#        print s_arr['rel_east'], s_arr['rel_north']
#        print mid_east, mid_north
#        print '-'*30
#
#    if new_data_fn is None:
#        new_dfn = '{0}{1}'.format(data_fn[:-4], '_center.dat')
#    else:
#        new_dfn=new_data_fn
#
#    d_obj.write_data_file(save_path=os.path.dirname(new_dfn),
#                          fn_basename=os.path.basename(new_dfn),
#                          compute_error=False,
#                          fill=False,
#                          elevation=True)
#
#    return new_dfn
#
