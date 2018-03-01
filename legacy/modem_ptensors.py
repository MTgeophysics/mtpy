#!/bin/env python
"""
Description:
    Class for loading MODEM data and plotting phase tensors.

References:
 
CreationDate:   11/24/17
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     11/24/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""
import os
import mtpy.modeling.modem as modem
import numpy as np
import logging, traceback

from mtpy.analysis.pt import ResidualPhaseTensor

from matplotlib.patches import Ellipse

class ModEM_ptensors:
    def __init__(self, data_fn, resp_fn=None, read_mode=1):
        '''
        Class for loading MODEM data and plotting phase tensors.

        :param data_fn: MODEM data file name
        '''

        self._data_fn = data_fn
        self._resp_fn = resp_fn
        
        try:
            self._data_obj = modem.Data()
            self._data_obj.read_data_file(self._data_fn)
            self._modem_obj = modem.Data()
            self._modem_obj.read_data_file(self._data_fn)
        except Exception as err:
            print 'Failed to read %s' % (self._data_fn)
            logging.error(traceback.format_exc())
            exit(-1)

        if self._resp_fn is not None:
            try:
                self._resp_obj = modem.Data()
                self._resp_obj.read_data_file(self._resp_fn)
            except Exception as err:
                print 'Failed to read %s' % (self._resp_fn)
                logging.error(traceback.format_exc())
                exit(-1)
        else:
            self._resp_obj = None


        self._plot_period = self._data_obj.period_list.copy()
        self._mt_obj_list = [self._data_obj.mt_dict[key]
                             for key in self._data_obj.mt_dict.keys()]

        # Read data
        self._pt_dict = {}
        self._ptol = 0.05
        
        self.read_mode = read_mode
        if read_mode==1:
            self.get_pt()
        else:
            self.get_pt2()
        
    
        
    def get_pt(self):
        
        for p_index,plot_per in enumerate(self._plot_period):
            self._pt_dict[plot_per] = []
            for mt_obj in self._mt_obj_list:
                p_index = [ff for ff, f2 in enumerate(1. / mt_obj.Z.freq)
                           if (f2 > plot_per * (1 - self._ptol)) and
                           (f2 < plot_per * (1 + self._ptol))][0]

                s_ind = np.where(self._modem_obj.station_locations.station==mt_obj.station)
                rel_east = self._modem_obj.station_locations.rel_east[s_ind]
                rel_north = self._modem_obj.station_locations.rel_north[s_ind]
                
                pt_tuple = (mt_obj.station, mt_obj.lon, mt_obj.lat,
                            rel_east,rel_north,
                            mt_obj.pt.phimin[p_index],
                            mt_obj.pt.phimax[p_index],
                            mt_obj.pt.azimuth[p_index],
                            mt_obj.pt.beta[p_index],
                            2 * mt_obj.pt.beta[p_index],
                            mt_obj.pt.ellipticity[p_index])

                self._pt_dict[plot_per].append(pt_tuple)
            # end for
            self._pt_dict[plot_per] = np.array(self._pt_dict[plot_per],
                                               dtype=[('station', '|S15'),
                                                      ('lon', np.float),
                                                      ('lat', np.float),
                                                      ('rel_east', np.float),
                                                      ('rel_north', np.float),
                                                      ('phimin', np.float),
                                                      ('phimax', np.float),
                                                      ('azimuth', np.float),
                                                      ('skew', np.float),
                                                      ('n_skew', np.float),
                                                      ('ellipticity', np.float)])


    def get_pt2(self):
        """
        put pt parameters into something useful for plotting
        """

        ns = len(self._data_obj.mt_dict.keys())
        nf = len(self._data_obj.period_list)

        data_pt_arr = np.zeros((nf, ns), dtype=[('lon', np.float),
                                                ('lat', np.float),
                                                ('phimin', np.float),
                                                ('phimax', np.float),
                                                ('skew', np.float),
                                                ('azimuth', np.float),
                                                ('rel_east', np.float),
                                                ('rel_north', np.float)])
        if self._resp_obj is not None:
            model_pt_arr = np.zeros((nf, ns), dtype=[('lon', np.float),
                                                ('lat', np.float),
                                                ('phimin', np.float),
                                                ('phimax', np.float),
                                                ('skew', np.float),
                                                ('azimuth', np.float),
                                                ('rel_east', np.float),
                                                ('rel_north', np.float)])

            res_pt_arr = np.zeros((nf, ns), dtype=[('lon', np.float),
                                                   ('lat', np.float),
                                                   ('phimin', np.float),
                                                   ('phimax', np.float),
                                                   ('skew', np.float),
                                                   ('azimuth', np.float),
                                                   ('rel_east', np.float),
                                                   ('rel_north', np.float),
                                                   ('geometric_mean', np.float)])

        for key in self._data_obj.mt_dict.keys():
            mt_obj = self._data_obj.mt_dict[key]
            ii = np.where(self._data_obj.station_locations.station==mt_obj.station)[0][0]
            
            east = self._data_obj.mt_dict[key].grid_east
            north = self._data_obj.mt_dict[key].grid_north
            lon = self._data_obj.mt_dict[key].lon
            lat = self._data_obj.mt_dict[key].lat

            dpt = self._data_obj.mt_dict[key].pt
            data_pt_arr[:, ii]['lon'] = lon
            data_pt_arr[:, ii]['lat'] = lat
            data_pt_arr[:, ii]['rel_east'] = east
            data_pt_arr[:, ii]['rel_north'] = north
            data_pt_arr[:, ii]['phimin'] = dpt.phimin
            data_pt_arr[:, ii]['phimax'] = dpt.phimax
            data_pt_arr[:, ii]['azimuth'] = dpt.azimuth
            data_pt_arr[:, ii]['skew'] = dpt.beta
            if self._resp_obj is not None:
                mpt = self._resp_obj.mt_dict[key].pt
                try:
                    rpt = ResidualPhaseTensor(pt_object1=dpt,
                                                   pt_object2=mpt)
                    rpt = rpt.residual_pt
                    res_pt_arr[:, ii]['lon'] = lon
                    res_pt_arr[:, ii]['lat'] = lat
                    res_pt_arr[:, ii]['rel_east'] = east
                    res_pt_arr[:, ii]['rel_north'] = north
                    res_pt_arr[:, ii]['phimin'] = rpt.phimin
                    res_pt_arr[:, ii]['phimax'] = rpt.phimax
                    res_pt_arr[:, ii]['azimuth'] = rpt.azimuth
                    res_pt_arr[:, ii]['skew'] = rpt.beta
                    res_pt_arr[:, ii]['geometric_mean'] = np.sqrt(abs(rpt.phimin[0] * \
                                                                      rpt.phimax[0]))
                except:
                    pass
                
                model_pt_arr[:, ii]['lon'] = lon
                model_pt_arr[:, ii]['lat'] = lat
                model_pt_arr[:, ii]['rel_east'] = east
                model_pt_arr[:, ii]['rel_north'] = north
                model_pt_arr[:, ii]['phimin'] = mpt.phimin
                model_pt_arr[:, ii]['phimax'] = mpt.phimax
                model_pt_arr[:, ii]['azimuth'] = mpt.azimuth
                model_pt_arr[:, ii]['skew'] = mpt.beta

        # make these attributes
        self.pt_data_arr = data_pt_arr
        if self._resp_obj is not None:
            self.pt_resp_arr = model_pt_arr
            self.pt_resid_arr = res_pt_arr
                                                      
                                                      
        # end for
    # end func

    def get_period_attributes(self, periodIdx, key, ptarray='data'):
        '''
        Returns, for a given period, a list of attribute values for key
        (e.g. skew, phimax, etc.).

        :param periodIdx: index of period; print out _plot_period for periods available
        :param key: attribute key
        :return: numpy array of attribute values
        '''
        assert (periodIdx >= 0 and periodIdx < len(self._plot_period)), \
            'Error: Index for plot-period out of bounds.'

        if self.read_mode == 1:
            pk = self._pt_dict.keys()[periodIdx]
            try:
                vals = self._pt_dict[pk][key]
                return vals
            except Exception as err:
                logging.error(traceback.format_exc())
        else:
            pk = periodIdx
            vals = getattr(self,'pt_'+ptarray+'_arr')[pk][key]
            
            return vals


    # end func

    def plot(self, ax, m, periodIdx, ellipse_size_factor=10000,
             cvals=None, **kwargs):

        '''
        Plots phase tensors for a given period index.

        :param ax: plot axis
        :param m: basemap instance
        :param periodIdx: period index
        :param ellipse_size_factor: factor to control ellipse size
        :param cvals: list of colour values for colouring each ellipse; must be of
                      the same length as the number of tuples for each period
        :param kwargs: list of relevant matplotlib arguments (e.g. zorder, alpha, etc.)
        '''

        assert (periodIdx >= 0 and periodIdx < len(self._plot_period)), \
            'Error: Index for plot-period out of bounds.'

        k = self._pt_dict.keys()[periodIdx]
        for i in range(len(self._pt_dict[k])):
            lon = self._pt_dict[k]['lon'][i]
            lat = self._pt_dict[k]['lat'][i]
            phimax = self._pt_dict[k]['phimax'][i] / self._pt_dict[k]['phimax'].max()
            phimin = self._pt_dict[k]['phimin'][i] / self._pt_dict[k]['phimax'].max()
            az = self._pt_dict[k]['azimuth'][i]
            nskew = self._pt_dict[k]['n_skew'][i]

            # print az
            if (phimax > 0 and phimin > 0):
                c = None
                if (cvals is not None): c = cvals[i]
                if (c is not None): kwargs['facecolor'] = c

                x, y = m(lon, lat)
                
                e = Ellipse([x, y],
                            phimax * ellipse_size_factor,
                            phimin * ellipse_size_factor,
                            az, **kwargs)
                ax.add_artist(e)
            # end if
        # end for
    # end func

    def plot2(self, ax, m, periodIdx, param='data' ,ellipse_size_factor=10000,
              cvals=None, map_scale='m', centre_shift=[0,0], **kwargs):
    
        '''
        Plots phase tensors for a given period index.

        :param ax: plot axis
        :param m: basemap instance
        :param periodIdx: period index
        :param ellipse_size_factor: factor to control ellipse size
        :param cvals: list of colour values for colouring each ellipse; must be of
                      the same length as the number of tuples for each period
        :param map_scale: map length scale
        :param kwargs: list of relevant matplotlib arguments (e.g. zorder, alpha, etc.)
        '''

        assert (periodIdx >= 0 and periodIdx < len(self._plot_period)), \
            'Error: Index for plot-period out of bounds.'

        if self.read_mode == 2:
            k = periodIdx
            pt_array = getattr(self,'pt_'+param+'_arr')
        else:
            k = self._pt_dict.keys()[periodIdx]
            pt_array = self._pt_dict

        for i in range(len(pt_array[k])):
            lon = pt_array[k]['lon'][i]
            lat = pt_array[k]['lat'][i]
            phimax = pt_array[k]['phimax'][i] #/ pt_array[k]['phimax'].max()
            phimin = pt_array[k]['phimin'][i] #/ pt_array[k]['phimax'].max()
            az = pt_array[k]['azimuth'][i]
            if param == 'resid':
                phimin = np.abs(phimin)
            nskew = pt_array[k]['skew'][i]

            # print az
            if (phimax > 0 and phimin > 0):
                c = None
                if (cvals is not None): c = cvals[i]
                if (c is not None): kwargs['facecolor'] = c

                if m is None:
                    x = pt_array[k]['rel_east'][i]
                    y = pt_array[k]['rel_north'][i]
                    if map_scale == 'km':
                        x /= 1e3
                        y /= 1e3
                else:
                    x, y = m(lon, lat)

                e = Ellipse([x, y],
                            phimax * ellipse_size_factor,
                            phimin * ellipse_size_factor,
                            az, **kwargs)
                ax.add_artist(e)
            # end if
        # end for
    # end func
# end class

def main():
    """
    define main function
    :return:
    """

    imaging2 = os.path.dirname(__file__)
    mtpy = os.path.dirname(imaging2)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    ModEM_files = os.path.join(data, 'ModEM_files')
    data_fn = os.path.join(ModEM_files, 'ModEM_Data_im2.dat')

    pmp = ModEM_ptensors(data_fn=data_fn)

    return
# end


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    main()