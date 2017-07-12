# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from mtpy.gui.SmartMT.gui.plot_parameter import FrequencySingle, Ellipse, FrequencyTolerance, ColorBar
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps


class PhaseTensorMap(VisualizationBase):
    @staticmethod
    def plot_description():
        return """
<p>Plot phase tensor map in Lat-Lon Coordinate System.</p>
        """

    def update_ui(self):
        self._frequency_ui.set_data(self._mt_objs)

    @staticmethod
    def plot_name():
        return "Phase Tensor Map"

    def plot(self):
        # get data
        # NOTE: this is a hack because the existing bug(s) in the PlotPhaseTensorMaps class, that is the
        # constructor of the class only populate all the necessary information correctly when reads from file
        # this is the only way before this bug(s) is fixed
        file_list = []
        for mt_obj in self._mt_objs:
            file_list.append(mt_obj.fn)

        params = {
            'fn_list': file_list,
            'plot_freq': self._frequency_ui.get_frequency(),
            'ftol': self._tolerance_ui.get_tolerance_in_float(),
            'ellipse_dict': self._ellipse_ui.get_ellipse_dict(),
            'mapscale': 'deg',  # deg or m, or km
            'xpad': 0.4,  # plot margin; change according to lat-lon in edifiles
            'ypad': 0.4,  # ~ 2* ellipse size
            'plot_tipper': 'yr',
            'arrow_dict': {
                'size': 0.5,
                'lw': 0.2,
                'head_width': 0.04,
                'head_length': 0.04,
                'threshold': 0.8,
                'direction': 0}
        }

        cb_dict = self._colorbar_ui.get_colorbar_dict()
        if cb_dict is not None:
            params['cb_dict'] = cb_dict

        self._plotting_object = PlotPhaseTensorMaps(**params)
        self._plotting_object.plot(show=False)
        self._fig = self._plotting_object.fig

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        # set up ui
        self._frequency_ui = FrequencySingle(self._parameter_ui)
        self._frequency_ui.setTitle("Frequency (Hz)")
        self._parameter_ui.add_parameter_groubox(self._frequency_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._tolerance_ui)

        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._ellipse_ui)

        self._colorbar_ui = ColorBar(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._colorbar_ui)

        # resize
        self._parameter_ui.resize(self._parameter_ui.width(),
                                  self._parameter_ui.sizeHint().height())

        self.update_ui()
