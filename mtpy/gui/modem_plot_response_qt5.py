# -*- coding: utf-8 -*-
"""
ModEM data and response visualization with a gui.

The user will be able to choose from stations within the data to look at
in either impedance or apparent resistivity and phase.

The functionality is quite simple at the moment

JP 2016
"""
#
# ==============================================================================
# Imports
# ==============================================================================
import sys
from pathlib import Path
import numpy as np

try:
    from PyQt5 import QtCore, QtWidgets
except ImportError:
    raise ImportError("This version needs PyQt5")

from mtpy.gui.modem_plot_response_gui import PlotResponses
from mtpy.gui.response_plot_settings import PlotSettings
from mtpy.gui.get_stations import GetStations
from mtpy.gui.plot_stations import PlotStations

# ==============================================================================


class ModEMPlotResponse(QtWidgets.QMainWindow):
    """
    main window
    """

    def __init__(self):
        super(ModEMPlotResponse, self).__init__()

        self.data_fn = None
        self.resp_fn = None
        self.modem_data = None

        self.setup_ui()

    def setup_ui(self):
        self.setWindowTitle("Plot ModEM Responses")
        self.setWindowState(QtCore.Qt.WindowMaximized)
        screen_shape = QtWidgets.QDesktopWidget().screenGeometry()

        # create a menu bar on the window with 4 different items
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, screen_shape.width(), 38))

        # add a tab for File --> open, close, save
        self.menu_data_file = QtWidgets.QMenu(self.menubar)
        self.menu_data_file.setTitle("Data File")

        self.menu_resp_file = QtWidgets.QMenu(self.menubar)
        self.menu_resp_file.setTitle("Response File")

        # add a tab for chaning the display
        self.menu_display = QtWidgets.QMenu(self.menubar)
        self.menu_display.setTitle("Display")

        # add a tab for help
        self.menu_help = QtWidgets.QMenu(self.menubar)
        self.menu_help.setTitle("Help")

        self.setMenuBar(self.menubar)

        # set the actions for the data file menu item
        # set an open option that on click opens a modem file
        self.data_action_open = QtWidgets.QAction(self)
        self.data_action_open.setText("&Open")
        self.data_action_open.setShortcut("Ctrl+o")
        self.data_action_open.triggered.connect(self.get_data_file)

        # set an open option that on click opens a modem file
        self.data_action_new = QtWidgets.QAction(self)
        self.data_action_new.setText("&New")
        self.data_action_new.setShortcut("Ctrl+n")
        self.data_action_new.triggered.connect(self.new_data_file)

        # set a close that closes the main window
        self.data_action_close = QtWidgets.QAction(self)
        self.data_action_close.setText("Close")
        self.data_action_close.setShortcut("Ctrl+x")
        self.data_action_close.triggered.connect(self.close)

        # set a save option that will eventually save the masked data
        self.data_action_save = QtWidgets.QAction(self)
        self.data_action_save.setText("&Save Edits")
        self.data_action_save.setShortcut("Ctrl+s")
        self.data_action_save.triggered.connect(self.save_edits)

        # add station(s)
        self.data_action_add = QtWidgets.QAction(self)
        self.data_action_add.setText("Add Stations")
        self.data_action_add.triggered.connect(self.add_station)

        # remove station(s)
        self.data_action_remove = QtWidgets.QAction(self)
        self.data_action_remove.setText("Remove Stations")
        self.data_action_remove.triggered.connect(self.remove_station)

        # add the action on the menu tab
        self.menu_data_file.addAction(self.data_action_open)
        self.menu_data_file.addAction(self.data_action_new)
        self.menu_data_file.addAction(self.data_action_add)
        self.menu_data_file.addAction(self.data_action_remove)
        self.menu_data_file.addAction(self.data_action_close)
        self.menu_data_file.addAction(self.data_action_save)
        self.menubar.addAction(self.menu_data_file.menuAction())

        # set the action items for the response file
        self.action_resp_open = QtWidgets.QAction(self)
        self.action_resp_open.setText("Open")
        self.action_resp_open.triggered.connect(self.get_resp_fn)
        self.menu_resp_file.addAction(self.action_resp_open)
        self.menubar.addAction(self.menu_resp_file.menuAction())
        #
        # adding options for display plot type
        self.menu_plot_type = QtWidgets.QMenu(self)
        self.menu_plot_type.setTitle("Plot Type")
        self.menu_display.addMenu(self.menu_plot_type)
        self.menubar.addAction(self.menu_display.menuAction())

        # set plot impedance or resistivity and phase
        self.action_plot_z = QtWidgets.QAction(self)
        self.action_plot_z.setText("Impedance")
        self.action_plot_z.setCheckable(True)
        self.menu_plot_type.addAction(self.action_plot_z)
        self.action_plot_z.toggled.connect(self.status_checked_ptz)

        self.action_plot_rp = QtWidgets.QAction(self)
        self.action_plot_rp.setText("Resistivity-Phase")
        self.action_plot_rp.setCheckable(True)
        self.menu_plot_type.addAction(self.action_plot_rp)
        self.action_plot_rp.toggled.connect(self.status_checked_ptrp)

        self.action_plot_settings = QtWidgets.QAction(self)
        self.action_plot_settings.setText("Settings")
        self.action_plot_settings.triggered.connect(self.show_settings)
        self.menu_display.addAction(self.action_plot_settings)
        self.menubar.addAction(self.menu_display.menuAction())

        self.menu_display.addAction(self.menu_plot_type.menuAction())

        self.action_help = QtWidgets.QAction(self)
        self.action_help.setText("Help Documentation")
        self.action_help.triggered.connect(self.disp_help)
        self.menu_help.addAction(self.action_help)
        self.menubar.addAction(self.menu_help.menuAction())

        self.plot_response = PlotResponses(self.data_fn, self.resp_fn)
        self.setCentralWidget(self.plot_response)

        # self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(self)

    def status_checked_ptz(self, toggled):
        """
        be sure that only one plot style is checked
        """

        self.plot_response.plot_z = toggled
        if toggled == True:
            untoggled = False

        elif toggled == False:
            untoggled = True

        self.action_plot_z.setChecked(toggled)
        self.action_plot_rp.setChecked(untoggled)

    def status_checked_ptrp(self, toggled):
        """
        be sure that only one plot style is checked
        """

        if toggled == True:
            untoggled = False
            self.plot_response.plot_z = False
        elif toggled == False:
            untoggled = True
            self.plot_response.plot_z = True

        self.action_plot_z.setChecked(untoggled)
        self.action_plot_rp.setChecked(toggled)

    def get_data_file(self):
        """
        get the filename from a file dialogue

        """

        fn_dialog = QtWidgets.QFileDialog()
        fn = Path(
            str(
                fn_dialog.getOpenFileName(
                    caption="Choose ModEM data file", filter="(*.dat);; (*.data)"
                )[0]
            )
        )

        self.plot_response.data_fn = fn
        self.dir_path = fn.parent

        self.station_plot = PlotStations(
            self.plot_response.modem_data.station_locations
        )
        self.station_plot.plot()

        self.station_plot.show()
        self.station_plot.stationChanged.connect(self.station_picked)

        self.plot_response.list_widget.currentItemChanged.connect(
            self.update_station_map
        )

    def update_station_map(self, widget_item):
        self.station_plot.previous_index = int(self.station_plot.current_index)
        self.station_plot.current_index = int(
            np.where(
                self.plot_response.modem_data.station_locations.station
                == self.plot_response.station
            )[0][0]
        )
        self.station_plot.plot_new_station()

    def station_picked(self):
        self.plot_response.station = self.station_plot.current_station
        self.plot_response.plot()

    def save_edits(self):
        """
        save edits to another file
        """
        fn_dialog = QtWidgets.QFileDialog()
        save_fn = Path(
            str(
                fn_dialog.getSaveFileName(
                    caption="Choose File to save", filter="*.dat"
                )[0]
            )
        )

        self.plot_response.modem_data.write_data_file(
            save_path=save_fn.parent,
            fn_basename=save_fn.name,
            compute_error=False,
            fill=True,
            elevation=self.plot_response.modem_data.topography,
        )

    def new_data_file(self):
        """build a new data file from scratch"""
        pass

    def add_station(self):
        """
        Add a station or list of stations from files
        """
        extensions = "EDI (*.edi);;EMTFXML (*.xml);;ZMM (*.zmm);;J (*.j)"
        fn_dialog = QtWidgets.QFileDialog()
        fn_names = fn_dialog.getOpenFileNames(
            caption="Choose ModEM data file", filter=extensions
        )

        fn_list = []
        for ii in range(0, len(fn_names), 2):
            fn_list += fn_names[ii]
        fn_list = [Path(fn) for fn in fn_list]

        new_array, new_dict = self.plot_response.modem_data.add_station(fn_list)
        self.plot_response.modem_data.data_array = new_array
        self.plot_response.modem_data.mt_dict = new_dict

        # fill list of stations
        station_list = list(sorted(self.plot_response.modem_data.mt_dict.keys()))
        self.plot_response.list_widget.clear()
        for station in station_list:
            self.plot_response.list_widget.addItem(station)

        if self.plot_response.station is None:
            self.plot_response.station = station_list[0]

        self.plot_response.plot()

        self.station_plot.redraw_plot(self.plot_response.modem_data.station_locations)

    def remove_station(self):
        """
        Remove stations.

        :return: DESCRIPTION
        :rtype: TYPE

        """

        rs = GetStations(
            stations=list(self.plot_response.modem_data.station_locations.station)
        )
        rs.exec_()

        new_data, new_mtdict = self.plot_response.modem_data.remove_station(
            rs.checked_stations
        )
        self.plot_response.modem_data.data_array = new_data
        self.plot_response.modem_data.mt_dict = new_mtdict

        # fill list of stations
        station_list = list(sorted(self.plot_response.modem_data.mt_dict.keys()))
        self.plot_response.list_widget.clear()
        for station in station_list:
            self.plot_response.list_widget.addItem(station)

        if self.plot_response.station not in station_list:
            self.plot_response.station = station_list[0]

        self.plot_response.plot()

        self.station_plot.redraw_plot(self.plot_response.modem_data.station_locations)

    def get_resp_fn(self):
        """
        get response file name
        """

        fn_dialog = QtWidgets.QFileDialog(directory=self.dir_path.as_posix())
        fn = Path(
            str(
                fn_dialog.getOpenFileName(
                    caption="Choose ModEM response file", filter="(*.dat);; (*.data)"
                )[0]
            )
        )

        self.plot_response.resp_fn = fn

    def show_settings(self):
        self.settings_window = PlotSettings(**self.__dict__)
        self.settings_window.show()
        self.settings_window.settings_updated.connect(self.update_settings)

    def update_settings(self):

        for attr in sorted(self.settings_window.__dict__.keys()):
            setattr(self, attr, self.settings_window.__dict__[attr])

        self.plot()

    def disp_help(self):
        """
        display a help dialogue
        """
        ll = [
            "This GUI will allow you to edit your data by masking points",
            "and adding error bars to dodgy data points.  Only the top row",
            "row is editable for now. However, all edits to the top row ",
            "are applied to the bottom row (real and imaginary parts).\n",
            "   * Left-Click the mouse to mask a point this will mask both",
            "     the real and imaginary part of that component.\n",
            "   * Right-Click the mouse to add error bars to the data point",
            "     again this will apply to both real and imaginary parts of",
            "     the selected component. Current it goes up by 5%\n",
            "   * To save your masking, go to Data File -> Save Edits",
        ]

        help_string = "\n".join(ll)

        QtWidgets.QMessageBox.information(self.centralWidget, "Help", help_string)


# ==============================================================================
# Def Main
# ==============================================================================
def main():
    app = QtWidgets.QApplication(sys.argv)
    ui = ModEMPlotResponse()
    ui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
