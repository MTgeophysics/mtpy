# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 22:12:37 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""

import sys
from PyQt5 import QtCore, QtWidgets


class GetStations(QtWidgets.QDialog):
    def __init__(self, stations=None, parent=None):
        self.stations = stations
        self.checked_stations = []
        if self.stations is None:
            self.stations = ["mt01", "m102", "mt03"]
        super().__init__(parent)

        self.list_widget = QtWidgets.QListWidget()

        self.setup_ui()

    def setup_ui(self):
        for ss in self.stations:
            item = QtWidgets.QListWidgetItem(ss)
            # could be Qt.Unchecked; setting it makes the check appear
            item.setCheckState(QtCore.Qt.Unchecked)
            self.list_widget.addItem(item)

        done_button = QtWidgets.QPushButton("Done")
        done_button.clicked.connect(self.get_stations)

        cancel_button = QtWidgets.QPushButton("Cancel")
        cancel_button.clicked.connect(self.close)

        h_layout = QtWidgets.QHBoxLayout()
        h_layout.addWidget(self.list_widget, 1)

        buttons_layout = QtWidgets.QHBoxLayout()
        buttons_layout.addStretch(1)
        buttons_layout.addWidget(done_button)
        buttons_layout.addWidget(cancel_button)

        main_layout = QtWidgets.QVBoxLayout()
        main_layout.addLayout(h_layout)
        main_layout.addSpacing(12)
        main_layout.addLayout(buttons_layout)

        self.setLayout(main_layout)
        self.setWindowTitle("Choose Stations")
        self.show()

    def get_stations(self):
        self.checked_stations = []
        for item in self.list_widget.findItems("", QtCore.Qt.MatchContains):
            if item.checkState() > 0:
                self.checked_stations.append(item.text())
        self.close()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    dialog = GetStations()
    sys.exit(app.exec_())
