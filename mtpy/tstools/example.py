import sys

from PyQt5.QtWidgets import QWidget

from PyQt5.QtWidgets import QGraphicsView
from PyQt5.QtWidgets import QGraphicsScene
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QTreeWidget
from PyQt5.QtWidgets import QTreeWidgetItem
from PyQt5.QtWidgets import QListView
from PyQt5.QtWidgets import QSplitter
from PyQt5.QtWidgets import QDateTimeEdit
from PyQt5.QtCore import QDateTime
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QLineEdit

from PyQt5.QtGui import QStandardItemModel
from PyQt5.QtGui import QStandardItem
from PyQt5.QtWidgets import QCheckBox
from PyQt5 import QtCore

from tsscene import TSScene
from tsdata import TSData
from tswavetree import TSWaveTree


import datetime


class TSWindow(QWidget):
    def __init__(self):
        super(TSWindow, self).__init__()

        # time edit
        startlabel = QLabel('Start time')
        self.starttime = QLineEdit()
        endlabel = QLabel('End time')
        self.endtime = QLineEdit()
        buttonApply = QPushButton("Apply")
        buttonApply.clicked.connect(self.applytime)

        timeLayout = QHBoxLayout()
        timeLayout.addWidget(startlabel)
        timeLayout.addWidget(self.starttime)
        timeLayout.addWidget(endlabel)
        timeLayout.addWidget(self.endtime)
        timeLayout.addWidget(buttonApply)
        timeWidget = QWidget()
        timeWidget.setLayout(timeLayout)

        # view
        self.scene = TSScene(self)
        self.scene.starttimechanged.connect(self.starttime.setText)
        self.scene.endtimechanged.connect(self.endtime.setText)

        # gap mark
        gapmarkcheckbox = QCheckBox('mark gaps')
        gapmarkcheckbox.stateChanged.connect(self.scene.togglegap)

        viewLayout = QVBoxLayout()
        viewLayout.addWidget(timeWidget)
        viewLayout.addWidget(gapmarkcheckbox)
        viewLayout.addWidget(QGraphicsView(self.scene))
        viewWidget = QWidget()
        viewWidget.setLayout(viewLayout)

        # control
        self.buttonOpenFile = QPushButton("open file")
        self.buttonOpenFile.clicked.connect(self.openfile)
        self.waveTree = TSWaveTree()
        self.waveTree.header().hide()
        self.waveTree.itemClicked.connect(self.clickwave)
        self.waveTree.viewsegments.connect(self.scene.getsegments)
        self.buttonExportMeta = QPushButton("Export Meta Data")
        self.buttonExportMeta.clicked.connect(self.exportmeta)
        self.buttonExportWave = QPushButton("Export Waveforms")
        self.buttonExportWave.clicked.connect(self.exportwave)


        controlLayout = QVBoxLayout()
        controlLayout.addWidget(self.buttonOpenFile)
        controlLayout.addWidget(self.waveTree)
        controlLayout.addWidget(self.buttonExportMeta)
        controlLayout.addWidget(self.buttonExportWave)

        controlWidget = QWidget()
        controlWidget.setLayout(controlLayout)

        # put together
        split = QSplitter()
        split.addWidget(viewWidget)
        split.addWidget(controlWidget)

        layout = QHBoxLayout()
        layout.addWidget(split)
        self.setLayout(layout)
        self.setWindowTitle("TSView")

    def applytime(self):
        self.scene.applytime(self.starttime.text(), self.endtime.text())
        return

    def clickwave(self, wave: QTreeWidgetItem):
        if wave.childCount()==0:
            self.scene.togglewave(wave.text(0))








    def exportmeta(self):
        fname = QFileDialog.getSaveFileName(self,
                                            'Save as',
                                            '/g/data1a/ge3/yuhang/tmp', 'Text files (*.txt)')
        if len(fname[0]) > 0:
            self.scene.exportmetadata(fname)

    def openfile(self):
        fname = QFileDialog.getOpenFileName(self,
                                            'Open file',
                                            '/g/data/ha3/Passive/_AusArray/OA/ASDF_BU/OA.h5', 'asdf file (*.h5)')
                                            #'/g/data/ha3/rakib/ausLAMP/Data/Output/fixed/au.vic.h5', 'asdf file (*.h5)')
        if len(fname[0]) > 0:
            selecteditems = [i.text(0) for i in self.waveTree.selectedItems()]
            self.scene.loadfile(fname[0])
            self.waveTree.settree(self.scene.getlist(), selecteditems)



    def exportwave(self):
        fname = QFileDialog.getSaveFileName(self,
                                            'Save as',
                                            '/g/data1a/ge3/yuhang/tmp', 'MiniSEED (*.MSEED);; Text files (*.txt)')
        if len(fname[0]) > 0:
            self.scene.exportwaveform(fname)




if __name__ == "__main__":
    app = QApplication(sys.argv)

    widget = TSWindow()
    widget.resize(1680, 1050)
    widget.show()


    sys.exit(app.exec_())
