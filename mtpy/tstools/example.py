import sys

from PyQt5.QtWidgets import QWidget

from PyQt5.QtWidgets import QGraphicsView
from PyQt5.QtWidgets import QGraphicsScene
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QAbstractItemView
from PyQt5.QtWidgets import QTreeWidget
from PyQt5.QtWidgets import QTreeWidgetItem
from PyQt5.QtWidgets import QListView
from PyQt5.QtWidgets import QSplitter
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QStandardItemModel
from PyQt5.QtGui import QStandardItem

from tswaveitem import TSWaveItem

from tsscene import TSScene
from tsdata import TSData

import datetime


class TSWindow(QWidget):
    def __init__(self):
        super(TSWindow, self).__init__()

        # view
        self.scene = TSScene()
        viewWidget = QGraphicsView(self.scene)

        # control
        self.buttonOpenFile = QPushButton("open file")
        self.buttonOpenFile.clicked.connect(self.openfile)
        self.waveTree = QTreeWidget()
        self.waveTree.itemClicked.connect(self.showwave)

        controlLayout = QVBoxLayout()
        controlLayout.addWidget(self.buttonOpenFile)
        controlLayout.addWidget(self.waveTree)

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

    def showwave(self, wave):
        if wave.channelitem is None:
            print("empty channelitem")
            return
        else:
            self.scene.togglewave(wave)




    def exportwave(self, wave):
        print("export",wave.wavename)
        fname = QFileDialog.getSaveFileName(self, 'Save to','/g/data1a/ge3/yuhang/code/mtpy/mtpy/tstools')
        self.scene.exportwaveform(wave.wavename, fname[0])

    # set up wave tree in control region
    def setlist(self):
        item = self.waveTree.invisibleRootItem()
        self.fillitem(item, self.scene.getlist())
        self.waveTree.setSelectionMode(QAbstractItemView.MultiSelection)
        self.waveTree.show()

    # build wave tree
    def fillitem(self, item, value, parent=None):
        item.setExpanded(False)
        if type(value) is dict:
            for key, val in sorted(value.items()):
                child = TSWaveItem()
                child.setText(0, str(key))
                item.addChild(child)
                self.fillitem(child, val, str(key))
        elif type(value) is list:
            for val in value:
                child = TSWaveItem()
                item.addChild(child)
                if type(val) is dict:
                    child.setText(0, '[dict]')
                    self.fillitem(child, val, str(key))
                elif type(val) is list:
                    child.setText(0, '[list]')
                    self.fillitem(child, val, str(key))
                else:
                    print(type(val), "!!here")
                    child.setText(0, str(val))
                    child.channelitem = val
                    child.wavename = parent
                child.setExpanded(True)
        else:
            child = TSWaveItem()
            child.setText(0, str(value))
            child.channelitem = value
            child.wavename = parent
            item.addChild(child)

    def openfile(self):
        fname = QFileDialog.getOpenFileName(self,
                                            'Open file',
                                            '/g/data/ha3/Passive/_AusArray/OA/ASDF_BU/OA.h5', 'asdf file (*.h5)')

        if len(fname[0]) > 0:
            self.scene.setdata(fname[0])
            self.setlist()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = TSWindow()
    widget.resize(1024, 768)
    widget.show()
    sys.exit(app.exec_())
