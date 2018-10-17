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

class TSWindow(QWidget):
    def __init__(self):
        super(TSWindow,self).__init__()

        self.scene = TSScene()
        self.view = QGraphicsView(self.scene)


        self.filebutton = QPushButton("open file")
        self.filebutton.clicked.connect(self.openfile)




        layout = QHBoxLayout()
        split = QSplitter()
        split.addWidget(self.view)

        buttonlayout = QVBoxLayout()
        buttonlayout.addWidget(self.filebutton)

        # self.wavelist = QListView()
        # self.wavelist.clicked.connect(self.showWave)
        # self.wavelistmodel = QStandardItemModel()
        # buttonlayout.addWidget(self.wavelist)

        self.wavetree = QTreeWidget()
        self.wavetree.itemClicked.connect(self.showWave)

        buttonlayout.addWidget(self.wavetree)
        buttonregion = QWidget()
        buttonregion.setLayout(buttonlayout)

        split.addWidget(buttonregion)
        layout.addWidget(split)

        self.setLayout(layout)
        self.setWindowTitle("TSView")

    def showWave(self, wave):
        print("here",wave.wavename)
        if wave.channelitem is None:
            print("empty channelitem")
            return
        else:
            print("showwave")
            self.scene.togglewave(wave)


            print(wave.channelitem.start_date, wave.channelitem.end_date)
        #wavename = self.wavelistmodel.itemFromIndex(index).text()
        #self.scene.togglewave(wavename, index.row())




    def setList(self):
        item = self.wavetree.invisibleRootItem()
        self.fill_item(item, self.scene.getList())
        self.wavetree.setSelectionMode(QAbstractItemView.MultiSelection)
        self.wavetree.show()

    def fill_item(self, item, value, parent=None):
        item.setExpanded(False)
        if type(value) is dict:
            for key, val in sorted(value.items()):
                child = TSWaveItem()
                child.setText(0, str(key))
                item.addChild(child)
                self.fill_item(child, val, str(key))
        elif type(value) is list:
            for val in value:
                child = TSWaveItem()
                item.addChild(child)
                if type(val) is dict:
                    child.setText(0, '[dict]')
                    self.fill_item(child, val, str(key))
                elif type(val) is list:
                    child.setText(0, '[list]')
                    self.fill_item(child, val, str(key))
                else:
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
        fname = QFileDialog.getOpenFileName(self, 'Open file','/g/data/ha3/Passive/_AusArray/OA/ASDF_BU/OA.h5','asdf file (*.h5)')

        if len(fname[0])>0:
            self.scene.setdata(fname[0])
            self.setList()





if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = TSWindow()
    widget.resize(1024, 768)
    widget.show()
    sys.exit(app.exec_())




