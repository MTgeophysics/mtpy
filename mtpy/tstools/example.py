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
from PyQt5.QtWidgets import QListView
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QStandardItemModel
from PyQt5.QtGui import QStandardItem


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
        layout.addWidget(self.view)

        buttonlayout = QVBoxLayout()
        buttonlayout.addWidget(self.filebutton)

        self.wavelist = QListView()
        self.wavelist.clicked.connect(self.showWave)
        #self.wavelist.selectionChanged.connect(self.showWave)
        self.wavelistmodel = QStandardItemModel()
        buttonlayout.addWidget(self.wavelist)







        layout.addLayout(buttonlayout)


        self.setLayout(layout)
        self.setWindowTitle("TSView")

    def showWave(self, index):
        wavename = self.wavelistmodel.itemFromIndex(index).text()
        self.scene.togglewave(wavename, index.row())




    def setList(self,wlist):
        for w in wlist:
            item = QStandardItem(w)
            item.setSelectable(True)
            self.wavelistmodel.appendRow(item)
            self.wavelist.setModel(self.wavelistmodel)
            self.wavelist.setSelectionMode(QAbstractItemView.MultiSelection)


    def openfile(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file','/g/data1/ha3/rakib/ausLAMP/Data/Output','asdf_file (*.h5)')

        if len(fname[0])>0:
            self.scene.setdata(fname[0])

            self.setList(self.scene.getList())





if __name__ == "__main__":
    app = QApplication(sys.argv)
    widget = TSWindow()
    widget.resize(1024, 768)
    widget.show()
    sys.exit(app.exec_())




