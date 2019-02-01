from PyQt5.QtWidgets import QTreeWidget
from PyQt5.QtWidgets import QTreeWidgetItem
from PyQt5.QtCore import Qt
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QAbstractItemView
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QListWidget
from PyQt5.QtWidgets import QListView
from PyQt5.QtWidgets import QDialog
from PyQt5.QtWidgets import QMenu
from PyQt5.QtWidgets import QVBoxLayout

class TSWaveTree(QTreeWidget):
    viewsegments = pyqtSignal(str, object)

    def __init__(self):
        super(TSWaveTree, self).__init__()
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.custommenu)


    # set up wave tree in control region
    def settree(self, wavelist, selecteditems):
        item = self.invisibleRootItem()

        for c in reversed(range(item.childCount())):
            item.removeChild(item.child(c))

        self.fillitem(item, wavelist, selecteditems)
        self.setSelectionMode(QAbstractItemView.MultiSelection)
        self.show()

    # build wave tree
    def fillitem(self, node: QTreeWidgetItem, value: object, selecteditems):
        node.setExpanded(False)
        if type(value) is dict:
            for key, val in sorted(value.items()):
                child = QTreeWidgetItem()
                child.setText(0, str(key))
                node.addChild(child)
                child.setFlags(child.flags() & ~Qt.ItemIsSelectable)
                self.fillitem(child, val, selecteditems)
        elif type(value) is list:
            for idx, val in enumerate(value):
                child = QTreeWidgetItem()
                child.setText(0, val)
                # child.setFlags(Qt.ItemIsSelectable | Qt.ItemIsUserCheckable)
                node.addChild(child)
                if val in selecteditems:
                    child.setSelected(True)


    def custommenu(self, pos):
        currentitem = self.itemAt(pos)
        if currentitem.childCount()==0:
            menu = QMenu(self)
            actionviewwave = menu.addAction("view/hide waveform")
            actionviewsegment = menu.addAction("view segments")
            action = menu.exec_(self.mapToGlobal(pos))

            if action == actionviewwave:
                self.itemClicked.emit(currentitem, self.currentColumn())
                currentitem.setSelected(not currentitem.isSelected())
            elif action == actionviewsegment:
                print("view segments")
                self.viewsegments.emit(currentitem.text(0), currentitem)
                #
                #
                # wavelist = QListWidget()
                # wavelist.addItem("123")
                # wavelist.addItem("456")
                # wavelist.addItem("789")
                #
                # wavelistwindowlayout = QVBoxLayout()
                # wavelistwindowlayout.addWidget(wavelist)
                #
                # wavelistwindow = QDialog(self)
                # wavelistwindow.setWindowTitle('segments in '+currentitem.parent().text(0))
                # wavelistwindow.setLayout(wavelistwindowlayout)
                # wavelistwindow.show()



    # def contextMenuEvent(self, event: QtGui.QContextMenuEvent):
    #     menu  = QMenu()
    #     listaction = menu.addAction("view segments")
    #     action = menu.exec_(self.mapToGlobal(event.pos()))
    #     if action == listaction:
    #         print("test successful")
    #
    #
    # def eventFilter(self, source: QtCore.QObject, event: QtCore.QEvent):
    #     if event.type() == QtCore.QEvent.ContextMenu:
    #         print("here")
    #         return False
    #     else:
    #         print(source)
    #         print(event)
    #         return super(TSWaveTree, self).eventFilter(source, event)
