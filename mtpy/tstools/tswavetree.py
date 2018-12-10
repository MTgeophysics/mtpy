from PyQt5.QtWidgets import QTreeWidget
from PyQt5.QtWidgets import QTreeWidgetItem
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QAbstractItemView
from PyQt5.QtWidgets import QMenu
from PyQt5 import QtGui
from PyQt5 import QtCore

class TSWaveTree(QTreeWidget):
    def __init__(self):
        super(TSWaveTree, self).__init__()
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.custommenu)


    # set up wave tree in control region
    def settree(self, wavelist):
        item = self.invisibleRootItem()

        for c in reversed(range(item.childCount())):
            item.removeChild(item.child(c))

        self.fillitem(item, wavelist)
        self.setSelectionMode(QAbstractItemView.MultiSelection)
        self.show()

    # build wave tree
    def fillitem(self, node: QTreeWidgetItem, value: object):
        node.setExpanded(False)
        if type(value) is dict:
            for key, val in sorted(value.items()):
                child = QTreeWidgetItem()
                child.setText(0, str(key))
                node.addChild(child)
                child.setFlags(child.flags() & ~Qt.ItemIsSelectable)
                self.fillitem(child, val)
        elif type(value) is list:
            for idx, val in enumerate(value):
                child = QTreeWidgetItem()
                child.setText(0, val)
                # child.setFlags(Qt.ItemIsSelectable | Qt.ItemIsUserCheckable)
                node.addChild(child)


    def custommenu(self, pos):
        currentitem = self.itemAt(pos)
        if currentitem.childCount()==0:
            menu = QMenu(self)
            actionviewwave = menu.addAction("view/hide waveform")
            actionviewsegment = menu.addAction("view segments")
            action = menu.exec_(self.mapToGlobal(pos))

            if action == actionviewwave:
                self.itemClicked.emit(currentitem, self.currentColumn())
            elif action == actionviewsegment:
                print("view segments")



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
