from PyQt5.QtWidgets import QTreeWidgetItem


class TSWaveItem(QTreeWidgetItem):
    def __init__(self):
        super(TSWaveItem, self).__init__()

        self.channelitem = None
        self.wavename = None
        self.network = None
        self.station = None




