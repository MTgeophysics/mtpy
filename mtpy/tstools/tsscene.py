
from PyQt5.QtWidgets import QGraphicsScene

from PyQt5.QtGui import QMouseEvent
from PyQt5.QtGui import QPen

from PyQt5.QtCore import Qt

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qt5agg import FigureCanvas
import matplotlib.pyplot as plt

from tsdata import TSData
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from PyQt5.QtCore import pyqtSignal

from obspy.core.utcdatetime import UTCDateTime
from datetime import datetime

class TSScene(QGraphicsScene):

    starttimechanged = pyqtSignal(str)
    endtimechanged = pyqtSignal(str)

    def __init__(self, width=14, height=9, numofchannel=4):
        super(TSScene, self).__init__()

        # set waveform windows
        figure = Figure()
        figure.set_size_inches(width, height)
        self.graphwidth = figure.dpi * width
        self.canvas = FigureCanvas(figure)
        self.addWidget(self.canvas)
        self.axesavailability = [True for i in range(numofchannel)]
        self.axes = []
        for i in range(numofchannel):
            self.axes.append(figure.add_subplot(str(numofchannel)+'1'+str(i+1)))

        # set backend data model
        self.data = None
        self.visibleWave = {}
        self.starttime = None
        self.endtime = None

        # prepare for user input
        self.downx = None
        self.wheelactive = False
        self.rect = None

    def applytime(self, start: str, end: str):
        if self.data is None:
            return
        self.starttime = UTCDateTime(start)
        self.endtime = UTCDateTime(end)

        tmplist = self.visibleWave.copy()
        for wave in tmplist:
            self.togglewave(wave)
            self.togglewave(wave, tmplist[wave][2])

    def setdata(self, filename: str):
        self.data = TSData(filename)

    def getlist(self):
        return self.data.getlist()

    def togglewave(self, wave: str, colorcode:int=0, samplerate: int=1000):
        if wave in self.visibleWave:
            axes = self.visibleWave[wave][0]
            lines = self.visibleWave[wave][1]
            self.removewave(axes, lines)
            self.visibleWave.pop(wave, None)
            self.axesavailability[self.axes.index(axes)] = True

        else:
            waveform, wavename, starttime, endtime = self.data.getwaveform(wave, self.starttime, self.endtime,samplerate)
            axes, lines = self.displaywave(wavename, waveform)
            if axes is not None:
                self.visibleWave[wave] = (axes, lines, colorcode, starttime, endtime)

    def displaywave(self, wavename: str, waveform: Trace, colorcode: int=None):
        if True not in self.axesavailability:
            return None, None
        else:
            location = self.axesavailability.index(True)
            axes = self.axes[location]
            self.axesavailability[location] = False
            if wavename is not None and waveform is not None:
                if colorcode is None:
                    colorcode = 'C'+str(location%10)

                times = [waveform.meta['starttime']+t for t in waveform.times()]
                span = round(len(times)/4)
                if span<1:
                    span = 1
                axes.set_xticks(times[::span])
                axes.set_xticklabels([t.strftime("%Y-%m-%d %H:%M:%S") for t in times[::span]])
                lines = axes.plot(times, waveform.data,linestyle="-", label=wavename, color=colorcode)
                axes.legend()
                self.downx = None

                self.canvas.draw()

                self.starttime = waveform.meta['starttime']
                self.endtime = waveform.meta['endtime']

                self.starttimechanged.emit(self.starttime.strftime("%Y-%m-%d %H:%M:%S"))
                self.endtimechanged.emit(self.endtime.strftime("%Y-%m-%d %H:%M:%S"))
                return axes, lines
            else:
                lines = None
                axes.legend([wavename])

            return axes, lines


    def removewave(self, axes: Axes, lines: Line2D):
        if lines is not None:
            lines.pop(0).remove()
        axes.relim()
        axes.autoscale_view(True, True, True)
        axes.legend()
        self.canvas.draw()

    def timeshift(self, shift: float):
        if self.starttime is None:
            return

        shift = (self.endtime-self.starttime)*shift

        starttime = self.starttime + shift
        endtime = self.endtime + shift

        for wave in self.visibleWave:
            if starttime<self.visibleWave[wave][3]:
                starttime = self.starttime
            if endtime>self.visibleWave[wave][4]:
                endtime = self.endtime

        if starttime!=self.starttime and endtime!=self.endtime:
            self.starttime = starttime
            self.endtime = endtime
            tmplist = self.visibleWave.copy()
            for wave in tmplist:
                self.togglewave(wave)
                self.togglewave(wave, tmplist[wave][2])


    def timescale(self, delta: float):
        if self.starttime is None:
            return

        shift = (self.endtime - self.starttime) * -delta*0.1

        starttime = self.starttime + shift
        endtime = self.endtime - shift


        for wave in self.visibleWave:
            if starttime<self.visibleWave[wave][3]:
                starttime = self.starttime
            if endtime>self.visibleWave[wave][4]:
                endtime = self.endtime


        if endtime-starttime<0.1:
            pass
        elif starttime==self.starttime and endtime==self.endtime:
            pass
        else:
            self.starttime = starttime
            self.endtime = endtime
            tmplist = self.visibleWave.copy()
            for wave in tmplist:
                self.togglewave(wave)
                self.togglewave(wave, tmplist[wave][1])

        self.wheelactive = False




    def mousePressEvent(self, event: QMouseEvent):
        super(TSScene, self).mousePressEvent(event)
        if self.starttime is None:
            return
        self.downx = event.scenePos().x()
        self.downbutton = event.button()


    def mouseMoveEvent(self, event: QMouseEvent):
        if self.starttime is None:
            return
        if self.downx is not None:
            if self.downbutton == Qt.LeftButton:
                self.upx = event.scenePos().x()
                shift = float(self.downx - self.upx) / self.graphwidth
                self.timeshift(shift)
                self.downx=self.upx
            elif self.downbutton == Qt.RightButton:
                self.removeItem(self.rect)
                self.rect = self.addRect(self.downx,0, event.screenPos().x()-self.downx, self.height(), pen=QPen(Qt.red))

    def mouseReleaseEvent(self, event: QMouseEvent):
        super(TSScene, self).mousePressEvent(event)
        if self.starttime is None:
            return
        if event.button() == Qt.RightButton:
            left = 225
            right = 1215
            start = self.starttime+(self.downx-left)/(right-left)*(self.endtime-self.starttime)
            end = self.starttime+(event.scenePos().x()-left)/(right-left)*(self.endtime-self.starttime)
            self.applytime(start, end)
        self.downx = None
        self.downbutton = None
        self.removeItem(self.rect)
        self.rect = None

    def wheelEvent(self, event: QMouseEvent):
        super(TSScene, self).wheelEvent(event)

        delta = -event.delta() / 8 / 15

        if self.wheelactive==False:
            self.wheelactive = True
            self.timescale(delta)




    def exportwaveform(self, filename: tuple):
        traces = []
        for wave in self.visibleWave:
            waveform, wavename, starttime, endtime = self.data.getwaveform(wave, self.starttime, self.endtime, 0)
            traces.append(waveform)

        stream = Stream(traces=traces)
        if 'MSEED' in filename[1]:
            stream.write(filename[0] + ".mseed", format='MSEED')
        elif 'txt' in filename[1]:
            stream.write(filename[0] + ".txt", format='TSPAIR')


    def gettimeboundary(self):
        return self.starttime, self.endtime