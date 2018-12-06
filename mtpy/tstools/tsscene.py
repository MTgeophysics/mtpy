
from PyQt5.QtWidgets import QGraphicsScene

from PyQt5.QtGui import QMouseEvent
from PyQt5.QtGui import QPen
from PyQt5 import QtGui

from PyQt5.QtCore import QTimer
from PyQt5.QtCore import Qt
from PyQt5 import QtCore

from PyQt5 import QtWidgets

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
from multiprocessing import Queue

import numpy as np




class TSScene(QGraphicsScene):

    starttimechanged = pyqtSignal(str)
    endtimechanged = pyqtSignal(str)

    def __init__(self, width=14, height=12, numofchannel=6):
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

        self.states = {'ready':0, 'busy':1}

        self.state = self.states['ready']
        self.installEventFilter(self)
        self.showgap = False

    def togglegap(self):
        self.showgap = ~self.showgap

        tmplist = self.visibleWave.copy()
        for wave in tmplist:
            self.togglewave(wave)
            self.togglewave(wave, tmplist[wave][1])

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

    def togglewave(self, wave: str, colorcode:int=0):
        if wave in self.visibleWave:
            axes, lines, _, _, _, _ = self.visibleWave[wave]
            self.removewave(axes, lines)
            self.visibleWave.pop(wave, None)
            channelid = self.axes.index(axes)
            self.axesavailability[channelid] = True
        else:
            # print(wave)
            waveform, wavename, starttime, endtime, gaps = self.data.getwaveform(wave, self.starttime, self.endtime)
            axes, lines = self.displaywave(wavename, waveform, gaps)
            if axes is not None:
                self.visibleWave[wave] = (axes, lines, colorcode, starttime, endtime, gaps)
                #print("togglewave:", starttime, endtime)


    def displaywave(self, wavename: str, waveform: np.array, gaps, colorcode: int=None):
        # print(gaps)
        if True not in self.axesavailability:
            return None, None
        else:
            # print(waveform.shape,'!'*10,'displaywave')
            location = self.axesavailability.index(True)
            axes = self.axes[location]
            self.axesavailability[location] = False
            if wavename is not None and waveform is not None:
                if colorcode is None:
                    colorcode = 'C'+str(location%10)
                # print(waveform.shape,'='*8)
                times = waveform[0,:]
                span = round(len(times)/4)
                print(UTCDateTime(times[0]),UTCDateTime(times[-1]),'out')
                # print(span)
                # if span<1:
                #     span = 1
                axes.set_xticks(times[::span])
                #axes.set_xticklabels([datetime(int(t)).strftime("%Y-%m-%d %H:%M:%S") for t in times[::span]])
                axes.set_xticklabels([UTCDateTime(t).strftime("%Y-%m-%d %H:%M:%S") for t in times[::span]])
                #print([UTCDateTime(t).strftime("%Y-%m-%d %H:%M:%S") for t in times[::span]])
                #print([UTCDateTime(t) for t in times[::span]])
                #print(times[::span])
                lines = axes.plot(times, waveform[1,:],linestyle="-", label=wavename, color=colorcode)
                #lines = axes.plot(range(len(times)),times, linestyle="-", label=wavename, color=colorcode)
                if self.showgap:
                    for g in gaps:
                        #print(g[4],g[5])
                        if g[4].timestamp>=times[0] and g[5].timestamp<times[-1]:
                            axes.axvspan(g[4],g[5],facecolor='0.2',alpha=0.5)
                axes.legend()
                self.downx = None

                self.canvas.draw()

                if self.endtime is not None and self.starttime is not None:
                    timewindow = self.endtime-self.starttime
                    if abs(times[0]-times[-1]-timewindow)/timewindow<0.1:
                        self.starttime = UTCDateTime(times[0])
                        self.endtime = self.starttime + timewindow
                else:
                    self.starttime = UTCDateTime(times[0])
                    self.endtime = UTCDateTime(times[-1])



                self.starttimechanged.emit(self.starttime.strftime("%Y-%m-%d %H:%M:%S"))
                self.endtimechanged.emit(self.endtime.strftime("%Y-%m-%d %H:%M:%S"))
                return axes, lines
            else:
                lines = None
                axes.legend([wavename])

            return axes, lines

    def displaywave_(self, wavename: str, waveform: Trace, colorcode: int=None):
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
        #axes.get_legend().remove()
        axes.clear()
        self.canvas.draw()

    def timeshift(self, shift: float):
        # print("shift")
        self.state = False

        if self.starttime is None:
            return

        #print('shift=',self.starttime- self.endtime)

        shift = (self.endtime-self.starttime)*shift

        starttime = self.starttime + shift
        endtime = self.endtime + shift

        #print('+'*10,self.endtime- self.starttime )

        #print(starttime, endtime, self.starttime, self.endtime)

        for wave in self.visibleWave:
            if starttime<self.visibleWave[wave][3]:
                starttime = self.starttime
            if endtime>self.visibleWave[wave][4]:
                endtime = self.endtime

        #print(starttime, endtime, self.starttime, self.endtime,'!!!!!!')


        if starttime!=self.starttime or endtime!=self.endtime:
            # print('update'*10)
            self.starttime = starttime
            self.endtime = endtime

            tmplist = self.visibleWave.copy()

            for wave in tmplist:
                self.togglewave(wave)
                self.togglewave(wave, tmplist[wave][2])



        return

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




    def mousePressEvent(self, event: QMouseEvent):
        super(TSScene, self).mousePressEvent(event)
        if self.starttime is None:
            return
        self.downx = event.scenePos().x()
        self.downbutton = event.button()
        self.state = self.states['ready']


    def performshift(self):
        if self.state == self.states['busy']:
            shift = float(self.downx - self.currentx) / self.graphwidth
            self.timeshift(shift)
            self.downx = self.currentx
            self.state = self.states['ready']


    def mouseMoveEvent(self, event: QMouseEvent):
        #print(self.state)
        if self.starttime is None:
            return
        if self.downx is not None:

            if self.downbutton == Qt.LeftButton and self.state==self.states['ready']:
                self.state = self.states['busy']
                QTimer.singleShot(0, self.performshift)
            elif self.downbutton == Qt.LeftButton:
                # print('ignored=============================')
                pass


            elif self.downbutton == Qt.RightButton:
                self.removeItem(self.rect)
                self.rect = self.addRect(self.downx,0, event.scenePos().x()-self.downx, self.height(), pen=QPen(Qt.red))

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
        # self.downx = None
        self.downbutton = None
        self.removeItem(self.rect)
        self.rect = None

    def wheelEvent(self, event: QMouseEvent):
        super(TSScene, self).wheelEvent(event)

        delta = -event.delta() / 8 / 15

        if self.wheelactive==False:
            self.wheelactive = True
            self.timescale(delta)
            self.wheelactive = False




    def exportwaveform(self, filename: tuple):
        traces = []
        for wave in self.visibleWave:
            fill_value = 'last'
            waveform, wavename, starttime, endtime, gaps = self.data.readdisc(wave, self.starttime, self.endtime, resample=False, fill_value=fill_value)
            traces.append(waveform)

        stream = Stream(traces=traces)
        if 'MSEED' in filename[1]:
            stream.write(filename[0] + ".mseed", format='MSEED')
        elif 'txt' in filename[1]:
            stream.write(filename[0] + ".txt", format='TSPAIR')


    def gettimeboundary(self):
        return self.starttime, self.endtime


    def eventFilter(self, source, event):
        #print(event.type())
        if event.type() == 155:
            #print(event)
            if event.buttons() == QtCore.Qt.LeftButton:
                self.currentx = event.scenePos().x()
            else:
                pass # do other stuff

        return False