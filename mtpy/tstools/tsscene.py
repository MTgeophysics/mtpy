
from PyQt5.QtWidgets import QGraphicsScene

from PyQt5.QtGui import QPixmap

from PyQt5.QtCore import Qt

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas
import matplotlib.pyplot as plt

from tsdata import TSData


class TSScene(QGraphicsScene):
    def __init__(self, width=14, height=12, numofchannel=4):
        super(TSScene, self).__init__()



        self.pixmap = QPixmap()

        figure = Figure()
        figure.set_size_inches(width, height)
        figure.tight_layout()
        self.canvas = FigureCanvas(figure)
        self.plothandle=self.addWidget(self.canvas)

        self.graphwidth = figure.dpi * width


        self.line = None

        self.downx = None
        self.data = None

        self.axes = []
        for i in range(numofchannel):
            self.axes.append(figure.add_subplot(str(numofchannel)+'1'+str(i+1)))

        self.axesavailability = [True for i in range(numofchannel)]


        self.visibleWave = {}

        self.starttime = None
        self.endtime = None

        self.wheelactive = False

    def togglewave(self, wave, colorcode=0):
        if wave in self.visibleWave:
            axes = self.visibleWave[wave][0]
            handle = self.visibleWave[wave][1]
            self.removewave(axes, handle)
            self.visibleWave.pop(wave, None)
            self.axesavailability[self.axes.index(axes)] = True

        else:
            stream, wavename, starttime, endtime = self.data.getwaveform(wave, self.starttime, self.endtime)
            waveform = stream[0]
            axes, handle = self.displaywave(wavename, waveform, colorcode)
            self.visibleWave[wave] = (axes, handle, colorcode, starttime, endtime)

    def displaywave(self, wavename, waveform, colorcode):
        if True not in self.axesavailability:
            pass
        else:
            location = self.axesavailability.index(True)
            axes = self.axes[location]
            self.axesavailability[location] = False

            colorcode = 'C'+str(colorcode%10)

            times = [waveform.meta['starttime']+t for t in waveform.times()]
            handle = axes.plot(times, waveform.data,linestyle="-", label=wavename, color=colorcode)
            axes.legend()
            self.downx = None

            self.canvas.draw()

            self.starttime = waveform.meta['starttime']
            self.endtime = waveform.meta['endtime']


            return axes, handle


    def timeshift(self, shift):
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


    def timescale(self, delta):
        shift = (self.endtime - self.starttime) * -delta*0.1

        starttime = self.starttime + shift
        endtime = self.endtime - shift

        print(starttime, endtime,'='*8)

        for wave in self.visibleWave:
            if starttime<self.visibleWave[wave][3]:
                starttime = self.starttime
            if endtime>self.visibleWave[wave][4]:
                endtime = self.endtime

        print(starttime, endtime,'!'*8)

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

    def removewave(self, axes, handle):
        handle.pop(0).remove()
        axes.relim()
        axes.autoscale_view(True, True, True)
        axes.legend()
        self.canvas.draw()


    def mousePressEvent(self, event):
        super(TSScene, self).mousePressEvent(event)
        self.downx = event.scenePos().x()


    def mouseMoveEvent(self, event):
        if self.downx is not None:
            self.upx = event.scenePos().x()
            shift = float(self.downx - self.upx) / self.graphwidth
            self.timeshift(shift)
            self.downx=self.upx

    def mouseReleaseEvent(self, event):
        super(TSScene, self).mousePressEvent(event)
        self.downx = None

    def wheelEvent(self, event):
        super(TSScene, self).wheelEvent(event)

        delta = -event.delta() / 8 / 15

        if self.wheelactive==False:
            self.wheelactive = True
            self.timescale(delta)





    def setdata(self, filename):
        self.data = TSData(filename)

    def getlist(self):
        return self.data.getlist()


    def exportwaveform(self, wavename, filename):
        print(wavename)
        print(list(self.visibleWave))
        if wavename in self.visibleWave:
            wave = self.visibleWave[wavename][0]
            stream = self.data.getwaveform(wave, self.starttime, self.endtime)
            print(type(stream),'type of stream')
            stream.write(filename+".mseed", format='MSEED', encoding=3, reclen=256)
        else:
            pass

