# Import PySide classes
import sys
from PySide.QtCore import *
from PySide.QtGui import *

# Create a Qt application
app = QApplication(sys.argv)

# Create a Window
mywindow = QWidget()
mywindow.resize(320, 240)
mywindow.setWindowTitle('Hello World!')

# Create a label and display it all together
mylabel = QLabel(mywindow)
mylabel.setText('Hello World!')
mylabel.setGeometry(QRect(130, 110, 60, 10))
mywindow.show()

# Enter Qt application main loop
sys.exit(app.exec_())
