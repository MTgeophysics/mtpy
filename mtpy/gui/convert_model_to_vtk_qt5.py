# -*- coding: utf-8 -*-
"""
Created on Fri Nov 04 14:56:12 2016

@author: jpeacock

Gui to transform model data into vtk
"""

#==============================================================================
# Imports
#==============================================================================
# standard packages
import os
import sys
import subprocess

# 3rd party packages
from PyQt5 import QtCore, QtGui, QtWidgets

import mtpy.modeling.modem
from mtpy.gui.my_stream import MyStream
import mtpy.modeling.modem as modem
#==============================================================================
class ConvertModel2VTK(QtWidgets.QWidget):
    """
    
    """
    
    def __init__(self):
        super(ConvertModel2VTK, self).__init__()
        
        self.model_fn = None
        self.resp_fn = None
        self.cwd = os.getcwd()
        
        self._model_list = ['ModEM', 'WS']
        self.model_type = self._model_list[0]
        self._mfn_list = []
        self._rfn_list = []
        
        self.vtk_model_fn = None
        self.vtk_station_fn = None
        
        
        self.setup_ui()
        
    def setup_ui(self):
        self.setWindowTitle('Convert Model To VTK')
        

        self.cwd_button = QtWidgets.QPushButton("Working Directory")
        self.cwd_button.pressed.connect(self.get_cwd)
        self.cwd_edit = QtWidgets.QLineEdit()
        self.cwd_edit.editingFinished.connect(self.set_cwd)
        
        self.model_fn_label = QtWidgets.QLabel("Model File")
        self.model_fn_combo = QtWidgets.QComboBox()
        self.model_fn_combo.setMinimumWidth(500)
        self.model_fn_combo.currentIndexChanged.connect(self.set_model_fn)
        
        self.resp_fn_label = QtWidgets.QLabel("Response File")
        self.resp_fn_combo = QtWidgets.QComboBox()
        self.resp_fn_combo.setMinimumWidth(500)
        self.resp_fn_combo.currentIndexChanged.connect(self.set_resp_fn)
        
        self.vtk_model_fn_label = QtWidgets.QLabel("VTK Model File")
        self.vtk_model_fn_edit = QtWidgets.QLineEdit()
        self.vtk_model_fn_edit.editingFinished.connect(self.set_vtk_model_fn)
        
        self.vtk_station_fn_label = QtWidgets.QLabel("VTK Station File")
        self.vtk_station_fn_edit = QtWidgets.QLineEdit()
        self.vtk_station_fn_edit.editingFinished.connect(self.set_vtk_station_fn)
        
        self.make_vtk_button = QtWidgets.QPushButton("Make VTK Files")
        self.make_vtk_button.setStyleSheet("background-color: #ed939f")
        self.make_vtk_button.pressed.connect(self.make_vtk_files)
        
        self.model_type_combo = QtWidgets.QComboBox()
        self.model_type_combo.addItems(self._model_list)
        self.model_type_combo.currentIndexChanged.connect(self.set_model_type)
        
        ## add in model center to put into geographic coordinates
        
        self.output_box = QtWidgets.QTextEdit()        
        
        self.my_stream = MyStream()
        self.my_stream.message.connect(self.normal_output)
        
        sys.stdout = self.my_stream
        

        # --> layout        
        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.cwd_button, 0, 0)
        layout.addWidget(self.cwd_edit, 0, 1)
        layout.addWidget(self.model_fn_label, 1, 0)
        layout.addWidget(self.model_fn_combo, 1, 1)
        layout.addWidget(self.resp_fn_label, 2, 0)
        layout.addWidget(self.resp_fn_combo, 2, 1)
        layout.addWidget(self.vtk_model_fn_label, 3, 0)
        layout.addWidget(self.vtk_model_fn_edit, 3, 1)
        layout.addWidget(self.vtk_station_fn_label, 4, 0)
        layout.addWidget(self.vtk_station_fn_edit, 4, 1)
        layout.addWidget(self.model_type_combo, 5, 0)
        layout.addWidget(self.make_vtk_button, 5, 1)
        
        final_layout = QtWidgets.QVBoxLayout()
        final_layout.addLayout(layout)
        final_layout.addWidget(self.output_box)
        
        self.setLayout(final_layout)
        
        QtCore.QMetaObject.connectSlotsByName(self)

    def set_cwd(self):
        self.cwd = os.path.abspath(str(self.cwd_edit.text()))
        os.chdir(self.cwd)
        if not os.path.exists(self.cwd):
            print(('Path doesnt exist {0}'.format(self.cwd)))
            
        self.set_fn_lists()
    
    def get_cwd(self):
        dir_dialog = QtWidgets.QFileDialog()
        self.cwd = os.path.abspath(str(dir_dialog.getExistingDirectory(caption='Choose Model Directory')))
        self.cwd_edit.setText(self.cwd)
        os.chdir(self.cwd)
        self.set_fn_lists()
                            
    def set_fn_lists(self):
        
        if self.model_type == 'ModEM':
            self._mfn_list = [fn for fn in os.listdir(self.cwd) if fn.endswith('.rho')]
            self._rfn_list = [fn for fn in os.listdir(self.cwd) if fn.endswith('.dat')]
        elif self.model_type == 'WS':
            self._mfn_list = [fn for fn in os.listdir(self.cwd) if fn.find('model')>=0]
            self._rfn_list = [fn for fn in os.listdir(self.cwd) if fn.find('resp')>=0]
            
        self.model_fn_combo.clear()
        self.model_fn_combo.addItems(self._mfn_list)
        self.resp_fn_combo.clear()
        self.resp_fn_combo.addItems(self._rfn_list)
        
    def set_model_fn(self, index):
        self.model_fn = self._mfn_list[index]
        
        if self.model_type == 'ModEM':
            self.vtk_model_fn = '{0}_res'.format(self.model_fn[:-4])
        elif self.model_type == 'WS':
            self.vtk_model_fn = '{0}_res'.format(self.model_fn[:self.model_fn.find('.')])

        self.vtk_model_fn_edit.setText(self.vtk_model_fn)
 
    def set_resp_fn(self, index):
        self.resp_fn = self._rfn_list[index] 
        
        if self.model_type == 'ModEM':
            self.vtk_station_fn = '{0}_stations'.format(self.resp_fn[:-4])
        elif self.model_type == 'WS':
            self.vtk_station_fn = '{0}_stations'.format(self.resp_fn[:self.model_fn.find('.')])
          
        self.vtk_station_fn_edit.setText(self.vtk_station_fn)
        
    def set_vtk_model_fn(self):
        self.vtk_model_fn = str(self.vtk_model_fn_edit.text())
        
    def set_vtk_station_fn(self):
        self.vtk_station_fn = str(self.vtk_station_fn_edit.text())
        
    def set_model_type(self, index):
        self.model_type = self._model_list[index]
        
    def make_vtk_files(self):
        if self.model_type == 'ModEM':
            m_obj = mtpy.modeling.modem.Model()
            m_obj.read_model_file(os.path.join(self.cwd, self.model_fn))
            m_obj.write_vtk_file(vtk_save_path=self.cwd,
                                 vtk_fn_basename=self.vtk_model_fn)
                                 
            d_obj = mtpy.modeling.modem.Data()
            d_obj.read_data_file(self.resp_fn)
            d_obj.write_vtk_station_file(vtk_save_path=self.cwd,
                                         vtk_fn_basename=self.vtk_station_fn)
            program = 'modem2vtk'
        elif self.model_type == 'WS':
            program = 'ws2vtk'
            
            subprocess.call([program, 
                             self.model_fn, 
                             self.resp_fn,
                             self.vtk_model_fn,
                             self.vtk_station_fn])
        
    @QtCore.pyqtSlot(str)
    def normal_output(self, message):
        self.output_box.moveCursor(QtGui.QTextCursor.End)
        self.output_box.insertPlainText(message)

#==============================================================================
# create main            
#==============================================================================
def main():
    app = QtWidgets.QApplication(sys.argv)
    ui = ConvertModel2VTK()
    ui.show()
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    main()
