# -*- coding: utf-8 -*-

import sys 
from PyQt4 import QtGui
from guiView import Ui_MainWindow as View
from GuiControll import GuiControll

class RnaEditor(QtGui.QMainWindow, View): 
    def __init__(self):
        self.control = GuiControll(self) #create controller class
        QtGui.QDialog.__init__(self) 
        self.setupUi(self) #create user interface
        
app = QtGui.QApplication(sys.argv) 
dialog = RnaEditor() 
dialog.show() 
sys.exit(app.exec_())