# -*- coding: utf-8 -*-

import sys 
from PyQt4 import QtGui
from guiView import Ui_MainWindow as View

class MyDialog(QtGui.QMainWindow, View): 
    def __init__(self): 
        QtGui.QDialog.__init__(self) 
        self.setupUi(self)

app = QtGui.QApplication(sys.argv) 
dialog = MyDialog() 
dialog.show() 
sys.exit(app.exec_())