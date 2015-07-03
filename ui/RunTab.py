# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ressources/InputTab.ui'
#
# Created: Mon Nov  4 15:28:17 2013
#      by: PyQt4 UI code generator 4.10.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
from PyQt4.Qt import QSizePolicy
from PyQt4.QtGui import QGridLayout, QVBoxLayout
import os



class RunTab(QtGui.QWidget):
    
    def __init__(self,control):
        super(RunTab,self).__init__()
        
        self.connrol = control
        
        self.createMenu()
        self.createComponents()
        self.createLayout()
        self.createConnects()
    
        
    def createMenu(self):
        pass
    
    def createComponents(self):

        self.centralLayout = QtGui.QVBoxLayout()
        self.commandBox = QtGui.QTextEdit()
        self.processBar = QtGui.QProgressBar()

    def createLayout(self):
        
        self.commandBox.setReadOnly(True)
        
        self.centralLayout.addWidget(self.commandBox,1)
        self.centralLayout.addWidget(self.processBar,2)
        self.setLayout(self.centralLayout)
    
    def createConnects(self):
        pass
          
if __name__ == '__main__':
    import sys, os
    print(os.getcwd())
    app = QtGui.QApplication(sys.argv) 
    mainWindow = RunTab()
    mainWindow.show() 
    sys.exit(app.exec_())
   
