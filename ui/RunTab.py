# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ressources/InputTab.ui'
#
# Created: Mon Nov  4 15:28:17 2013
#      by: PyQt4 UI code generator 4.10.3
#
# WARNING! All changes made in this file will be lost!

import os

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.Qt import QSizePolicy
from PyQt5.QtWidgets import QGridLayout, QVBoxLayout


class RunTab(QtWidgets.QWidget):
    
    def __init__(self,control):
        super(RunTab,self).__init__()
        
        self.control = control
        
        self.createMenu()
        self.createComponents()
        self.createLayout()
        self.createConnects()
    
        
    def createMenu(self):
        pass
    
    def createComponents(self):

       
        self.commandBox = QtWidgets.QTextEdit()
        self.stopButton = QtWidgets.QPushButton("Cancel Analysis!!!")
        self.processBar = QtWidgets.QProgressBar()

    def createLayout(self):
        self.centralLayout = QGridLayout()
        
        self.centralLayout.addWidget(self.commandBox,1,1,1,3)
        #self.centralLayout.addWidget(self.processBar,2)
        self.centralLayout.addWidget(self.stopButton,2,3)
        
        self.commandBox.setReadOnly(True)
        self.stopButton.setStyleSheet("background-color: red;")
        #self.stopButton.setMaximumSize(50, 25)
        
        
       
        self.setLayout(self.centralLayout)
    
    def createConnects(self):
        self.stopButton.clicked.connect(self.control.stopAssay)
          
if __name__ == '__main__':
    import sys, os
    print(os.getcwd())
    app = QtGui.QApplication(sys.argv) 
    mainWindow = RunTab()
    mainWindow.show() 
    sys.exit(app.exec_())
   
