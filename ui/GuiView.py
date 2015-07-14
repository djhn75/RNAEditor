# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ressources/guiView.ui'
#
# Created: Sat Mar 15 22:03:33 2014
#      by: PyQt4 UI code generator 4.10.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
from ui.InputTab import InputTab
from ui.RunTab import RunTab
from PyQt4.QtGui import QSizePolicy
from PyQt4.Qt import QMenu, QString


class GuiView(QtGui.QMainWindow):
    def __init__(self,control):
        self.control = control
        super(GuiView, self).__init__()
        self.createMenu()
        self.createComponents()
        self.createLayout()
        self.createConnects()
        
        
        self.setGeometry(0, 0, 400, 400)
        
        self.tabMainWindow.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(self)

    def createMenu(self):

        
        exitAction = QtGui.QAction(self)        
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)
        exitAction.setText("Exit")
        
        self.menubar = self.menuBar()
        fileMenu = self.menubar.addMenu('File')
        fileMenu.addAction(exitAction)
        """"
        fileMenu.addAction('Open File')


        self.statusBar()
        """

        
    def createComponents(self):
        self.centralWidget = QtGui.QWidget()

        self.gridLayout = QtGui.QGridLayout()
        
        self.tabMainWindow = QtGui.QTabWidget()
        self.tabMainWindow.setEnabled(True)
        self.tabMainWindow.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        #self.tabMainWindow.setMinimumSize(900, 600)
        
        self.inputTab = InputTab(self.control)
        self.tabMainWindow.addTab(self.inputTab,self.tr("InputTab"))
        
        

    def createConnects(self):
        self.tabMainWindow.tabCloseRequested.connect(self.control.closeTab)
        closeTabAction = QtGui.QAction(self.tabMainWindow)
        closeTabAction.setShortcut('Ctrl+W')
        self.connect(closeTabAction,QtCore.SIGNAL('triggered()'),self.tabMainWindow,QtCore.SLOT('close()'))
        
    def createLayout(self):
        #self.resize(679, 417)
        
        
        self.setMinimumSize(QtCore.QSize(900, 600))
        self.setStyleSheet("""
            .QWidget{border: 1px solid black}
            .DropListWidget{border: 1px solid black; background-color: white; background-image: url(ui/icons/RNAeditor_small.png); background-repeat: no-repeat; background-position:center; }
            """)
        
        

        self.tabMainWindow.setTabPosition(QtGui.QTabWidget.North)        
        self.tabMainWindow.setTabsClosable(True)
        self.gridLayout.addWidget(self.tabMainWindow, 2, 2)
        self.gridLayout.setColumnStretch(2,1)
        self.centralWidget.setLayout(self.gridLayout)
        
        self.setCentralWidget(self.centralWidget)
    
    
        