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


        exitAction = QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)        
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(QtGui.qApp.quit)
        
        menubar = self.menuBar()
        fileMenu = menubar.addMenu('File')
        fileMenu.addAction(exitAction)

        """        
        self.menuBar = QtGui.QMenuBar()
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 679, 222))
        self.menuRnaEditor = QtGui.QMenu(self.menuBar)
        self.menuRnaEditor.setObjectName(self.tr("menuRnaEditor"))
        self.menuFile = QtGui.QMenu(self.menuBar)
        
        self.menuRessource = QtGui.QMenu(self.menuBar)"""
        
        """        
        self.setMenuBar(self.menuBar)
        self.mainToolBar = QtGui.QToolBar(self)
        
        self.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.toolBar = QtGui.QToolBar(self)
        
        self.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionAbout_RnaEditor = QtGui.QAction(self)
        
        self.actionQuit_RnaEditor_2 = QtGui.QAction(self)
        self.actionNew_Tab = QtGui.QAction(self)
        self.actionOpen_File = QtGui.QAction(self)
        self.actionOpen = QtGui.QAction(self)
        
        self.actionDownload_Ressources = QtGui.QAction(self)
        
        self.actionSet_Genome = QtGui.QAction(self)
        self.actionSet_SNP_Ressource = QtGui.QAction(self)
        self.actionSet_Alu_Regions = QtGui.QAction(self)
        self.actionSet_Simple_Repeats = QtGui.QAction(self)
        self.actionSet_Ressource_Folder = QtGui.QAction(self)

        self.actionNew_Assay = QtGui.QAction(self)

        self.menuRnaEditor.addAction(self.actionAbout_RnaEditor)
        self.menuRnaEditor.addSeparator()
        self.menuRnaEditor.addAction(self.actionQuit_RnaEditor_2)
        self.menuFile.addAction(self.actionNew_Assay)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuRessource.addAction(self.actionDownload_Ressources)
        self.menuRessource.addSeparator()
        self.menuRessource.addAction(self.actionSet_Ressource_Folder)
        self.menuRessource.addSeparator()
        self.menuRessource.addAction(self.actionSet_Genome)
        self.menuRessource.addAction(self.actionSet_SNP_Ressource)
        self.menuRessource.addAction(self.actionSet_Alu_Regions)
        self.menuRessource.addAction(self.actionSet_Simple_Repeats)
        self.menuBar.addAction(self.menuRnaEditor.menuAction())
        self.menuBar.addAction(self.menuFile.menuAction())
        self.menuBar.addAction(self.menuRessource.menuAction())
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
    
    
        