# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ressources/guiView.ui'
#
# Created: Sun Nov  3 23:02:12 2013
#      by: PyQt4 UI code generator 4.10.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui
from ui import InputTab

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(679, 417)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMinimumSize(QtCore.QSize(300, 300))
        MainWindow.setStyleSheet(_fromUtf8(""))
        self.centralWidget = QtGui.QWidget(MainWindow)
        self.centralWidget.setObjectName(_fromUtf8("centralWidget"))
        self.gridLayoutWidget = QtGui.QWidget(self.centralWidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(130, -1, 511, 341))
        self.gridLayoutWidget.setObjectName(_fromUtf8("gridLayoutWidget"))
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.tabMainWindow = QtGui.QTabWidget(self.gridLayoutWidget)
        self.tabMainWindow.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.tabMainWindow.sizePolicy().hasHeightForWidth())
        self.tabMainWindow.setSizePolicy(sizePolicy)
        self.tabMainWindow.setTabPosition(QtGui.QTabWidget.North)
        self.tabMainWindow.setElideMode(QtCore.Qt.ElideRight)
        self.tabMainWindow.setUsesScrollButtons(False)
        self.tabMainWindow.setTabsClosable(True)
        self.tabMainWindow.setMovable(True)
        self.tabMainWindow.setObjectName(_fromUtf8("tabMainWindow"))
        self.tab = InputTab.Ui_Form()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.tabMainWindow.addTab(self.tab, _fromUtf8(""))
        self.gridLayout.addWidget(self.tabMainWindow, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtGui.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 679, 22))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuRnaEditor = QtGui.QMenu(self.menuBar)
        self.menuRnaEditor.setObjectName(_fromUtf8("menuRnaEditor"))
        self.menuFile = QtGui.QMenu(self.menuBar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuRessource = QtGui.QMenu(self.menuBar)
        self.menuRessource.setObjectName(_fromUtf8("menuRessource"))
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtGui.QToolBar(MainWindow)
        self.mainToolBar.setObjectName(_fromUtf8("mainToolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.toolBar = QtGui.QToolBar(MainWindow)
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionAbout_RnaEditor = QtGui.QAction(MainWindow)
        self.actionAbout_RnaEditor.setObjectName(_fromUtf8("actionAbout_RnaEditor"))
        self.actionQuit_RnaEditor_2 = QtGui.QAction(MainWindow)
        self.actionQuit_RnaEditor_2.setObjectName(_fromUtf8("actionQuit_RnaEditor_2"))
        self.actionNew_Tab = QtGui.QAction(MainWindow)
        self.actionNew_Tab.setObjectName(_fromUtf8("actionNew_Tab"))
        self.actionOpen_File = QtGui.QAction(MainWindow)
        self.actionOpen_File.setObjectName(_fromUtf8("actionOpen_File"))
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.actionDownload_Ressources = QtGui.QAction(MainWindow)
        self.actionDownload_Ressources.setObjectName(_fromUtf8("actionDownload_Ressources"))
        self.actionSet_Genome = QtGui.QAction(MainWindow)
        self.actionSet_Genome.setObjectName(_fromUtf8("actionSet_Genome"))
        self.actionSet_SNP_Ressource = QtGui.QAction(MainWindow)
        self.actionSet_SNP_Ressource.setObjectName(_fromUtf8("actionSet_SNP_Ressource"))
        self.actionSet_Alu_Regions = QtGui.QAction(MainWindow)
        self.actionSet_Alu_Regions.setObjectName(_fromUtf8("actionSet_Alu_Regions"))
        self.actionSet_Simple_Repeats = QtGui.QAction(MainWindow)
        self.actionSet_Simple_Repeats.setObjectName(_fromUtf8("actionSet_Simple_Repeats"))
        self.actionSet_Ressource_Folder = QtGui.QAction(MainWindow)
        self.actionSet_Ressource_Folder.setObjectName(_fromUtf8("actionSet_Ressource_Folder"))
        self.actionNew_Assay = QtGui.QAction(MainWindow)
        self.actionNew_Assay.setObjectName(_fromUtf8("actionNew_Assay"))
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

        self.retranslateUi(MainWindow)
        self.tabMainWindow.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.tabMainWindow.setTabText(self.tabMainWindow.indexOf(self.tab), _translate("MainWindow", "Tab 1", None))
        self.menuRnaEditor.setTitle(_translate("MainWindow", "RnaEditor", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuRessource.setTitle(_translate("MainWindow", "Ressource", None))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar", None))
        self.actionAbout_RnaEditor.setText(_translate("MainWindow", "About RnaEditor", None))
        self.actionQuit_RnaEditor_2.setText(_translate("MainWindow", "Quit RnaEditor", None))
        self.actionNew_Tab.setText(_translate("MainWindow", "New Tab", None))
        self.actionOpen_File.setText(_translate("MainWindow", "open File", None))
        self.actionOpen.setText(_translate("MainWindow", "Open FastQ... ", None))
        self.actionDownload_Ressources.setText(_translate("MainWindow", "Download Ressources...", None))
        self.actionSet_Genome.setText(_translate("MainWindow", "Set Genome...", None))
        self.actionSet_SNP_Ressource.setText(_translate("MainWindow", "Set SNP Ressource...", None))
        self.actionSet_Alu_Regions.setText(_translate("MainWindow", "Set Alu-Regions...", None))
        self.actionSet_Simple_Repeats.setText(_translate("MainWindow", "Set Simple Repeats...", None))
        self.actionSet_Ressource_Folder.setText(_translate("MainWindow", "Set Ressource Folder...", None))
        self.actionNew_Assay.setText(_translate("MainWindow", "New Assay", None))

