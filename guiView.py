# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'guiView.ui'
#
# Created: Fri Oct 11 16:25:19 2013
#      by: PyQt4 UI code generator 4.10.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

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
        MainWindow.resize(715, 407)
        self.centralWidget = QtGui.QWidget(MainWindow)
        self.centralWidget.setObjectName(_fromUtf8("centralWidget"))
        self.lineEdit = QtGui.QLineEdit(self.centralWidget)
        self.lineEdit.setGeometry(QtCore.QRect(90, 280, 431, 21))
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.pushButton = QtGui.QPushButton(self.centralWidget)
        self.pushButton.setGeometry(QtCore.QRect(530, 270, 114, 32))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.tableView = QtGui.QTableView(self.centralWidget)
        self.tableView.setGeometry(QtCore.QRect(90, 40, 431, 221))
        self.tableView.setObjectName(_fromUtf8("tableView"))
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QtGui.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 715, 22))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuRnaEditor = QtGui.QMenu(self.menuBar)
        self.menuRnaEditor.setObjectName(_fromUtf8("menuRnaEditor"))
        self.menuFile = QtGui.QMenu(self.menuBar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QtGui.QToolBar(MainWindow)
        self.mainToolBar.setObjectName(_fromUtf8("mainToolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtGui.QStatusBar(MainWindow)
        self.statusBar.setObjectName(_fromUtf8("statusBar"))
        MainWindow.setStatusBar(self.statusBar)
        self.actionAbout_RnaEditor = QtGui.QAction(MainWindow)
        self.actionAbout_RnaEditor.setObjectName(_fromUtf8("actionAbout_RnaEditor"))
        self.actionQuit_RnaEditor_2 = QtGui.QAction(MainWindow)
        self.actionQuit_RnaEditor_2.setObjectName(_fromUtf8("actionQuit_RnaEditor_2"))
        self.actionNew_Tab = QtGui.QAction(MainWindow)
        self.actionNew_Tab.setObjectName(_fromUtf8("actionNew_Tab"))
        self.actionOpen_File = QtGui.QAction(MainWindow)
        self.actionOpen_File.setObjectName(_fromUtf8("actionOpen_File"))
        self.menuRnaEditor.addAction(self.actionAbout_RnaEditor)
        self.menuRnaEditor.addSeparator()
        self.menuRnaEditor.addAction(self.actionQuit_RnaEditor_2)
        self.menuFile.addAction(self.actionNew_Tab)
        self.menuFile.addSeparator()
        self.menuBar.addAction(self.menuRnaEditor.menuAction())
        self.menuBar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.pushButton, QtCore.SIGNAL(_fromUtf8("clicked()")), self.view.button_pressed)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
    
    def test(self):
        print "test"
    
    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.pushButton.setText(_translate("MainWindow", "X", None))
        self.menuRnaEditor.setTitle(_translate("MainWindow", "RnaEditor", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.actionAbout_RnaEditor.setText(_translate("MainWindow", "About RnaEditor", None))
        self.actionQuit_RnaEditor_2.setText(_translate("MainWindow", "Quit RnaEditor", None))
        self.actionNew_Tab.setText(_translate("MainWindow", "New Tab", None))
        self.actionOpen_File.setText(_translate("MainWindow", "open File", None))

