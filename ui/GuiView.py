from PyQt5 import QtCore, QtGui, QtWidgets
from ui.InputTab import InputTab
from PyQt5.QtWidgets import QSizePolicy
from ui.GuiControll import GuiControll


class GuiView(QtWidgets.QMainWindow):
    def __init__(self):
        self.control = GuiControll(self)
        super(GuiView, self).__init__()
        
        self.createComponents()
        self.createLayout()
        self.createMenu()        
        self.createConnects()
        
        
        self.setWindowTitle("RnaEditor")
        self.setGeometry(0, 0, 400, 400)
        
        self.tabMainWindow.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(self)

    def createMenu(self):
        self.exitAction.setShortcut('Ctrl+Q')
        self.exitAction.setStatusTip('Exit application')
        self.exitAction.triggered.connect(QtWidgets.qApp.quit)
        self.exitAction.setText("Exit")
        
        self.menubar = self.menuBar()
        self.fileMenu = self.menubar.addMenu('File') 
        self.fileMenu.addAction(self.exitAction)
        
        self.openAnalysisAction = QtWidgets.QAction("Open Analysis",self)
        self.fileMenu.addAction(self.openAnalysisAction)
        

        self.statusBar()
 
    def createComponents(self):
        self.centralWidget = QtWidgets.QWidget()
        self.gridLayout = QtWidgets.QGridLayout()
        self.tabMainWindow = QtWidgets.QTabWidget()
        self.exitAction = QtWidgets.QAction(self)        
        
        
        self.tabMainWindow.setEnabled(True)
        self.tabMainWindow.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
        #self.tabMainWindow.setMinimumSize(900, 600)
        
        self.inputTab = InputTab(self.control)
        self.tabMainWindow.addTab(self.inputTab,self.tr("InputTab"))
        
        self.inputTab.createDefaults()

    def createConnects(self):
       
        self.closeTabAction = QtWidgets.QAction("Close Tab", shortcut=QtGui.QKeySequence("Ctrl+w"), triggered=lambda:self.control.closeTab)
        self.tabMainWindow.addAction(self.closeTabAction)
        self.tabMainWindow.tabCloseRequested.connect(self.control.closeTab)
        self.closeTabAction.triggered.connect(self.tabMainWindow.close)

        #self.connect(closeTabAction,QtCore.SIGNAL('triggered()'),self.tabMainWindow,QtCore.SLOT('close()'))        
        self.openAnalysisAction.triggered.connect(self.control.openAnalysis)
        
    def createLayout(self):
        #self.resize(679, 417)
        
        
        self.setMinimumSize(QtCore.QSize(900, 600))
        self.setStyleSheet("""
            .QWidget{border: 1px solid grey}
            .DropListWidget{border: 1px solid black; background-color: #f2f2f2; background-image: url(inputTab_icon.png); background-repeat: no-repeat; background-position:center; background-size:cover;}
            """)
        
        

        self.tabMainWindow.setTabPosition(QtWidgets.QTabWidget.North)        
        self.tabMainWindow.setTabsClosable(True)
        self.gridLayout.addWidget(self.tabMainWindow, 2, 2)
        self.gridLayout.setColumnStretch(2,1)
        self.centralWidget.setLayout(self.gridLayout)
        
        self.setCentralWidget(self.centralWidget)
    
    
        