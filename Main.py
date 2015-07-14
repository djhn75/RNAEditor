# -*- coding: utf-8 -*-

import sys 
from PyQt4 import QtGui, QtCore
from ui.GuiView import GuiView
from ui.GuiControll import GuiControll

def main(argv):
    app = QtGui.QApplication(argv) 
    mainWindow = RnaEditor()
    
    app.setApplicationName("RNAEditor")
    app.setApplicationVersion("0.1")
    
    app_icon = QtGui.QIcon()
    app_icon.addFile('ui/icons/rnaEditor_16x16.png', QtCore.QSize(16,16))
    app_icon.addFile('ui/icons/rnaEditor_24x24.png', QtCore.QSize(24,24))
    app_icon.addFile('ui/icons/rnaEditor_32x32.png', QtCore.QSize(32,32))
    app_icon.addFile('ui/icons/rnaEditor_48x48.png', QtCore.QSize(48,48))
    app_icon.addFile('ui/icons/rnaEditor_256x256.png', QtCore.QSize(256,256))
    app_icon.addFile('ui/icons/rnaEditor_512x512.png', QtCore.QSize(512,512))
    app_icon.addFile('ui/icons/rnaEditor_1024x1024.png', QtCore.QSize(1024,1024))
    
    app.setWindowIcon(app_icon)

     
    mainWindow.show() 
    sys.exit(app.exec_())

class RnaEditor(GuiView): 
    assayList=[]
    
    def __init__(self):
        self.control = GuiControll(self) #create controller class
        super(RnaEditor, self).__init__(self.control)
        
        #set window title
        self.setWindowTitle("RnaEditor")
        

        #set default Values
        self.inputTab.createDefaults()
if __name__ == "__main__":
    main(sys.argv)
