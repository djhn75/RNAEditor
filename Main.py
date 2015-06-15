# -*- coding: utf-8 -*-

import sys 
from PyQt4 import QtGui
from ui.GuiView import GuiView
from ui.GuiControll import GuiControll

def main(argv):
    app = QtGui.QApplication(argv) 
    mainWindow = RnaEditor() 
    mainWindow.show() 
    sys.exit(app.exec_())

class RnaEditor(GuiView): 
    assayList=[]
    
    def __init__(self):
        self.control = GuiControll(self) #create controller class
        super(RnaEditor, self).__init__(self.control)


        #set default Values
        self.inputTab.createDefaults()
if __name__ == "__main__":
    main(sys.argv)
