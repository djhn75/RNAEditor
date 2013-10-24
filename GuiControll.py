'''
Created on Oct 21, 2013

@author: david
'''
import guiView
from ui.InputTab import Ui_Form as InputTab
from PyQt4 import QtCore, QtGui

class GuiControll(object):
    '''
    classdocs
    '''
    

    def __init__(self, v):
        '''
        Constructor
        '''
        self.view=v
        
    def button_pressed(self):
        print "button1"
        
    def newAssay(self):
        print "newAssay"
        self.tab = InputTab()
        self.view.tabMainWindow.addTab(self.tab, "TEST")