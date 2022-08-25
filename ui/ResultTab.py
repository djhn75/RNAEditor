import sys

from PyQt5 import QtGui
from PyQt5.QtCore import QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView


class ResultTab(QWebEngineView):
    
    def __init__(self,control,site):
        super(ResultTab,self).__init__()
        
        self.control = control
        self.site = str(site)
        self.createMenu()
        self.createComponents()
        self.createLayout()
        self.createConnects()
    
        
    def createMenu(self):
        pass
    
    def createComponents(self):
        with open(self.site,'r') as f:
            html=f.read()
            self.setHtml(html)
       
     
    def _result_available(self,ok):
        pass
        
           
    def createLayout(self):
        pass
        
    
    def createConnects(self):
        pass
          
if __name__ == '__main__':
    import sys, os
    print(os.getcwd())
    app = QtGui.QApplication(sys.argv) 
    mainWindow = ResultTab(None,'http://google.com')
    mainWindow.show() 
    sys.exit(app.exec_())
   
