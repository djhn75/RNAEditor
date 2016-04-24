import sys

from PyQt4 import QtGui
from PyQt4.QtCore import QUrl
from PyQt4.QtWebKit import QWebView


class ResultTab(QWebView):
    
    def __init__(self,control,site):
        super(ResultTab,self).__init__()
        
        self.control = control
        self.site = site
        self.createMenu()
        self.createComponents()
        self.createLayout()
        self.createConnects()
    
        
    def createMenu(self):
        pass
    
    def createComponents(self):
        self.load(QUrl(self.site))
       
     
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
   
