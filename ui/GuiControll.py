'''
Created on Oct 21, 2013

@author: david
'''

from ui.InputTab import InputTab
from ui.RunTab import RunTab
from PyQt4 import QtCore, QtGui
import os
from Helper import Parameters, Helper
from RnaEdit import RnaEdit
import subprocess

class GuiControll(object):
    '''
    classdocs
    '''
    
    def __init__(self, v):
        '''
        Constructor
        '''
        self.view=v
        

    @QtCore.pyqtSlot()
    def newAssay(self,):
        '''
        Function wich starts a new analysis
        '''
        
        currentInputTab = self.view.tabMainWindow.widget(0)
        fastq = currentInputTab.dropList.dropFirstItem()
        
        
        
        #check if droplist returned a value
        if fastq == None:
            QtGui.QMessageBox.information(self.view,"Warning","Warning:\nNo Sequencing Files found!!!\n\nDrop FASTQ-Files to the drop area!")
            return
        print fastq.text()
        
        #update Parameters 
        self.getParametersFromInputTab()
        self.runTab = RunTab(self)
        
        currentIndex =  self.view.tabMainWindow.count()
        self.view.tabMainWindow.addTab(self.runTab, "Analysis"+ str(currentIndex))
        
        currentTab = self.view.tabMainWindow.widget(currentIndex)
        currentTab.commandBox.append("This is your " + str(currentIndex) + " Analysis")
        
        
        fastqFiles = [str(fastq.text())]
                
        Helper.runningAssaysTabs[currentIndex]=currentTab
        """currentAssay = RnaEdit(fastqFiles, Parameters.refGenome, Parameters.dbSNP, 
                               Parameters.hapmap, Parameters.omni, Parameters.esp, 
                               Parameters.aluRegions, Parameters.gtfFile, 
                               Parameters.output, Parameters.binary, Parameters.threads, 
                               Parameters.maxDiff, Parameters.seedDiff, Parameters.paired, 
                               Parameters.standCall, Parameters.standEmit, Parameters.edgeDistance, 
                               Parameters.keepTemp, Parameters.overwrite,currentIndex)"""

        parameterList=['-i'," ".join(fastqFiles), '-r', Parameters.refGenome, '-s',Parameters.dbSNP, 
                               '-m',Parameters.hapmap, '-g',Parameters.omni, '-e',Parameters.esp, 
                               '-a',Parameters.aluRegions, '-G',Parameters.gtfFile, 
                               '-o',Parameters.output, '-d',Parameters.binary, '-t',str(Parameters.threads), 
                               '-n',str(Parameters.maxDiff), '--seedDiff',str(Parameters.seedDiff), '-sc',str(Parameters.standCall), 
                               '-se',str(Parameters.standEmit), '-ed',str(Parameters.edgeDistance), 
                               '--index',str(0)]
        
        
        
        if Parameters.paired:
            parameterList= parameterList + ["-p"]
        if Parameters.overwrite:
            parameterList = parameterList + ["--overwrite"]
        if Parameters.keepTemp:
            parameterList = parameterList + ["--keepTemp"]
        
        cmd =['python', 'RnaEdit.py'] + parameterList
       
        currentAssay = subprocess.Popen(cmd)
        
        Helper.runningAssays[currentIndex]=currentAssay
        
    def getParametersFromInputTab(self):
        '''
        get the Parameters and update the default Parameters from the Default class 
        '''
        currentInputTab = self.view.tabMainWindow.widget(0)
        
        Parameters.refGenome = str(currentInputTab.refGenomeTextBox.text())
        Parameters.dbSNP = str(currentInputTab.dbsnpTextBox.text())
        Parameters.hapmap = str(currentInputTab.hapmapTextBox.text())
        Parameters.omni = str(currentInputTab.omniTextBox.text())
        Parameters.esp = str(currentInputTab.espTextBox.text())
        Parameters.aluRegions = str(currentInputTab.aluRegionsTextBox.text())
        
        Parameters.maxDiff = str(currentInputTab.maxDiffSpinBox.value())
        Parameters.seedDiff = str(currentInputTab.seedSpinBox.value())
        Parameters.standCall = str(currentInputTab.standCallSpinBox.value())
        Parameters.standEmit = str(currentInputTab.standEmitSpinBox.value())
        Parameters.paired = currentInputTab.pairedCheckBox.isChecked()
        Parameters.overwrite = currentInputTab.overwriteCheckBox.isChecked()
        Parameters.keepTemp = currentInputTab.keepTempCheckBox.isChecked()
        
    @QtCore.pyqtSlot()
    def openFileDialog(self,textBox):
        #TODO: make single function for all necessary files
        print "open File Dialog"
        fileName = QtGui.QFileDialog.getOpenFileName(self.view.centralWidget,'Open file', QtCore.QDir.homePath())
        
        textBox.setText(fileName)
    
    @QtCore.pyqtSlot()    
    def openFolderDialog(self,textBox):
        folderName = str(QtGui.QFileDialog.getExistingDirectory(self.view.centralWidget, "Select Directory"))
        textBox.setText(folderName)
        
    @QtCore.pyqtSlot()    
    def fileDropped(self, l):
        for url in l:
            if os.path.exists(url):
                print(url)                
                icon = QtGui.QIcon(url)
                pixmap = icon.pixmap(72, 72)                
                icon = QtGui.QIcon(pixmap)
                item = QtGui.QListWidgetItem(url, self.view.inputTab.dropList)
                item.setIcon(icon)        
                item.setStatusTip(url)     
                
    @QtCore.pyqtSlot()            
    def closeTab (self, currentIndex):
        if currentIndex != 0:
            
            currentQWidget = self.view.tabMainWindow.widget(currentIndex)
            currentQWidget.deleteLater()
            self.view.tabMainWindow.removeTab(currentIndex) 