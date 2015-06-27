'''
Created on Oct 21, 2013

@author: david
'''

from ui.InputTab import InputTab
from ui.RunTab import RunTab
from PyQt4 import QtCore, QtGui
import os, sys, time
from Helper import Parameters, Helper
from RnaEdit import RnaEdit
import subprocess
import traceback



class WorkThread(QtCore.QThread):
    def __init__(self,fastqFiles,currentIndex,textField):
        QtCore.QThread.__init__(self)
        self.fastqFiles=fastqFiles
        self.textField=textField
        self.currentIndex=currentIndex
        try:
            self.currentAssay = RnaEdit(self.fastqFiles, Parameters.refGenome, Parameters.dbSNP, 
                                   Parameters.hapmap, Parameters.omni, Parameters.esp, 
                                   Parameters.aluRegions, Parameters.gtfFile, 
                                   Parameters.output, Parameters.binary, Parameters.threads, 
                                   Parameters.maxDiff, Parameters.seedDiff, Parameters.paired, 
                                   Parameters.standCall, Parameters.standEmit, Parameters.edgeDistance, 
                                   Parameters.keepTemp, Parameters.overwrite,self.textField)
        except Exception:
            print(traceback.format_exc())
            Helper.error("creating rnaEditor Object Failed!", self.currentAssay.logFile, self.textField)
            #print "rnaEditor Failed!"
            
        
    def __del__(self):
        self.wait()
    
    def run(self):
        #print "Start Thread" + str(self.fastqFiles)

        try:
            Helper.assays[self.currentIndex].start()
        except Exception:
            print(traceback.format_exc())
            Helper.assayIsRunning[self.currentIndex]=False
            Helper.error("rnaEditor Failed!", self.currentAssay.logFile, self.textField)
            
            
        

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
        Helper.assayCount+=1
        currentInputTab = self.view.tabMainWindow.widget(0)
        fastq = currentInputTab.dropList.dropFirstItem()
        
        if fastq == None:
            QtGui.QMessageBox.information(self.view,"Warning","Warning:\nNo Sequencing Files found!!!\n\nDrop FASTQ-Files to the drop area!")
            return
        sampleName = Helper.getSampleName(str(fastq.text()))
        if sampleName == None:
            QtGui.QMessageBox.information(self.view,"Warning","Warning:\nNo valid Sequencing File!!!\n\nDrop FASTQ-Files to the drop area!")
            return
        
        #check if droplist returned a value
        
        
        #update Parameters 
        self.getParametersFromInputTab()
        self.runTab = RunTab(self)
        
        currentIndex =  self.view.tabMainWindow.count()

        #self.view.tabMainWindow.addTab(self.runTab, "Analysis"+ str(Helper.assayCount))
        self.view.tabMainWindow.addTab(self.runTab, sampleName + " "+ str(Helper.assayCount) + " " + str(currentIndex))
        
        currentTab = self.view.tabMainWindow.widget(currentIndex)
        fastqFiles = [str(fastq.text())]
                
        Helper.runningAssaysTabs.append(currentTab)
        
        #initialize new Thread with new assay
        try:
            workThread=WorkThread(fastqFiles,Helper.assayCount,self.runTab.commandBox)
        except:
            raise
        
        Helper.runningAssaysThreads.append(workThread)
        Helper.assays.append(workThread.currentAssay)
        
        Helper.runningAssaysThreads[Helper.assayCount].start()
        Helper.assayIsRunning.append(True)
        
        Helper.runningCommand.append(False)
        
        
        
        
    def getParametersFromInputTab(self):
        '''
        get the Parameters and update the default Parameters from the Default class 
        '''
        currentInputTab = self.view.tabMainWindow.widget(0)
        
        
        Parameters.refGenome = str(currentInputTab.refGenomeTextBox.text())
        Parameters.gtfFile = str(currentInputTab.gtfFileTextBox.text())
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
        
        
        print "tab close %i" % currentIndex
        if currentIndex != 0:
            #check if rnaEditor is still running
            if len(Helper.assayIsRunning) > currentIndex and Helper.assayIsRunning[currentIndex]== True:
                
                quitMessage = "Are you sure you want to quit the running Sample %s?" % str(self.view.tabMainWindow.tabText(currentIndex))
                reply = QtGui.QMessageBox.question(self.view, 'Message', quitMessage, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
    
                if reply == QtGui.QMessageBox.Yes:
                    self.deleteAssay(currentIndex)
            else:
                self.deleteAssay(currentIndex)
                
    def deleteAssay(self,currentIndex):
        '''
        Finally removes an Tab and deletes the assay from all global Arrays
        :param currentIndex: Index of the Tab and Assay which should be removed
        '''
        print "deleteAssay" + str(currentIndex)
        Helper.assayCount-=1
        currentQWidget = self.view.tabMainWindow.widget(currentIndex) 
        if len(Helper.runningCommand)>currentIndex:
            if Helper.runningCommand[currentIndex] != False:
                print "kill:" + str(Helper.runningCommand[currentIndex])
                Helper.runningCommand[currentIndex].kill()
        
        #wait for the thread to finish
        if len(Helper.runningAssaysThreads) > currentIndex:
            Helper.runningAssaysThreads[currentIndex].wait()
            del Helper.runningAssaysThreads[currentIndex]    
        if len(Helper.assays) > currentIndex:
            del Helper.assays[currentIndex]
        if len(Helper.runningAssaysTabs) > currentIndex:
            del Helper.runningAssaysTabs[currentIndex]

        if len(Helper.runningCommand) > currentIndex:
            del Helper.runningCommand[currentIndex]
        self.view.tabMainWindow.removeTab(currentIndex)
        currentQWidget.deleteLater()
