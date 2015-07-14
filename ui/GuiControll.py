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
    def __init__(self,fastqFiles,parameters,textField):
        QtCore.QThread.__init__(self)

        try:
            self.assay = RnaEdit(fastqFiles, parameters,textField)
        except Exception:
            
            Helper.error("creating rnaEditor Object Failed!" ,textField=textField)
            #print "rnaEditor Failed!"

    def run(self):
        #print "Start Thread" + str(self.fastqFiles)

        try:
            self.assay.start()
        except Exception:
            Helper.error("RnaEditor Failed")


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
        inputTab = self.view.tabMainWindow.widget(0)
        
        #get Parameters 
        parameters=Parameters(inputTab)
        if parameters.paired==True:
            #TODO: add support for paired-end
            fastqs=inputTab.dropList.dropFirstTwoItems()
        else:
            fastqs = inputTab.dropList.dropFirstItem()        
        
        """
        check if droplist returned a value
        """
        if fastqs == None:
            QtGui.QMessageBox.information(self.view,"Warning","Warning:\nNo Sequencing Files found!!!\n\nDrop FASTQ-Files to the drop area!")
            return
        sampleName = Helper.getSampleName(str(fastqs[0].text()))
        if sampleName == None:
            QtGui.QMessageBox.information(self.view,"Warning","Warning:\nNo valid Sequencing File!!!\n\nDrop FASTQ-Files to the drop area!")
            return
        
        fastqFiles=[]
        for fastq in fastqs:
            fastqFiles.append(str(fastq.text()))

        
        
        runTab = RunTab(self)
        currentIndex =  self.view.tabMainWindow.count()

        #self.view.tabMainWindow.addTab(self.runTab, "Analysis"+ str(Helper.assayCount))
        self.view.tabMainWindow.addTab(runTab, sampleName + " "+ str(currentIndex))
        
        
        #initialize new Thread with new assay
        try:
            workThread=WorkThread(fastqFiles,parameters,runTab.commandBox)
        except:
            return
        
        Helper.runningThreads.append(workThread)
        
        workThread.start()
        
        
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
            currentThreat=Helper.runningThreads[currentIndex]
            currentQWidget=self.view.tabMainWindow.widget(currentIndex)
            #check if rnaEditor is still running
            if currentThreat.assay.runningCommand != False:
                
                
                quitMessage = "Are you sure you want to quit the running Sample %s?" % str(self.view.tabMainWindow.tabText(currentIndex))
                reply = QtGui.QMessageBox.question(self.view, 'Message', quitMessage, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
    
                if reply == QtGui.QMessageBox.Yes:
                    self.deleteAssay(currentThreat.assay)
                    currentThreat.wait()
                    self.view.tabMainWindow.removeTab(currentIndex)
                    currentQWidget.deleteLater()
                    del Helper.runningThreads[currentIndex]
                    
            else:
                self.deleteAssay(currentThreat.assay)
                currentThreat.quit()
                #currentThreat.wait()
                self.view.tabMainWindow.removeTab(currentIndex)
                currentQWidget.deleteLater()
                del Helper.runningThreads[currentIndex]
    def deleteAssay(self,assay):
        '''
        Finally removes an Tab and deletes the assay from all global Arrays
        :param currentIndex: Index of the Tab and Assay which should be removed
        '''
        print "deleteAssay" + str(assay)
        
        if assay.runningCommand != False:
            print "kill:" + str(assay.runningCommand)
            assay.runningCommand.kill()
        
        del assay
