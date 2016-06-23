'''
Created on Oct 21, 2013

@author: david
'''

from ui.InputTab import InputTab
from ui.RunTab import RunTab
from PyQt4 import QtCore, QtGui
import os, sys, time
from Helper import Parameters, Helper
from RNAEditor import RnaEdit
import subprocess
import traceback

import gc
from ui.ResultTab import ResultTab
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
    def newAssay(self):
        '''
        Function wich starts a new analysis
        '''
        
        inputTab = self.view.tabMainWindow.widget(0)
        
        #get Parameters 
        parameters=Parameters(inputTab)
        if parameters.paired==True:
            #fastqs=inputTab.dropList.dropFirstTwoItems()
            fastqs = inputTab.dropList.dropFirstItem()
            if fastqs[0]!=None:
                if not str(fastqs[0].text()).endswith(".bam"):
                    fastqs+=inputTab.dropList.dropFirstItem()
        else:
            fastqs = inputTab.dropList.dropFirstItem()        
        
        """
        check if droplist returned a value
        """
        if parameters.paired==True:
            if fastqs[-1] == None:
                QtGui.QMessageBox.information(self.view,"Warning","Warning:\nNot enough Sequencing Files for paired-end sequencing!!!\n\nDrop FASTQ-Files to the drop area!")
                return
        if fastqs[0] == None:
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

        
        
        #initialize new Thread with new assay
        try:
            assay = RnaEdit(fastqFiles, parameters,runTab.commandBox)
        except Exception as err:
            QtGui.QMessageBox.information(self.view,"Error", str(err)+"Cannot start Analysis!")
            Helper.error(str(err) + "\n creating rnaEditor Object Failed!", textField=runTab.commandBox)
        currentIndex = self.view.tabMainWindow.count()

        # self.view.tabMainWindow.addTab(self.runTab, "Analysis"+ str(Helper.assayCount))
        self.view.tabMainWindow.addTab(runTab, sampleName + " " + str(currentIndex))
        Helper.runningThreads.append(assay)
        
        assay.start()
        
        self.view.connect(assay, QtCore.SIGNAL("taskDone"), self.openAnalysis)
        
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
                if url.endswith(".txt"):
                    self.view.inputTab.createDefaults(url)
                else:
                    print(url)                
                    icon = QtGui.QIcon(url)
                    pixmap = icon.pixmap(72, 72)                
                    icon = QtGui.QIcon(pixmap)
                    item = QtGui.QListWidgetItem(url, self.view.inputTab.dropList)
                    item.setIcon(icon)        
                    item.setStatusTip(url)     
            
    def openAnalysis(self,fileName=None):
        if fileName==None:
            fileName = str(QtGui.QFileDialog.getOpenFileName(self.view.centralWidget,'Open Result HTML file', QtCore.QDir.homePath(), filter = QtCore.QString("*html")))
        
        resultTab = ResultTab(self,fileName)
        self.view.tabMainWindow.addTab(resultTab,fileName[fileName.rfind("/")+1:fileName.rfind(".html")])
        Helper.runningThreads.append(resultTab)
        print(fileName)
            
    @QtCore.pyqtSlot()            
    def closeTab(self, currentIndex):
        print ("tab close %i" % currentIndex)
        if currentIndex != 0:
            currentThread=Helper.runningThreads[currentIndex]
            currentQWidget=self.view.tabMainWindow.widget(currentIndex)
            
            #check if Tab is a result Tab
            if isinstance(currentQWidget, ResultTab):
                self.view.tabMainWindow.removeTab(currentIndex)
                currentQWidget.deleteLater()
                del Helper.runningThreads[currentIndex]
                return

            #check if rnaEditor is still running or if it is finished it can be deleteted immediately
            if currentThread.isTerminated == False:
                    quitMessage = "Are you sure you want to cancel the running Sample %s?" % str(self.view.tabMainWindow.tabText(currentIndex))
                    reply = QtGui.QMessageBox.question(self.view, 'Message', quitMessage, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
                    if reply== QtGui.QMessageBox.Yes:
                        self.view.tabMainWindow.removeTab(currentIndex)
                        currentQWidget.deleteLater()
                        del Helper.runningThreads[currentIndex]
                        currentThread.stopImmediately()
            else:
                self.view.tabMainWindow.removeTab(currentIndex)
                currentQWidget.deleteLater()
                del Helper.runningThreads[currentIndex]

    def stopAssay(self,index=False):
        if index == False:
            index = self.view.tabMainWindow.currentIndex()
        if index != 0:
            currentThreat=Helper.runningThreads[index]
            #check if Thread was already terminated
            if currentThreat.isTerminated==True:
                return True

            quitMessage = "Are you sure you want to cancel the running Sample %s?" % str(self.view.tabMainWindow.tabText(index))
            reply = QtGui.QMessageBox.question(self.view, 'Message', quitMessage, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.Yes:
                currentThreat.stopImmediately()