# -*- coding: utf-8 -*-
#!/usr/bin/python


import sys
from sys import platform
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from pylab import *
from pylab import *
import scipy
from scipy import ndimage,optimize, signal
import scipy.ndimage.interpolation as spni
import scipy.fftpack as spf
from scipy.signal import *
from skimage.feature import match_template
from PIL import Image
import os
from os.path import isfile, join
import string
import pyqtgraph as pg
from PySide import QtGui, QtCore
from pyqtgraph import QtGui, QtCore
import h5py
import tomopy
import time
#import subpixelshift
#from subpixelshift import *


class Example(QtGui.QMainWindow):

      def __init__(self):
            super(Example, self).__init__()
            self.initUI()
        
      def initUI(self):
              
 
##            textEdit = QtGui.QTextEdit()
##            self.setCentralWidget(textEdit)
            
            self.ImageTag = "exchange"
            self.thetaPos = 657
            self.thetaPos = 8


            exitAction = QtGui.QAction('Exit', self)
            exitAction.triggered.connect(self.close)
            exitAction.setShortcut('Ctrl+Q')

            closeAction = QtGui.QAction('Quit', self)
            closeAction.triggered.connect(sys.exit)
            closeAction.setShortcut('Ctrl+X')

            openFileAction =QtGui.QAction('Open File',self)
            openFileAction.triggered.connect(self.openfile)

            openFolderAction = QtGui.QAction('Open Folder', self)
            openFolderAction.triggered.connect(self.openfolder)
            
            openTiffFolderAction = QtGui.QAction("Open Tiff Folder", self)
            openTiffFolderAction.triggered.connect(self.openTiffFolder)
            
            sinogramAction = QtGui.QAction('Sinogram', self)
            sinogramAction.triggered.connect(self.sinogram)

            saveImageAction = QtGui.QAction('Save Projections', self)
            saveImageAction.triggered.connect(self.saveImage)

            testAction =QtGui.QAction('test',self)
            testAction.triggered.connect(self.test)

            selectElementAction=QtGui.QAction('Select Element',self)
            selectElementAction.triggered.connect(self.selectElement)

            selectFilesAction=QtGui.QAction('Select Files',self)
            selectFilesAction.triggered.connect(self.selectFilesShow)

            saveThetaTxtAction = QtGui.QAction("Save Theta Postion as txt",self)
            saveThetaTxtAction.triggered.connect(self.saveThetaTxt)

            convertAction = QtGui.QAction('Save data in memory', self)
            convertAction.triggered.connect(self.convert)

            saveSinogramAction = QtGui.QAction('Save Sinogram',self)
            saveSinogramAction.triggered.connect(self.saveSinogram)

            runReconstructAction = QtGui.QAction("Reconstruction", self)
            runReconstructAction.triggered.connect(self.runReconstruct)

            selectImageTagAction = QtGui.QAction("Select Image Tag", self)
            selectImageTagAction.triggered.connect(self.selectImageTag)

            xCorAction = QtGui.QAction("Cross Correlation", self)
            xCorAction.triggered.connect(self.CrossCorrelation_test)

            phaseXCorAction = QtGui.QAction("Phase Correlation", self)
            phaseXCorAction.triggered.connect(self.CrossCorrelation_test)

            alignFromTextAction = QtGui.QAction("Alignment from Text", self)
            alignFromTextAction.triggered.connect(self.alignFromText)

            alignFromText2Action = QtGui.QAction("Alignment from Text2", self)
            alignFromText2Action.triggered.connect(self.alignFromText2)

            saveAlignToTextAction = QtGui.QAction("Save Alignment information to text", self)
            saveAlignToTextAction.triggered.connect(self.saveAlignToText)

            restoreAction = QtGui.QAction("Restore", self)
            restoreAction.triggered.connect(self.restore)

            readConfigAction = QtGui.QAction("Read configuration file",self)
            readConfigAction.triggered.connect(self.readConfigFile)

            runCenterOfMassAction = QtGui.QAction("run center of mass action",self)
            runCenterOfMassAction.triggered.connect(self.centerOfMassWindow)

            alignCenterOfMassAction = QtGui.QAction("Align by fitting center of mass position into sine curve",self)
            alignCenterOfMassAction.triggered.connect(self.alignCenterOfMass)

            matcherAction = QtGui.QAction("match template", self)
            matcherAction.triggered.connect(self.match_window)

            exportDataAction = QtGui.QAction("export data",self)
            exportDataAction.triggered.connect(self.export_data)

            runTransRecAction = QtGui.QAction("Transmission Recon", self)
            runTransRecAction.triggered.connect(self.runTransReconstruct)

            saveHotSpotPosAction = QtGui.QAction("Save Hot Spot Pos", self)
            saveHotSpotPosAction.triggered.connect(self.saveHotSpotPos)

            alignHotSpotPosAction = QtGui.QAction("Align Hot Spot pos", self)
            alignHotSpotPosAction.triggered.connect(self.alignHotSpotPos1)

            reorderAction = QtGui.QAction("Reorder", self)
            reorderAction.triggered.connect(self.reorder_matrix)

            wienerAction = QtGui.QAction("Wiener", self)
            wienerAction.triggered.connect(self.ipWiener)

            reorderAction = QtGui.QAction("Reorder", self)
            reorderAction.triggered.connect(self.reorder_matrix)

            externalImageRegAction= QtGui.QAction("External Image Registaration", self)
            externalImageRegAction.triggered.connect(self.externalImageReg)

            openPositionPanelAction = QtGui.QAction("Open Position Panel", self)
            openPositionPanelAction.triggered.connect(self.openPositionPanel)

            ###
            self.frame = QtGui.QFrame()
            self.vl = QtGui.QVBoxLayout()


            
            self.tab_widget = QtGui.QTabWidget()
            #self.tab_widget.addTab(self.createMessageWidget(), unicode("Message"))
            self.tab_widget.addTab(self.createImageProcessWidget(), unicode("Image Process"))
            self.tab_widget.addTab(self.createSaveHotspotWidget(),unicode("Alignment"))
            self.tab_widget.addTab(self.createSinoWidget(), unicode("Sinogram"))
            self.tab_widget.addTab(self.createReconWidget(), unicode("Reconstruction"))
            #self.tab_widget.addTab(self.sinoGroup, unicode("Sinogram"))



            self.tab_widget.currentChanged.connect(self.callW2)

            self.vl.addWidget(self.tab_widget)
            self.vl.addWidget(self.createMessageWidget())

##            self.win=pg.ImageView()
##            self.layout1 = QtGui.QVBoxLayout()
##            self.layout1.addWidget(self.win)
##            self.vl.addWidget(self.layout1)
            self.frame.setLayout(self.vl)
            self.setCentralWidget(self.frame)
            self.tab_widget.setDisabled(True)

            ###
            
##            self.lbl=QtGui.QLabel()
##            self.setCentralWidget(self.lbl)
##            self.lbl.setText("Starting")

            menubar = self.menuBar()
            fileMenu = menubar.addMenu('&File')
            fileMenu.addAction(openFileAction)
            fileMenu.addAction(openFolderAction)
            fileMenu.addAction(openTiffFolderAction)
            fileMenu.addAction(readConfigAction)
##            fileMenu.addAction(test1Action)
##            fileMenu.addAction(reorderAction)
            fileMenu.addAction(exitAction)
            fileMenu.addAction(closeAction)

            self.optionMenu = menubar.addMenu('Convert Option')
            self.optionMenu.addAction(convertAction)
            
            self.optionMenu.addAction(selectFilesAction)
            self.optionMenu.addAction(selectImageTagAction)
            self.optionMenu.addAction(openPositionPanelAction)
            self.optionMenu.setDisabled(True)

            self.alignmentMenu = menubar.addMenu("Alignment")
            self.alignmentMenu.addAction(saveAlignToTextAction)
            self.alignmentMenu.addAction(runCenterOfMassAction)
            self.alignmentMenu.addAction(alignCenterOfMassAction)
            self.alignmentMenu.addAction(xCorAction)
            self.alignmentMenu.addAction(phaseXCorAction)
            self.alignmentMenu.addAction(matcherAction)
            self.alignmentMenu.addAction(alignFromTextAction)
            self.alignmentMenu.addAction(alignFromText2Action)
            self.alignmentMenu.addAction(saveHotSpotPosAction)
            self.alignmentMenu.addAction(alignHotSpotPosAction)
            self.alignmentMenu.addAction(externalImageRegAction)
            self.alignmentMenu.addAction(restoreAction)
            self.alignmentMenu.setDisabled(True)

            self.afterConversionMenu = menubar.addMenu('After saving data in memory')
            self.afterConversionMenu.addAction(saveImageAction)
            self.afterConversionMenu.addAction(saveThetaTxtAction)
            #self.afterConversionMenu.addAction(selectElementAction)
            self.afterConversionMenu.addAction(saveSinogramAction)
            self.afterConversionMenu.addAction(runReconstructAction)
            self.afterConversionMenu.addAction(reorderAction)
            self.afterConversionMenu.setDisabled(True)

            toolbar = self.addToolBar('ToolBar')
            toolbar.addAction(exitAction)
            toolbar.addAction(openFileAction)
            toolbar.addAction(openFolderAction)
            toolbar.addAction(saveHotSpotPosAction)
            toolbar.addAction(alignHotSpotPosAction)
            toolbar.addAction(exportDataAction)
            toolbar.addAction(runTransRecAction)
            toolbar.addAction(runCenterOfMassAction)
            toolbar.addAction(matcherAction)
            toolbar.addAction(runReconstructAction)
            toolbar.addAction(selectElementAction)
            toolbar.addAction(convertAction)
            toolbar.addAction(saveSinogramAction)
            toolbar.setVisible(False)
  
            add =0
            if platform == "win32":
                add = 50
            self.setGeometry(add,add, 1100+add,500+add)
            self.setWindowTitle('Maps_To_Tomopy')    
            self.show()


#### TEST
      def externalImageReg(self):

            original_path=os.getcwd()
            fileName = QtGui.QFileDialog.getExistingDirectory(self, "Open Extension",
                        QtCore.QDir.currentPath())
            fileName=str(fileName)
            sys.path.append(os.path.abspath(fileName))

            import imageReg
            from imageReg import align
            x=imageReg.align()
            print x

            os.chdir(original_path)
      def callW2(self):
            self.manual=Manual()
            self.manual.show()
            if self.tab_widget.currentIndex()==1:
                  self.manual.setVisible(True)
            if self.tab_widget.currentIndex()!=1:
                  self.manual.setVisible(False)
#############################
## creating tab
      def createGridInputWidget(self):
            GridStartVal = '2'
            hBox = QtGui.QHBoxLayout()
            self.GridXIn = QtGui.QLineEdit()
            self.GridYIn = QtGui.QLineEdit()
            self.GridZIn = QtGui.QLineEdit()

            self.GridXIn.setText(GridStartVal)
            self.GridYIn.setText(GridStartVal)
            self.GridZIn.setText(GridStartVal)

            hBox.addWidget(QtGui.QLabel("X"))
            hBox.addWidget(self.GridXIn)
            hBox.addWidget(QtGui.QLabel("Y"))
            hBox.addWidget(self.GridYIn)
            hBox.addWidget(QtGui.QLabel("Z"))
            hBox.addWidget(self.GridZIn)

            gridGroup = QtGui.QGroupBox('Grid Size')
            gridGroup.setLayout(hBox)

            return gridGroup


      def createMessageWidget(self):
            GridStartVal = '2'
            vbox= QtGui.QVBoxLayout()
            hbox = QtGui.QHBoxLayout()
            hbox2=QtGui.QHBoxLayout()
            self.lbl = QtGui.QLineEdit()
            self.lbl.setText("Starting")
            self.lbl2= QtGui.QLineEdit()
            self.lbl2.setText(os.getcwd())
            self.directoryButton = QtGui.QPushButton("Change Directory")
            self.directoryButton.clicked.connect(self.changeDirectory)

            hbox.addWidget(QtGui.QLabel("Message"))
            hbox.addWidget(self.lbl)
            hbox2.addWidget(QtGui.QLabel("Set Directory"))
            hbox2.addWidget(self.lbl2)
            hbox2.addWidget(self.directoryButton)
            vbox.addLayout(hbox)
            vbox.addLayout(hbox2)
            messageGroup = QtGui.QGroupBox("Message Box")
            messageGroup.setLayout(vbox)

            return messageGroup

      def createSinoWidget(self):
            self.sino = QSelect2()
            self.sinoView = IView4()


            sinoBox = QtGui.QHBoxLayout()
            sinoBox.addWidget(self.sino)
            sinoBox.addWidget(self.sinoView,10)
            sinoGroup=QtGui.QGroupBox("Sinogram")
            sinoGroup.setLayout(sinoBox)

            return sinoGroup

      def createReconWidget(self):
            self.recon = QSelect3()
            self.recon.sld.setVisible(False)
            self.reconView = IView3()
            self.reconView.view.ROI.setVisible(False)

            reconBox = QtGui.QHBoxLayout()
            reconBox.addWidget(self.recon)
            reconBox.addWidget(self.reconView,10)
            reconGroup=QtGui.QGroupBox("Reconstruction")
            reconGroup.setLayout(reconBox)


            return reconGroup

      def createProjWidget(self):
            self.projection = QSelect2()
            self.projectionView =pg.ImageView()

            projectionBox = QtGui.QHBoxLayout()
            projectionBox.addWidget(self.projection)
            projectionBox.addWidget(self.projectionView,10)
            projectionGroup=QtGui.QGroupBox("Projections")
            projectionGroup.setLayout(projectionBox) 

            return projectionGroup

      def createSaveHotspotWidget(self):
            self.projViewControl = QSelect4()
            self.projView=IView3()
            self.boxSize=20

            projViewBox = QtGui.QHBoxLayout()
            projViewBox.addWidget(self.projViewControl)
            projViewBox.addWidget(self.projView,10)
            projViewGroup=QtGui.QGroupBox("Save HotSpot")
            projViewGroup.setLayout(projViewBox)
            return projViewGroup

      def createImageProcessWidget(self):
            self.imgProcessControl = imageProcess()
            self.imgProcess = IView3()

            imgProcessBox = QtGui.QHBoxLayout()
            imgProcessBox.addWidget(self.imgProcessControl)
            imgProcessBox.addWidget(self.imgProcess,10)
            imgProcessGroup=QtGui.QGroupBox("Image Process")
            imgProcessGroup.setLayout(imgProcessBox)
            return imgProcessGroup
            
#############################
#! Alignment
      def test(self):
            1

      def changeDirectory(self):
            folderName=QtGui.QFileDialog.getExistingDirectory(self,"Change Directory",
                                                                   QtCore.QDir.currentPath())
            os.chdir(str(folderName))
            self.lbl2.setText(folderName)
            
      def centerOfMassWindow(self):
            self.comer = QSelect3()
            self.comer.setWindowTitle("Center of Mass window")
            self.comer.numb= len(self.channelname)
            for j in arange(self.comer.numb):
                  self.comer.combo.addItem(self.channelname[j])
            self.comer.btn.setText("Center of Mass")
            self.comer.method.setVisible(False)
            self.comer.save.setVisible(True)
            self.comer.save.setText("Restore")
            self.comer.btn.clicked.connect(self.runCenterOfMass)
            self.comer.save.clicked.connect(self.restore)
            self.comer.show()

      def runCenterOfMass(self):
            self.centerOfMass()
            self.fitCenterOfMass(x=self.theta)
            self.lbl.setText("Center of Mass: "+str(self.p1[2]))
#### devel
      def runCenterOfMass2(self):
            self.com = zeros(self.projections)
            temp=zeros(self.data.shape[3])
            temp2=zeros(self.data.shape[3])
            self.comelem = self.sino.combo.currentIndex()
            for i in arange(self.projections):
                  temp=sum(self.data[self.comelem,i,self.sino.sld.value()-self.thickness/2:self.sino.sld.value()+self.thickness/2,:]-self.data[self.comelem,i,:10,:10].mean(), axis=0)
                  #temp=sum(self.data[self.comelem,i,:,:]-1, axis=0)
                  numb2=sum(temp)
                  for j in arange(self.data.shape[3]):
                        
                        temp2[j]=temp[j]*j
                  numb=float(sum(temp2))/numb2
                  self.com[i]=numb
            self.fitCenterOfMass(x=self.theta)
            self.lbl.setText("Center of Mass: "+str(self.p1[2]))
            self.alignCenterOfMass()
            self.sinogram()
            
#### devel
            
      def centerOfMass(self):
            self.com = zeros(self.projections)
            temp=zeros(self.data.shape[3])
            temp2=zeros(self.data.shape[3])
            self.comelem = self.comer.combo.currentIndex()
            for i in arange(self.projections):
                  temp=sum(self.data[self.comelem,i,:,:]-self.data[self.comelem,i,:10,:10].mean(), axis=0)
                  #temp=sum(self.data[self.comelem,i,:,:]-1, axis=0)
                  numb2=sum(temp)
                  for j in arange(self.data.shape[3]):
                        
                        temp2[j]=temp[j]*j
                  numb=float(sum(temp2))/numb2
                  self.com[i]=numb

            
      def fitCenterOfMass(self,x):
            self.fitfunc = lambda p,x: p[0]*sin(2*pi/360*(x-p[1]))+p[2]
            self.errfunc = lambda p,x,y: self.fitfunc(p,x)-y
            p0=[100,100,100]
            self.p1,success = optimize.leastsq(self.errfunc,p0[:],args=(x,self.com))
            self.centerOfMassDiff=self.fitfunc(self.p1,x)-self.com
            print "here", self.centerOfMassDiff

      def fitCenterOfMass2(self,x):
            fitfunc = lambda p,x: p[0]*sin(2*pi/360*(x-p[1]))+self.p1[2]
            errfunc = lambda p,x,y: fitfunc(p,x)-y
            p0=[100,100]
            p2,success = optimize.leastsq(errfunc,p0[:],args=(x,self.com))
            self.centerOfMassDiff=fitfunc(p2,x)-self.com
            print "here", self.centerOfMassDiff

      def alignCenterOfMass(self):
            for i in arange(self.projections):
                  self.xshift[i]+=int(self.centerOfMassDiff[i])
                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], int(round(self.xshift[i])), axis=2)
            self.lbl.setText("Alignment has been completed")
                  
      def alignCenterOfMass2(self):
            j=0

            for i in self.hotspotProj:

                  self.xshift[i]+=int(self.centerOfMassDiff[j])

                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], int(round(self.xshift[i])), axis=2)
                  j+=1
            self.lbl.setText("Alignment has been completed")
                  
      def CrossCorrelation_test(self):
            
            self.datacopy=zeros(self.data.shape)
            self.datacopy[...]=self.data[...]
            self.data[np.isnan(self.data)]=1
            self.xcor = AlignWindow()
            self.xcor.setWindowTitle("CrossCorrelation Window")
            self.xcor.numb = len(self.channelname)
            for j in arange(self.xcor.numb):
                  self.xcor.combo.addItem(self.channelname[j])
            self.xcor.btn.setText("Cross Correlation")
            self.xcor.btn2.setText("Restore")
            self.xcor.btn.clicked.connect(self.xCor)
            self.xcor.btn2.clicked.connect(self.restore)
            self.xcor.method.setVisible(False)
            self.xcor.show()
      def restore(self):
            self.xshift=zeros(self.projections,int)
            self.yshift=zeros(self.projections,int)
            self.data=zeros(self.oldData.shape)
            self.data[...]=self.oldData[...]
            self.projView.view.data=self.data[self.projViewElement,:,:,:]
      def prexCor(self):
            try:
                  self.xcor.savedir=QtGui.QFileDialog.getSaveFileName()
                  if not self.xcor.savedir[0]:
                        raise IndexError
                  self.xcor.savedir=self.xcor.savedir[0]
                  self.xCor()
            except IndexError:
                  print "type the header name"
      def match_window(self):
            self.matcher = AlignWindow()
            self.matcher.setWindowTitle("Match template window")
            self.matcher.numb= len(self.channelname)
            for j in arange(self.matcher.numb):
                  self.matcher.combo.addItem(self.channelname[j])
            self.matcher.btn.setText("Match Template")

            self.matcher.btn2.setText("Restore")
            self.matcher.btn.clicked.connect(self.match)
            self.matcher.btn2.clicked.connect(self.restore)
            self.xcor.method.setVisible(False)
            self.matcher.show()
      def match(self):
            self.matchElem=self.matcher.combo.currentIndex()
            for i in arange(self.projections-1):
                  img=self.data[self.matchElem,i,:,:]
                  img1=ones([img.shape[0]*2,img.shape[1]*2])*self.data[self.matchElem,i,:10,:10].mean()
                  img1[img.shape[0]/2:img.shape[0]*3/2,img.shape[1]/2:img.shape[1]*3/2]=img
                  img2=self.data[self.matchElem,i+1,:,:]
                  result=match_template(img1,img2)
                  result=np.where(result==np.max(result))
                  self.yshift[i+1]+=result[0][0]-img.shape[0]/2
                  self.xshift[i+1]+=result[1][0]-img.shape[1]/2
                  print self.xshift[i+1], self.yshift[i+1]
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], self.xshift[i+1], axis=2)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], self.yshift[i+1], axis=1)

      def xCor(self):
##            self.xcor.savedir="texting"
##            f=open(self.xcor.savedir+".txt",'w')
##            onlyfilename=self.fileNames[0].rfind("/")
##            f.write(self.fileNames[0][onlyfilename+1:]+" \n")
##            f.write("0 \n")
            self.xcorElement = self.xcor.combo.currentIndex()
            for i in arange(self.projections-1):
                  #onlyfilename=self.fileNames[i+1].rfind("/")
                  img1=self.data[self.xcorElement,i,:,:]
                  img2=self.data[self.xcorElement,i+1,:,:]
                  
                  self.t0,self.t1=self.xcorrelate(img1,img2)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:],self.t0,axis=1)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:],self.t1,axis=2)
                  self.xshift[i+1]+=self.t1
                  self.yshift[i+1]+=self.t0
                  

                  
##                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:],shift,axis=2)
##                  f.write(self.fileNames[i+1][onlyfilename+1:]+" \n")
##                  f.write(str(shift)+ "\n")
##                  print i
##            f.close()
      def alignFromText(self):
            try:
                  fileName = QtGui.QFileDialog.getOpenFileName(self, "Open File",
                        QtCore.QDir.currentPath(),"TXT (*.txt)")
                  ##### for future reference "All File (*);;CSV (*.csv *.CSV)"

                  f=open(fileName,'r')
                  read = f.readlines()
                  self.datacopy=zeros(self.data.shape)
                  self.datacopy[...]=self.data[...]
                  self.data[np.isnan(self.data)]=1
                  for i in arange(self.projections):
                        onlyfilename=self.selectedFiles[i].rfind("/")
                        for j in arange(len(read)):
                              if string.find(read[j],self.selectedFiles[i][onlyfilename+1:])!=-1:
                                    secondcol=read[j].rfind(",")
                                    firstcol=read[j][:secondcol].rfind(",")
                                    self.yshift[i]+=int(float(read[j][secondcol+1:-1]))
                                    self.xshift[i]+=int(float(read[j][firstcol+1:secondcol]))
                                    self.data[:,i,:,:]=np.roll(self.data[:,i,:,:],self.xshift[i],axis=2)
                                    self.data[:,i,:,:]=np.roll(self.data[:,i,:,:],self.yshift[i],axis=1)
                              if string.find(read[j],"rotation axis")!=-1:
                                    commapos=read[j].rfind(",")
                                    self.p1[2]=float(read[j][commapos+1:-1])
            
                  f.close()

                  self.lbl.setText("Alignment using values from Text has been completed")
                  self.updateImages()
            except IOError:
                  print "choose file please"
      def alignFromText2(self):
            try:
                  fileName = QtGui.QFileDialog.getOpenFileName(self, "Open File",
                        QtCore.QDir.currentPath(),"TXT (*.txt)")
                  ##### for future reference "All File (*);;CSV (*.csv *.CSV)"

                  f=open(fileName,'r')
                  read = f.readlines()
                  self.datacopy=zeros(self.data.shape)
                  self.datacopy[...]=self.data[...]
                  self.data[np.isnan(self.data)]=1
                  for i in arange(self.projections):
                        j=i+1
                        secondcol=read[j].rfind(",")
                        firstcol=read[j][:secondcol].rfind(",")
                        self.yshift[i]+=int(float(read[j][secondcol+1:-1]))
                        self.xshift[i]+=int(float(read[j][firstcol+1:secondcol]))
                        self.data[:,i,:,:]=np.roll(self.data[:,i,:,:],self.xshift[i],axis=2)
                        self.data[:,i,:,:]=np.roll(self.data[:,i,:,:],self.yshift[i],axis=1)

            
                  f.close()

                  self.lbl.setText("Alignment using values from Text has been completed")
                  self.updateImages()
            except IOError:
                  print "choose file please"

      def saveAlignToText(self):

            try:
                  self.alignFileName = QtGui.QFileDialog.getSaveFileName()
                  if string.rfind(str(self.alignFileName),".txt")==-1:
                        self.alignFileName=str(self.alignFileName)+".txt"
                  print str(self.alignFileName)
                  f=open(self.alignFileName,"w")
                  f.writelines("rotation axis, "+ str(self.p1[2])+"\n")
                  for i in arange(self.projections):
                        onlyfilename=self.selectedFiles[i].rfind("/")
                        print self.selectedFiles[i]
                        f.writelines(self.selectedFiles[i][onlyfilename+1:]+", "+str(self.xshift[i])+", "+str(self.yshift[i])+ "\n")
                        
                  f.close()
            except IOError:
                  print "choose file please"
                  

                        

#==========================
      def xcorrelate(self,a,b):
            fa = spf.fft2(a)
            fb = spf.fft2(b)

            shape = a.shape
            c = abs(spf.ifft2(fa*fb.conjugate()))
            t0, t1 = np.unravel_index(np.argmax(c), a.shape)
            if t0 > shape[0]//2:
                  t0 -= shape[0]
            if t1 > shape[1]//2:
                  t1 -= shape[1]
            
            return t0,t1

      def phasecorrelate(self,a,b):
            fa = spf.fft2(a)
            fb = spf.fft2(b)

            shape = a.shape
            c = abs(spf.ifft2(fa*fb.conjugate()/(abs(fa)*abs(fb))))
            t0, t1 = np.unravel_index(np.argmax(c), a.shape)
            if t0 > shape[0]//2:
                  t0 -= shape[0]
            if t1 > shape[1]//2:
                  t1 -= shape[1]
            return t0,t1
      
      def edgegauss(self,imagey,sigma=4):
            image=zeros(imagey.shape)
            image[...]=imagey[...]
            nx=image.shape[1]
            ny=image.shape[0]

            n_sigma = -log(10**-6)
            n_roll = max(int(1+sigma*n_sigma),2)
            exparg = float32(arange(n_roll)/float(sigma))
            rolloff = float32(1)-exp(-0.5*exparg*exparg)

            ## Top edge

            xstart = 0
            xstop = nx 
            iy= 0

            for i_roll in arange(n_roll):
                  image[iy,xstart:xstop] = image[iy, xstart:xstop] * rolloff[iy]
                  xstart = min(xstart+1,nx/2-1)
                  xstop = max(xstop - 1,nx/2)
                  iy = min(iy+1,ny-1)


            ## Bottom edge

            xstart = 0
            xstop = nx
            iy = ny-1

            for i_roll in arange(n_roll):
                  image[iy,xstart:xstop] = image[iy,xstart:xstop]*rolloff[ny-1-iy]
                  xstart = min(xstart +1,nx/2-1)
                  xstop = max(xstop-1,nx/2)
                  iy = max(iy-1, 0)

            ## Left edge

            ystart = 1
            ystop = ny - 1
            ix = 0

            for i_roll in arange(n_roll):
                  image[ystart:ystop,ix]=image[ystart:ystop,ix] * rolloff[ix]
                  ystart = min(ystart+1,ny/2-1)
                  ystop = max(ystop-1 , ny/2)
                  ix = min(ix+1,nx-1)

            ## Right edge

            ystart = 1
            ystop = ny-1
            ix = nx-1

            for i_roll in arange(n_roll):
                  image[ystart:ystop, ix] = image[ystart:ystop,ix] * rolloff[nx-1-ix]
                  ystart = min(ystart+1, ny/2-1)
                  ystop = max(ystop-1, ny/2)
                  ix = max(ix - 1,0)
                  
            return image

#==========================

      def showSaveHotSpotPos(self):
            self.tab_widget.removeTab(1)
            self.tab_widget.insertTab(1,self.createSaveHotspotWidget(),unicode("Alignment"))
            self.projViewControl.numb=len(self.channelname)
            for j in arange(self.projViewControl.numb):
                  self.projViewControl.combo.addItem(self.channelname[j])

            for k in arange(self.projections):
                  self.projViewControl.combo3.addItem(str(k+1))

            self.projViewControl.combo3.setVisible(False)

            
            self.projViewControl.combo.currentIndexChanged.connect(self.saveHotSpotPos)
            self.projViewControl.combo3.currentIndexChanged.connect(self.hotspotProjChanged)
            self.projViewControl.sld.setValue(20)
            self.projViewControl.sld.setRange(0,self.x/2)
            self.projViewControl.lcd.display(20)
            self.projViewControl.sld.valueChanged.connect(self.projViewControl.lcd.display)
            self.projViewControl.sld.valueChanged.connect(self.boxSizeChange)
            self.projViewControl.btn.clicked.connect(self.alignHotSpotPos3)
            self.projViewControl.btn2.clicked.connect(self.alignHotSpotPos4)
            self.projViewControl.btn3.clicked.connect(self.alignHotSpotY)
            self.projViewControl.btn4.clicked.connect(self.clearHotSpotData)
            self.projViewControl.combo2.currentIndexChanged.connect(self.hotSpotSetChanged)
            self.projViewControl.show()

            self.projView.view.hotSpotNumb=0
            self.projView.sld.setRange(0,self.projections-1)
            self.projView.sld.valueChanged.connect(self.projView.lcd.display)
            self.projView.sld.valueChanged.connect(self.hotSpotProjChanged)
            self.testtest=pg.ImageView()
#########

      def clearHotSpotData(self):
            self.projView.view.posMat[...]=zeros_like(self.projView.view.posMat)
      def hotspotProjChanged(self):
            1
      def hotSpotProjChanged(self):
            self.projView.view.hotSpotNumb=self.projView.sld.value()
            self.projView.view.projView.setImage(self.data[self.projViewElement,self.projView.view.hotSpotNumb,:,:])

#########
      
      def boxSizeChange(self):
            self.boxSize=self.projViewControl.sld.value()/2*2
            self.projView.view.ROI.setPos([int(round(self.projView.view.projView.iniX))-self.boxSize/2,-int(round(self.projView.view.projView.iniY))-self.boxSize/2])
            self.projView.view.ROI.setSize([self.boxSize,self.boxSize])
            self.projView.view.xSize=self.boxSize
            self.projView.view.ySize=self.boxSize
      def hotSpotSetChanged(self):
            self.projView.view.hotSpotSetNumb=self.projViewControl.combo2.currentIndex()

      def nextHotSpotPos(self):
            #self.projView.hotSpotNumb=self.projViewControl.sld.value()
            self.projView.view.projView.setImage(self.data[self.projViewElement,self.projView.view.hotSpotNumb,:,:])

            
      def saveHotSpotPos(self):
            #self.projView.view.hotSpotNumb=0
            self.projViewElement = self.projViewControl.combo.currentIndex()
            self.projView.view.data=self.data[self.projViewElement,:,:,:]
            self.projView.view.posMat=zeros([5,self.data.shape[1],2]) ## Later change 5 -> how many data are in the combo box.
            self.projView.view.projView.setImage(self.data[self.projViewElement,0,:,:])
      def alignHotSpotY(self):
            self.boxSize2=self.boxSize/2
            self.xPos=zeros(self.projections)
            self.yPos=zeros(self.projections)
            self.boxPos=zeros([self.projections,self.boxSize,self.boxSize])
            hotSpotSet=self.projViewControl.combo2.currentIndex()
            for i in arange(self.projections):
                  
                  self.yPos[i]=int(round(self.projView.view.posMat[hotSpotSet,i,0]))
                  self.xPos[i]=int(round(self.projView.view.posMat[hotSpotSet,i,1]))
            
                  if self.xPos[i]!=0 and self.yPos[i]!=0:
                        if self.yPos[i]<self.boxSize2:
                              self.yPos[i]=self.boxSize2
                        if self.yPos[i]>self.projView.view.data.shape[1]-self.boxSize2:
                              self.yPos[i]=self.projView.view.data.shape[1]-self.boxSize2
                        if self.xPos[i]<self.boxSize2:
                              self.xPos[i]=self.boxSize2
                        if self.xPos[i]>self.projView.view.data.shape[2]-self.boxSize2:
                              self.xPos[i]=self.projView.view.data.shape[2]-self.boxSize2
                        #self.boxPos[i,:,:]=self.projView.data[i,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2]
                        self.boxPos[i,:,:]=self.projView.view.data[i,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2]
            print self.boxPos.shape
            print self.xPos, self.yPos

            self.alignHotSpotY_next()

      def alignHotSpotY_next(self):
            self.hotSpotX=zeros(self.projections)
            self.hotSpotY=zeros(self.projections)
            self.newBoxPos=zeros(self.boxPos.shape)
            self.newBoxPos[...]=self.boxPos[...]
            ### need to change x and y's
            firstPosOfHotSpot=0
            add=1
            for i in arange(self.projections):
                  if self.xPos[i]==0 and self.yPos[i]==0:
                        firstPosOfHotSpot += add
                  if self.xPos[i]!=0 or self.yPos[i]!=0:
                        print self.xPos[i], self.yPos[i]
                        img=self.boxPos[i,:,:]
                        print img.shape
                        a,x,y,b,c=self.fitgaussian(img)
                        self.hotSpotY[i]=x
                        self.hotSpotX[i]=y
                        yshift=int(round(self.boxSize2-self.hotSpotY[i])) 
                        xshift=int(round(self.boxSize2-self.hotSpotX[i]))
                        self.newBoxPos[i,:,:]=np.roll(self.newBoxPos[i,:,:],xshift,axis=1)
                        self.newBoxPos[i,:,:]=np.roll(self.newBoxPos[i,:,:],yshift ,axis=0)
      ##                  subplot(211)
      ##                  plt.imshow(self.boxPos[i,:,:])
      ##                  subplot(212)
      ##                  plt.imshow(self.newBoxPos[i,:,:])
      ##                  show()
                        add=0


            add2=0
            for j in arange(self.projections):
                        
                  if self.xPos[j]!=0 and self.yPos[j]!=0:
                        yyshift=int(round(self.boxSize2-self.hotSpotY[j]-self.yPos[j]+self.yPos[firstPosOfHotSpot]))

                        print yyshift
                        self.data[:,j,:,:]=np.roll(self.data[:,j,:,:],
                                                   yyshift,axis=1)
##                        for l in arange(self.data.shape[0]):
##                              if yyshift>0:
##                                    self.data[l,j,:yyshift,:]=ones(self.data[l,j,:yyshift,:].shape)*self.data[l,j,:yyshift,:].mean()/2
##                              if yyshift<0:
##                                    self.data[l,j,yyshift:,:]=ones(self.data[l,j,yyshift:,:].shape)*self.data[l,j,-yyshift:,:].mean()/2
##                              if xxshift>0:
##                                    self.data[l,j,:,:xxshift]=ones(self.data[l,j,:,:xxshift].shape)*self.data[l,j,:xxshift,:].mean()/2
##                              if xxshift<0:
##                                    self.data[l,j,:,xxshift:]=ones(self.data[l,j,:,xxshift:].shape)*self.data[l,j,-xxshift:,:].mean()/2

                  if self.yPos[j]==0:
                        yyshift=0


                  self.yshift[j]+=yyshift

            print "align done"
            
      def alignHotSpotPos1(self):
            
            for i in arange(self.projections-1):
                  print "shifted"
                  
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], int(round(self.projView.view.posMat[0,0]-self.projView.view.posMat[i+1,0])), axis=2)
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], int(round(self.projView.view.posMat[0,1]-self.projView.view.posMat[i+1,1])), axis=1)            

      def alignHotSpotPos2(self):
            for i in arange(self.projections-1):
                  self.data[:,i+1,:,:]=np.roll(self.data[:,i+1,:,:], int(round(self.projView.view.posMat[0,1]-self.projView.view.posMat[i+1,1])), axis=1)     
            x=self.theta
            self.fitfunc = lambda p,x: p[0]*sin(2*pi/360*(x-p[1]))+p[2]
            self.errfunc = lambda p,x,y: self.fitfunc(p,x)-y
            p0=[100,100,100]
            self.p1,success = optimize.leastsq(self.errfunc,p0[:],args=(x,self.projView.view.posMat[:,0]))
            self.hotSpotPosDiff=self.fitfunc(self.p1,x)-self.projView.view.posMat[:,0]

            for i in arange(self.projections):
                  self.xshift[i]+=int(self.hotSpotPosDiff[i])
                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], self.xshift[i], axis=2)
            self.lbl.setText("Alignment has been completed")
            
      def alignHotSpotPos3(self):
            #self.projView.data2=self.data[7,:,:,:]
            self.boxSize2=self.boxSize/2
            self.xPos=zeros(self.projections)
            self.yPos=zeros(self.projections)
            self.boxPos=zeros([self.projections,self.boxSize,self.boxSize])
            hotSpotSet=self.projViewControl.combo2.currentIndex()
            for i in arange(self.projections):
                  
                  self.yPos[i]=int(round(self.projView.view.posMat[hotSpotSet,i,0]))
                  self.xPos[i]=int(round(self.projView.view.posMat[hotSpotSet,i,1]))
            
                  if self.xPos[i]!=0 and self.yPos[i]!=0:
                        if self.yPos[i]<self.boxSize2:
                              self.yPos[i]=self.boxSize2
                        if self.yPos[i]>self.projView.view.data.shape[1]-self.boxSize2:
                              self.yPos[i]=self.projView.view.data.shape[1]-self.boxSize2
                        if self.xPos[i]<self.boxSize2:
                              self.xPos[i]=self.boxSize2
                        if self.xPos[i]>self.projView.view.data.shape[2]-self.boxSize2:
                              self.xPos[i]=self.projView.view.data.shape[2]-self.boxSize2
                        #self.boxPos[i,:,:]=self.projView.data[i,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2]
                        self.boxPos[i,:,:]=self.projView.view.data[i,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2]
            print self.boxPos.shape
            print "x", self.xPos
            print "y", self.yPos


##            for i in arange(self.projections):
##                  j=Image.fromarray(self.boxPos[i,:,:].astype(np.float32))
##                  
##                  j.save("/Users/youngpyohong/Documents/Work/Python/2dfit/"+str(i)+".tif")


            self.alignHotSpotPos3_3()
            print "hotspot done"

      def alignHotSpotPos4(self):
            #self.projView.data2=self.data[7,:,:,:]
            self.boxSize2=self.boxSize/2
            self.xPos=zeros(self.projections)
            self.yPos=zeros(self.projections)
            self.boxPos=zeros([self.projections,self.boxSize,self.boxSize])
            hotSpotSet=self.projViewControl.combo2.currentIndex()
            for i in arange(self.projections):
                  
                  self.yPos[i]=int(round(self.projView.view.posMat[hotSpotSet,i,0]))
                  self.xPos[i]=int(round(self.projView.view.posMat[hotSpotSet,i,1]))
            
                  if self.xPos[i]!=0 and self.yPos[i]!=0:
                        if self.yPos[i]<self.boxSize2:
                              self.yPos[i]=self.boxSize2
                        if self.yPos[i]>self.projView.view.data.shape[1]-self.boxSize2:
                              self.yPos[i]=self.projView.view.data.shape[1]-self.boxSize2
                        if self.xPos[i]<self.boxSize2:
                              self.xPos[i]=self.boxSize2
                        if self.xPos[i]>self.projView.view.data.shape[2]-self.boxSize2:
                              self.xPos[i]=self.projView.view.data.shape[2]-self.boxSize2
                        #self.boxPos[i,:,:]=self.projView.data[i,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2]
                        self.boxPos[i,:,:]=self.projView.view.data[i,self.yPos[i]-self.boxSize2:self.yPos[i]+self.boxSize2,self.xPos[i]-self.boxSize2:self.xPos[i]+self.boxSize2]
            print self.boxPos.shape
            print self.xPos, self.yPos


##            for i in arange(self.projections):
##                  j=Image.fromarray(self.boxPos[i,:,:].astype(np.float32))
##                  
##                  j.save("/Users/youngpyohong/Documents/Work/Python/2dfit/"+str(i)+".tif")

            self.alignHotSpotPos3_4()
            print "hotspot done"



 


      def alignHotSpotPos3_1(self):
            self.newBoxPos=zeros(self.boxPos.shape)
            self.newBoxPos[0,:,:]=self.boxPos[0,:,:]
            self.boxShift=zeros([2,self.projections])
            for i in arange(self.projections-1):
                  #t0,t1=self.phasecorrelate(self.boxPos[i,:,:],self.boxPos[i+1,:,:])
                  t0,t1=self.xcorrelate(self.boxPos[i,:,:],self.boxPos[i+1,:,:])
                  print t0,t1
                  self.newBoxPos[i+1,:,:]=np.roll(np.roll(self.boxPos[i+1,:,:],t1,axis=1),t0,axis=0)
                  self.boxShift[0,i+1]=t0
                  self.boxShift[1,i+1]=t1
            self.testtest=pg.ImageView()
            self.testtest.setWindowTitle("rhyme")
            self.testtest.setImage(self.newBoxPos)
            self.testtest.show()

            for j in arange(self.projections-1):
                  self.data[:,j,:,:]=np.roll(np.roll(self.data[:,j,:,:],int(round(self.boxShift[1,j+1]-self.xPos[j+1]+self.xPos[0])),axis=2),int(round(self.boxShift[0,j+1]-self.yPos[j+1]+self.yPos[0])),axis=1)
            self.test2=pg.ImageView()
            self.test2.setWindowTitle("test2")
            self.test2.setImage(self.data[self.projViewElement,:,:,:])
            self.test2.show()

      def alignHotSpotPos3_2(self):
            self.newBoxPos=zeros(self.boxPos.shape)
            self.newBoxPos[0,:,:]=self.boxPos[0,:,:]
            self.boxShift=zeros([2,self.projections])
            for i in arange(self.projections-1):
                  
                  img1=self.boxPos[i,:,:]
                  img2=self.boxPos[i+1,:,:]
                  
                  corr=scipy.signal.correlate2d(img1,img2,fillvalue=np.average(img1[0:2,0:2]))
                  t0, t1 = np.unravel_index(np.argmax(corr), corr.shape)
                  self.newBoxPos[i+1,:,:]=np.roll(self.boxPos[i+1,:,:],40-t0,axis=0)
                  self.newBoxPos[i+1,:,:]=np.roll(self.boxPos[i+1,:,:],40-t1,axis=1)
                  print t0, t1
##            self.testtest=pg.ImageView()
##            self.testtest.setWindowTitle("rhyme")
##            self.testtest.setImage(self.newBoxPos)
##            self.testtest.show()

      def alignHotSpotPos3_3(self):
            self.hotSpotX=zeros(self.projections)
            self.hotSpotY=zeros(self.projections)
            self.newBoxPos=zeros(self.boxPos.shape)
            self.newBoxPos[...]=self.boxPos[...]
            ### need to change x and y's
            firstPosOfHotSpot=0
            add=1
            for i in arange(self.projections):
                  if self.xPos[i]==0 and self.yPos[i]==0:
                        firstPosOfHotSpot += add
                  if self.xPos[i]!=0 or self.yPos[i]!=0:
                        print self.xPos[i], self.yPos[i]
                        img=self.boxPos[i,:,:]
                        print img.shape
                        a,x,y,b,c=self.fitgaussian(img)
                        self.hotSpotY[i]=x
                        self.hotSpotX[i]=y
                        yshift=int(round(self.boxSize2-self.hotSpotY[i])) 
                        xshift=int(round(self.boxSize2-self.hotSpotX[i]))
                        self.newBoxPos[i,:,:]=np.roll(self.newBoxPos[i,:,:],xshift,axis=1)
                        self.newBoxPos[i,:,:]=np.roll(self.newBoxPos[i,:,:],yshift ,axis=0)
      ##                  subplot(211)
      ##                  plt.imshow(self.boxPos[i,:,:])
      ##                  subplot(212)
      ##                  plt.imshow(self.newBoxPos[i,:,:])
      ##                  show()
                        add=0

            add2=0
            for j in arange(self.projections):
                        
                  if self.xPos[j]!=0 and self.yPos[j]!=0:
                        yyshift=int(round(self.boxSize2-self.hotSpotY[j]-self.yPos[j]+self.yPos[firstPosOfHotSpot]))
                        xxshift=int(round(self.boxSize2-self.hotSpotX[j]-self.xPos[j]+self.xPos[firstPosOfHotSpot]))
                        print xxshift, yyshift
                        self.data[:,j,:,:]=np.roll(np.roll(self.data[:,j,:,:],xxshift,axis=2),
                                                   yyshift,axis=1)
##                        for l in arange(self.data.shape[0]):
##                              if yyshift>0:
##                                    self.data[l,j,:yyshift,:]=ones(self.data[l,j,:yyshift,:].shape)*self.data[l,j,:yyshift,:].mean()/2
##                              if yyshift<0:
##                                    self.data[l,j,yyshift:,:]=ones(self.data[l,j,yyshift:,:].shape)*self.data[l,j,-yyshift:,:].mean()/2
##                              if xxshift>0:
##                                    self.data[l,j,:,:xxshift]=ones(self.data[l,j,:,:xxshift].shape)*self.data[l,j,:xxshift,:].mean()/2
##                              if xxshift<0:
##                                    self.data[l,j,:,xxshift:]=ones(self.data[l,j,:,xxshift:].shape)*self.data[l,j,-xxshift:,:].mean()/2
                  if self.xPos[j]==0:
                        xxshift=0
                  if self.yPos[j]==0:
                        yyshift=0

                  self.xshift[j]+=xxshift
                  self.yshift[j]+=yyshift
                                                                                
                  
            self.p1[2]=self.xPos[0]

            print "align done"

      def alignHotSpotPos3_4(self):
            self.hotSpotX=zeros(self.projections)
            self.hotSpotY=zeros(self.projections)
            self.newBoxPos=zeros(self.boxPos.shape)
            self.newBoxPos[...]=self.boxPos[...]
            ### need to change x and y's
            firstPosOfHotSpot=0
            add=1
            for i in arange(self.projections):
                  if self.xPos[i]==0 and self.yPos[i]==0:
                        firstPosOfHotSpot += add
                  if self.xPos[i]!=0 or self.yPos[i]!=0:
                        print self.xPos[i], self.yPos[i]
                        img=self.boxPos[i,:,:]
                        print img.shape
                        a,x,y,b,c=self.fitgaussian(img)
                        self.hotSpotY[i]=x
                        self.hotSpotX[i]=y
                        yshift=int(round(self.boxSize2-self.hotSpotY[i])) 
                        xshift=int(round(self.boxSize2-self.hotSpotX[i]))
                        self.newBoxPos[i,:,:]=np.roll(self.newBoxPos[i,:,:],xshift,axis=1)
                        self.newBoxPos[i,:,:]=np.roll(self.newBoxPos[i,:,:],yshift ,axis=0)
      ##                  subplot(211)
      ##                  plt.imshow(self.boxPos[i,:,:])
      ##                  subplot(212)
      ##                  plt.imshow(self.newBoxPos[i,:,:])
      ##                  show()
                        add=0

            for j in arange(self.projections):
                        
                  if self.xPos[j]!=0 and self.yPos[j]!=0:
                        yyshift=int(round(self.boxSize2-self.hotSpotY[j]))
                        xxshift=int(round(self.boxSize2-self.hotSpotX[j]))

##                        for l in arange(self.data.shape[0]):
##                              if yyshift>0:
##                                    self.data[l,j,:yyshift,:]=ones(self.data[l,j,:yyshift,:].shape)*self.data[l,j,:yyshift,:].mean()/2
##                              if yyshift<0:
##                                    self.data[l,j,yyshift:,:]=ones(self.data[l,j,yyshift:,:].shape)*self.data[l,j,-yyshift:,:].mean()/2
##                              if xxshift>0:
##                                    self.data[l,j,:,:xxshift]=ones(self.data[l,j,:,:xxshift].shape)*self.data[l,j,:xxshift,:].mean()/2
##                              if xxshift<0:
##                                    self.data[l,j,:,xxshift:]=ones(self.data[l,j,:,xxshift:].shape)*self.data[l,j,-xxshift:,:].mean()/2
                  if self.xPos[j]==0:
                        xxshift=0
                  if self.yPos[j]==0:
                        yyshift=0

                  self.xshift[j]+=xxshift
                  self.yshift[j]+=yyshift

            add2=0

            global hotspotXPos, hotspotYPos
            hotspotXPos=zeros(self.projections)
            hotspotYPos=zeros(self.projections)
            for i in arange(self.projections):
                  hotspotYPos[i]=int(round(self.yPos[i]))
                  hotspotXPos[i]=int(round(self.xPos[i]))
            self.hotspotProj = np.where(hotspotXPos!=0)[0]
            print self.hotspotProj
            ## temp


            ## xfit
            print self.hotspotProj
            global a1,b4
            a1=self.theta
            b4=self.hotspotProj
            theta=self.theta[self.hotspotProj]
            print "theta", theta
            self.com=hotspotXPos[self.hotspotProj]
            if self.projViewControl.combo2.currentIndex()==0:
                  self.fitCenterOfMass(x=theta)
            else:
                  self.fitCenterOfMass2(x=theta)
            self.alignCenterOfMass2()

            ## yfit
            for i in self.hotspotProj:
                  self.yshift[i]+=int(hotspotYPos[self.hotspotProj[0]])-int(hotspotYPos[i])
                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], self.yshift[i],axis=1)
                  print int(hotspotYPos[0])-int(hotspotYPos[i])

            self.recon.sld.setValue(self.p1[2])
            print "align done"
            
#==========================
## This is for Image processing tab
      def showImageProcess(self):
            self.tab_widget.removeTab(0)
            self.tab_widget.insertTab(0,self.createImageProcessWidget(),unicode("Image Process"))
            self.imgProcessControl.numb=len(self.channelname)
            for j in arange(self.imgProcessControl.numb):
                  self.imgProcessControl.combo1.addItem(self.channelname[j])

            for k in arange(self.projections):
                  self.imgProcessControl.combo2.addItem(str(k+1))
            
            self.imgProcessControl.combo1.currentIndexChanged.connect(self.imgProcessProjShow)
            self.imgProcessControl.combo2.currentIndexChanged.connect(self.imgProcessProjShow)
            self.imgProcessControl.xUpBtn.clicked.connect(self.imgProcessBoxSizeChange)
            self.imgProcessControl.xDownBtn.clicked.connect(self.imgProcessBoxSizeChange)
            self.imgProcessControl.yUpBtn.clicked.connect(self.imgProcessBoxSizeChange)
            self.imgProcessControl.yDownBtn.clicked.connect(self.imgProcessBoxSizeChange)

            self.imgProcessControl.combo2.setVisible(False)
            

            self.imgProcessControl.bgBtn.clicked.connect(self.ipBg)
            self.imgProcessControl.delHotspotBtn.clicked.connect(self.ipDelHotspot)
            self.imgProcessControl.normalizeBtn.clicked.connect(self.ipNormalize)
            self.imgProcessControl.cutBtn.clicked.connect(self.ipCut)
            self.imgProcessControl.gaussian33Btn.clicked.connect(self.gauss33)
            self.imgProcessControl.gaussian33Btn.clicked.connect(self.gauss55)

            self.imgProcess.sld.setRange(0,self.projections-1)
            self.imgProcess.sld.valueChanged.connect(self.imgProcess.lcd.display)
            self.imgProcess.sld.valueChanged.connect(self.imgProcessProjChanged)
            self.testtest=pg.ImageView()
            
##            self.projViewControl.sld.setValue(20)
##            self.projViewControl.sld.setRange(0,self.x/2)
##            self.projViewControl.lcd.display(20)
##            self.projViewControl.sld.valueChanged.connect(self.projViewControl.lcd.display)
##            self.projViewControl.sld.valueChanged.connect(self.boxSizeChange)
##            self.projViewControl.btn.clicked.connect(self.alignHotSpotPos3)
##            self.projViewControl.btn2.clicked.connect(self.alignHotSpotPos4)
##            self.projViewControl.combo2.currentIndexChanged.connect(self.hotSpotSetChanged)
##            self.projViewControl.show()
##            self.testtest=pg.ImageView()

      def imgProcessProjChanged(self):
            element=self.imgProcessControl.combo1.currentIndex()
            self.imgProcessImg=self.data[element, self.imgProcess.sld.value(), : ,:]
            self.imgProcess.view.projView.setImage(self.imgProcessImg)
            
      def imgProcessProjShow(self):
            element=self.imgProcessControl.combo1.currentIndex()
            projection = self.imgProcessControl.combo2.currentIndex()
            self.imgProcessImg=self.data[element, projection, : ,:]
            self.imgProcess.view.projView.setImage(self.imgProcessImg)


      def imgProcessBoxSizeChange(self):
            xSize=self.imgProcessControl.xSize/2*2
            ySize=self.imgProcessControl.ySize/2*2
            self.imgProcess.view.ROI.setSize([xSize, ySize])
            self.imgProcess.view.ROI.setPos([int(round(self.imgProcess.view.projView.iniX))-xSize/2,-int(round(self.imgProcess.view.projView.iniY))-ySize/2])
            self.imgProcess.view.xSize=xSize
            self.imgProcess.view.ySize=ySize

      def ipBg(self):
            element=self.imgProcessControl.combo1.currentIndex()
            projection = self.imgProcessControl.combo2.currentIndex()
            xSize=self.imgProcessControl.xSize
            ySize=self.imgProcessControl.ySize
                  
            img=self.data[element,projection, int(round(self.imgProcess.view.projView.iniX))-xSize/2:int(round(self.imgProcess.view.projView.iniX))+xSize/2,
                          int(round(self.imgProcess.view.projView.iniY))-ySize/2:int(round(self.imgProcess.view.projView.iniY))+ySize/2]
            #if self.imgProcess.projView.iniX-xSize/2<##$#$#$#$#$#$#$#$
            self.bg = np.average(img)
            print self.bg
      def ipDelHotspot(self):
            element=self.imgProcessControl.combo1.currentIndex()
            projection = self.imgProcessControl.combo2.currentIndex()
            xSize=self.imgProcessControl.xSize
            ySize=self.imgProcessControl.ySize
            img=self.data[element,projection, int(round(self.imgProcess.view.projView.iniX))-xSize/2:int(round(self.imgProcess.view.projView.iniX))+xSize/2,
                          int(round(self.imgProcess.view.projView.iniY))-ySize/2:int(round(self.imgProcess.view.projView.iniY))+ySize/2]
            self.data[element,projection, int(round(self.imgProcess.view.projView.iniX))-xSize/2:int(round(self.imgProcess.view.projView.iniX))+xSize/2,
                          int(round(self.imgProcess.view.projView.iniY))-ySize/2:int(round(self.imgProcess.view.projView.iniY))+ySize/2]=ones(img.shape,dtype=img.dtype)*self.bg

            self.imgProcess.view.projView.setImage(self.data[element,projection,:,:])
      def ipNormalize(self):
            element=self.imgProcessControl.combo1.currentIndex()
            projection = self.imgProcessControl.combo2.currentIndex()
            normData = self.data[element,:,:,:]
            for i in arange(normData.shape[0]):
                  temp=normData[i,:,:][...]
                  tempMax=temp.max()
                  tempMin=temp.min()
                  temp=(temp-tempMin)/tempMax*1000
                  self.data[element,i,:,:]=temp


      def ipCut(self):
            element=self.imgProcessControl.combo1.currentIndex()
            projection = self.imgProcessControl.combo2.currentIndex()
            xSize=self.imgProcessControl.xSize
            ySize=self.imgProcessControl.ySize
            print xSize, ySize
            img=self.data[element,projection, int(round(self.imgProcess.view.projView.iniY))-ySize/2:int(round(self.imgProcess.view.projView.iniY))+ySize/2,
                          int(round(self.imgProcess.view.projView.iniX))-xSize/2:int(round(self.imgProcess.view.projView.iniX))+xSize/2]
            print img.shape,round(self.imgProcess.view.projView.iniX),round(self.imgProcess.view.projView.iniY)
            self.temp_data=zeros([len(self.channelname),self.projections,img.shape[0],img.shape[1]])
            print self.data.shape
            for i in arange(self.projections):
                  for j in arange(len(self.channelname)):

                        self.temp_data[j,i,:,:]= self.data[j,i, int(round(self.imgProcess.view.projView.iniY))-ySize/2:int(round(self.imgProcess.view.projView.iniY))+ySize/2,
                                int(round(self.imgProcess.view.projView.iniX))-xSize/2:int(round(self.imgProcess.view.projView.iniX))+xSize/2]
            print "done"
            self.data=self.temp_data

      def ipWiener(self):
            element=self.imgProcessControl.combo1.currentIndex()
            for i in arange(self.projections):
                  img=self.data[element, i,:,:]
                  img[np.where(img==0)]=10**(-8)
                  img_fin = wiener(img)
                  self.data[element, i,:,:]=img_fin


      def gauss2D(self,shape=(3,3),sigma=0.5):
            """
            2D gaussian mask - should give the same result as MATLAB's
            fspecial('gaussian',[shape],[sigma])
            """
            m,n = [(ss-1.)/2. for ss in shape]
            y,x = np.ogrid[-m:m+1,-n:n+1]
            h = np.exp( -(x*x + y*y) / (2.*sigma*sigma) )
            h[ h < np.finfo(h.dtype).eps*h.max() ] = 0
            sumh = h.sum()
            if sumh != 0:
              h /= sumh
            return h
      def gauss33(self):
            result=self.gauss2D(shape=(3,3),sigma=1.3)
            return result
      def gauss55(self):
            result=self.gauss2D(shape=(5,5),sigma=1.3)
            return result
            
#==========================
## Gaussian fit from wiki.scipy.org/Cookbook/FittingData
      def gaussian(self,height, center_x, center_y, width_x, width_y):
            """Returns a gaussian function with the given parameters"""
            width_x = float(width_x)
            width_y = float(width_y)
            return lambda x,y: height*exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

      def moments(self,data):
            """Returns (height, x, y, width_x, width_y)
            the gaussian parameters of a 2D distribution by calculating its
            moments """
            total = data.sum()
            X, Y = indices(data.shape)
            x = (X*data).sum()/total
            y = (Y*data).sum()/total
            col = data[:, int(y)]
            width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
            row = data[int(x), :]
            width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
            height = data.max()
            return height, x, y, width_x, width_y

      def fitgaussian(self,data):
            """Returns (height, x, y, width_x, width_y)
            the gaussian parameters of a 2D distribution found by a fit"""
            params = self.moments(data)
            errorfunction = lambda p: ravel(self.gaussian(*p)(*indices(data.shape)) - data)
            p, success = optimize.leastsq(errorfunction, params)
            return p

#==========================
      def export_data(self):
            a=h5py.File("export_data.h5")
            print np.where(self.data==inf)
            
            a.create_dataset("data",data=self.data,compression="gzip")
            a.close()

      
#==========================

      def runTransReconstruct(self):

            self.recon.numb=len(self.channelname)
            for j in arange(self.recon.numb):
                  self.recon.combo.addItem(self.channelname[j])
            self.recon.show()
            self.recon.btn.setText("Reconstruction")
            self.recon.btn.clicked.connect(self.reconstruct)
            self.recon.save.clicked.connect(self.saveRecTiff)
            self.recon.reconvalue=1


      def reconMultiply(self):
            self.rec=self.rec*10
            self.updateRecon()
      def reconDivide(self):
            self.rec=self.rec/10
            self.updateRecon()

      def updateRecon(self):
            self.reconProjNumb=self.reconView.sld.value()
            self.recon.maxText.setText(str(self.rec[self.reconProjNumb,:,:].max()))
            self.recon.minText.setText(str(self.rec[self.reconProjNumb,:,:].min()))
            
            self.reconView.view.projView.setImage(self.rec[self.reconProjNumb,:,:])




      
      def runReconstruct(self):
            self.tab_widget.removeTab(3)
            self.tab_widget.insertTab(3,self.createReconWidget(),unicode("Reconstruction"))
            self.recon.numb=len(self.channelname)
            for j in arange(self.recon.numb):
                  self.recon.combo.addItem(self.channelname[j])
            self.recon.show()
            self.recon.btn.setText("Reconstruction")
            self.recon.lcd.setText(str(self.p1[2]))
            self.recon.btn.clicked.connect(self.reconstruct)
            self.recon.save.clicked.connect(self.saveRecTiff)
            self.recon.reconvalue=0
            self.recon.mulBtn.clicked.connect(self.reconMultiply)
            self.recon.divBtn.clicked.connect(self.reconDivide)
            self.recon.maxText.setText("0")
            self.recon.minText.setText("0")
            self.reconView.sld.setRange(0,self.data.shape[2]-1)
            self.reconView.sld.valueChanged.connect(self.reconView.lcd.display)
            self.reconView.sld.valueChanged.connect(self.updateRecon)


            
      def reconstruct(self):
            self.recon.lbl.setText("Reconstruction is currently running")
            self.reconelement=self.recon.combo.currentIndex()
            self.recData=self.data[self.reconelement,:,:,:]
            self.recData[self.recData==inf]=0.00001
            self.recData[np.isnan(self.recData)]=0.00001
            #self.recData = tomopy.normalize_bg(self.recData)
            
            self.recCenter=float(self.recon.lcd.text())
            beta=float(self.recon.beta.text())
            delta=float(self.recon.delta.text())
            num_iter =int(self.recon.iters.text())
            print beta, delta, num_iter


            print "working fine"
            b=time.time()
            if self.recon.method.currentIndex()==0:
                  self.rec = tomopy.recon(self.recData, self.theta*np.pi/180, algorithm='mlem', #center=np.array(self.recCenter, dtype=float32),
                                          num_iter=num_iter, emission=True)
            elif self.recon.method.currentIndex()==1:
                  self.rec = tomopy.recon(self.recData, self.theta, algorithm='gridrec',
                                     emission=True)
            elif self.recon.method.currentIndex()==2:
                  self.rec = tomopy.recon(self.recData, self.theta*np.pi/180, algorithm='art',
                                     num_iter=num_iter, emission=True)
            elif self.recon.method.currentIndex()==3:
                  self.rec = tomopy.recon(self.recData, self.theta*np.pi/180, algorithm='pml_hybrid', #center=np.array(self.recCenter, dtype=float32),
                                     reg_par=np.array([beta,delta],dtype=np.float32), num_iter=num_iter, emission=True)
            elif self.recon.method.currentIndex()==4:
                  self.rec = tomopy.recon(self.recData, self.theta*np.pi/180, algorithm='pml_quad', center=None,#center=np.array(self.recCenter, dtype=float32),
                                     reg_par=np.array([beta,delta],dtype=np.float32), num_iter=num_iter, emission=True)

                  
            a=time.time()
            print a-b
            
            self.rec = tomopy.remove_nan(self.rec)
            self.reconProjNumb=self.projView.sld.value()
            self.reconView.view.projView.setImage(self.rec[self.reconProjNumb,:,:])
            self.reconView.setWindowTitle("Slices of reconstructed model")
            self.recon.lbl.setText("Done")
            self.recon.save.setHidden(False)

      def circular_mask(self):
            self.rec = tomopy.circ_mask(self.rec, axis=0)

      def saveRecTiff(self):
            try:
                  global debugging
                  self.savedir=QtGui.QFileDialog.getSaveFileName()
                  self.savedir=str(self.savedir)
                  if self.savedir=="":
                        raise IndexError
                  print self.savedir
                  self.circular_mask() ### temporary
                  tomopy.write_tiff_stack(self.rec,fname=self.savedir)
            except IndexError:
                  print "type the header name"
#=============================

      def reorder_matrix(self):
            argsorted= argsort(self.theta)
            print argsorted, self.theta[argsorted]
            self.data=self.data[:,argsorted,:,:]
            print "sorting done"

#=============================
                  
      def selectImageTag(self):
            self.sit = AlignWindow()
            self.sit.setWindowTitle("Seletect Image Tag from h5 file")
            self.sit.data=h5py.File(self.fileNames[0])
            self.sit.firstColumn = self.sit.data.items()
            self.sit.firstColumnNum = len(self.sit.firstColumn)
            for i in arange(self.sit.firstColumnNum):
                  self.sit.combo.addItem(self.sit.firstColumn[i][0])
            print self.sit.combo.currentIndex()

            self.sit.combo.setCurrentIndex(i)
            for i in arange(30):
                  self.sit.method.removeItem(0)
            self.sit.secondColumnName = self.sit.firstColumn[self.sit.combo.currentIndex()][0]
            
            self.sit.secondColumn = self.sit.data[self.sit.secondColumnName].items()
            self.sit.secondColumnNum = len(self.sit.secondColumn)
            for j in arange(self.sit.secondColumnNum):
                  self.sit.method.addItem(self.sit.secondColumn[j][0])


            
            self.sit.combo.currentIndexChanged.connect(self.selectImageTag_image)
            self.sit.btn.setText("Set")
            self.sit.btn2.setVisible(False)
            self.sit.btn.clicked.connect(self.setImageTag)

            self.sit.show()

      def selectImageTag_image(self):
            for i in arange(30):
                  self.sit.method.removeItem(0)
            self.sit.secondColumnName = self.sit.firstColumn[self.sit.combo.currentIndex()][0]
            
            self.sit.secondColumn = self.sit.data[self.sit.secondColumnName].items()
            self.sit.secondColumnNum = len(self.sit.secondColumn)
            for j in arange(self.sit.secondColumnNum):
                  self.sit.method.addItem(self.sit.secondColumn[j][0])

            self.sit.method.setHidden(False)

      def setImageTag(self):
            self.ImageTag = str(self.sit.combo.currentText())
            self.lbl.setText("Image Tag has been set to \""+self.ImageTag+"\"")
            print "Image Tag has been set to \"", self.ImageTag, "\""
            self.sit.setVisible(False)
#==============================


            
      def convert(self):
            self.afterConversionMenu.setDisabled(False)
            self.convert2array()
            self.x=W2(self)
            self.x.textedit.setText("All the files has been converted")
            self.x.btn.setText("OK")
            self.x.btn.clicked.connect(self.selectFilesHide)
            self.x.show()
      def selectElement(self):
            self.element=QSelect()
            f=h5py.File(os.path.abspath(self.fileNames[0]),"r")

            try:
                  self.channelnameTemp=list(f[self.ImageTag]["images_names"])
                  self.dataTag="images"

                  self.channels,self.y,self.x=f[self.ImageTag]["images"].shape
            except KeyError:
                  try:
                        self.dataTag="data"
                        self.channelnameTemp=list(f[self.ImageTag]["channel_names"])
                        self.channels,self.y,self.x=f[self.ImageTag]["data"].shape
                  except KeyError:
                        self.dataTag="XRF_roi"
                        global wha
                        self.channelnameTemp1=f[self.ImageTag]["channel_names"]

                        self.channelnameTemp2=f[self.ImageTag]["scaler_names"]
                        self.channels1,self.y,self.x=f[self.ImageTag]["XRF_roi"].shape
                        self.channels=self.channels1+len(self.channelnameTemp2)
                        self.channelnameTemp=list(self.channelnameTemp1) + list(self.channelnameTemp2)
                        
            for i in arange(len(self.channelnameTemp)):
                  self.element.button[i].setText(self.channelnameTemp[i])
                  self.element.button[i].setChecked(True)
            self.element.setWindowTitle("Select Element")
            self.element.setVisible(False)

            self.element.btn2.setText("Deselect All")
            self.element.btn.setText("set Element")
            self.element.btn.clicked.connect(self.setElement)
            self.element.btn2.clicked.connect(self.deselectAllElement)
            self.element.btn3.setVisible(False)

      def deselectAllElement(self):
            for i in arange(len(self.channelnameTemp)):
                  self.element.button[i].setText(self.channelnameTemp[i])
                  self.element.button[i].setChecked(False)
      def setElement(self):
            y=zeros(len(self.channelnameTemp),dtype=bool)
            k=arange(y.shape[0])
            for i in arange(len(self.channelnameTemp)):
                  y[i]=self.element.button[i].isChecked()
            self.channelname = [self.channelnameTemp[f] for f in k if y[f]]
            self.channelnamePos = zeros(len(self.channelname))
            for i in arange(len(self.channelname)):
                  self.channelnamePos[i] = self.channelnameTemp.index(self.channelname[i])
            self.element.setVisible(False)
      def selectElementShow(self):
            self.element.setVisible(True)

      def selectFilesHide(self):
            self.filecheck.setVisible(False)
                  

      def selectFiles(self):
            self.filecheck=QSelect()
            degree_sign= u'\N{DEGREE SIGN}'
            for i in arange(len(self.fileNames)):
                  self.fileNames[i]=str(self.fileNames[i])
                  f = h5py.File(os.path.abspath(self.fileNames[i]),"r")
                  thetatemp = f["MAPS"]["extra_pvs_as_csv"][self.thetaPos]
                  thetapos = string.rfind(thetatemp, ",")
                  theta = str(round(float(thetatemp[thetapos+1:])))
                  onlyfilename=self.fileNames[i].rfind("_")
                  if onlyfilename==-1:
                        onlyfilename=self.fileNames[i].rfind("/")
                  self.filecheck.button[i].setText(self.fileNames[i][onlyfilename+1:]+" ("+theta+degree_sign+")")
                  self.filecheck.button[i].setChecked(True)
            self.ImageTag=f.items()[-1][0]
            self.lbl.setText("Image Tag has been set to \""+self.ImageTag+"\"")
            self.filecheck.setWindowTitle("Select files")

            self.selectElement()

            self.optionMenu.setEnabled(True)
            self.filecheck.btn2.setVisible(True)
            self.filecheck.btn.clicked.connect(self.convert)
            self.filecheck.btn2.clicked.connect(self.selectImageTag)
            self.filecheck.btn3.clicked.connect(self.selectElementShow)

      def reorder_matrix(self):
            argsorted= argsort(self.theta)
            print argsorted, self.theta[argsorted]
            self.data=self.data[:,argsorted,:,:]
            print "sorting done"
            
      def selectFilesShow(self):
            self.filecheck.setVisible(True)

      def openfile(self):
            try:  
                  fileNametemp = QtGui.QFileDialog.getOpenFileNames(self, "Open File",
                                                                QtCore.QDir.currentPath(),filter="h5 (*.h5)")
                  self.fileNames=str(fileNametemp.join("\n")).split("\n")
                  if self.fileNames == [""]:
                        raise IndexError
                  RH2=self.fileNames
                  self.selectFiles()

            except IndexError:
                  print "no file has been selected"
            except IOError:
                  print "no file has been selected"


#================
      def openfolder(self):
            try:
                  folderName = QtGui.QFileDialog.getExistingDirectory(self,"Open Folder",
                                                                   QtCore.QDir.currentPath())
                  global RH,fileName
                  RH=folderName
                  folderName=str(folderName)
                  file_name_array = [ f for f in os.listdir(folderName) if isfile(join(folderName,f))]
                  file_name_array = [ f for f in os.listdir(folderName) if string.find(f,"h5")!=-1]
                  file_name_array = [folderName+"/"+f for f in file_name_array]

                  self.directory=1
                  self.fileNames=file_name_array
                  fileName=file_name_array
                  self.selectFiles()
            except IndexError:
                  print "no folder has been selected"
            except OSError:
                  print "no folder has been selected"

      def openTiffFolder(self):
            pass1= False
            try:
                  tiffFolderName = QtGui.QFileDialog.getExistingDirectory(self, "Open Tiff Folder",
                                                                          QtCore.QDir.currentPath())
                  file_name_array = [ f for f in os.listdir(tiffFolderName) if string.find(f,"tif")!=-1]
                  file_name_array = [tiffFolderName + "/" + f for f in file_name_array]

                  self.tiffNames = file_name_array

                  pass1= True
            except IndexError:
                  print "no folder has been selected"
            except OSError:
                  print "no folder has been selected"
            if pass1==True:
                  try:
                        angleTxt = QtGui.QFileDialog.getOpenFileName(self, "Open Angle Txt File",
                                                                      QtCore.QDir.currentPath(),filter="txt (*.txt)")
                        print str(angleTxt)
                        self.angleTxt =angleTxt
                        self.convertTiff2Array()
                  except IndexError:
                        print "no file has been selected"
                  except OSError:
                        print "no file has been selected"

      def convertTiff2Array(self):
            self.projections = len(self.tiffNames)
            print self.tiffNames
            self.y, self.x = np.asarray(Image.open(str(self.tiffNames[0])),dtype=float32).shape
            self.channelname = list()
            self.channelname.append("tiff")
            self.channelname.append("dummy")
            self.data =zeros([1,self.projections,self.y,self.x])
            self.theta=zeros(self.projections)
            f=open(str(self.angleTxt),'r')
            read = f.readlines()

            for i in arange(self.projections):
                  tempy,tempx = np.asarray(Image.open(str(self.tiffNames[i])),dtype=float32).shape
                  if tempy>self.y:
                          self.y=tempy
                  if tempx>self.x:
                          self.x=tempx
            self.data=zeros([1,self.projections,self.y,self.x])
            for i in arange(self.projections):
                  tempy,tempx=np.asarray(Image.open(str(self.tiffNames[i])),dtype=float32).shape
                  self.data[0,i,:tempy,:tempx]=np.asarray(Image.open(str(self.tiffNames[i])),dtype=float32)[...]
                  thetaPos=read[i].find(",")
                  self.theta[i] = float(read[i][thetaPos+1:])
            self.p1=[100,100,self.data.shape[3]/2]

            self.xshift=zeros(self.projections,int)
            self.yshift=zeros(self.projections,int)
            self.alignmentMenu.setEnabled(True)
            self.tab_widget.setEnabled(True)
            self.afterConversionMenu.setDisabled(False)
            
            self.oldData=zeros(self.data.shape)
            self.oldData[...]=self.data[...]

            self.showImageProcess()
            self.showSaveHotSpotPos()
            self.projView.view.hotSpotSetNumb=0
            self.showSinogram()
            self.sinogram()
            self.runReconstruct()

            self.selectedFiles=list()
            for i in arange(self.projections):
                  self.selectedFiles.append(str(self.tiffNames[i]))

            self.projView.view.hotSpotSetNumb=0
            self.tab_widget.setCurrentIndex(0)

            self.updateImages()           
      def convert2array(self):
            self.setElement()
            y=zeros(len(self.fileNames),dtype=bool)
            k=arange(y.shape[0])
            for i in arange(len(self.fileNames)):
                  y[i]=self.filecheck.button[i].isChecked()
      
            self.selectedFiles=[self.fileNames[f] for f in k if y[f]==True]

            ### From first data retrieve channel names, size of the image.            
            f=h5py.File(os.path.abspath(self.selectedFiles[0]),"r")
            try:
                  
                  self.dataTag="images"

                  self.channels,self.y,self.x=f[self.ImageTag]["images"].shape
            except KeyError:
                  try:
                        self.dataTag="data"
                        
                        self.channels,self.y,self.x=f[self.ImageTag]["data"].shape
                  except KeyError:
                        self.dataTag="XRF_roi"
                        self.channels1,self.y,self.x=f[self.ImageTag]["XRF_roi"].shape
                        self.channels=self.channels1+len(self.channelnameTemp2)
                        
            self.projections=len(self.selectedFiles)
            self.theta= zeros(self.projections)
            

            self.xshift=zeros(self.projections,int)
            self.yshift=zeros(self.projections,int)

            
            self.channels=len(self.channelname)

            #### check the size of the images
            x=0
            y=0
            for i in arange(self.projections):
                        file_name = os.path.abspath(self.selectedFiles[i])
                        f = h5py.File(file_name,"r")
                        dummy, tempY, tempX =f[self.ImageTag][self.dataTag][...].shape
                        if tempX > x:
                              x = tempX
                        if tempY > y:
                              y = tempY
            self.x=x
            self.y=y
                        
            if self.dataTag!="XRF_roi":
                  self.data=zeros([self.channels,self.projections,self.y,self.x])
                  for i in arange(self.projections):
                        file_name = os.path.abspath(self.selectedFiles[i])
                        f = h5py.File(file_name,"r")
                        thetatemp=f["MAPS"]["extra_pvs_as_csv"][self.thetaPos]

                        self.theta[i] = float(thetatemp[thetatemp.rfind(",")+1:])

      
                        for j in arange(len(self.channelnamePos)):
                              pos=self.channelnamePos[j]
                              imgY,imgX = f[self.ImageTag][self.dataTag][0,:,:].shape
                              self.data[j,i,:imgY,:imgX]=f[self.ImageTag][self.dataTag][pos,:,:]
                        print i+1, "projection(s) has/have been converted"
                  print "worked"

            if self.dataTag=="XRF_roi":
                  self.data=zeros([self.channels,self.projections,self.y,self.x])
                  for i in arange(self.projections):
                        file_name = os.path.abspath(self.selectedFiles[i])
                        f = h5py.File(file_name,"r")
                        thetatemp=f["MAPS"]["extra_pvs_as_csv"][self.thetaPos]

                        self.theta[i] = float(thetatemp[thetatemp.rfind(",")+1:])

                        for j in arange(len(self.channelnamePos)):
                              
                              if self.channelnamePos[j]<len(list(self.channelnameTemp1)):
                                    pos=self.channelnamePos[j]
                                    imgY,imgX = f[self.ImageTag][self.dataTag][0,:,:].shape
                                    self.data[j,i,:imgY,:imgX]=f[self.ImageTag][self.dataTag][pos,:,:]
                              else:
                                    pos = self.channelnamePos[j]-len(list(self.channelnameTemp1))
                                    imgY,imgX = f[self.ImageTag][self.dataTag][0,:,:].shape
                                    self.data[j,i,:imgY,:imgX]=f[self.ImageTag][self.dataTag][pos,:,:]

                        print i+1, "projection(s) has/have been converted"
                  print "worked"

            

            self.data[isnan(self.data)]=0.0001
            self.data[self.data==inf]=0.0001
            self.p1=[100,100,self.data.shape[3]/2]

            print self.p1



            self.alignmentMenu.setEnabled(True)
            self.tab_widget.setEnabled(True)
            
            self.oldData=zeros(self.data.shape)
            self.oldData[...]=self.data[...]

            global datadata
            datadata=self.oldData

            self.showImageProcess()
            self.showSaveHotSpotPos()
            self.projView.view.hotSpotSetNumb=0
            self.showSinogram()
            self.sinogram()

            self.runReconstruct()
            self.tab_widget.setCurrentIndex(0)

            self.updateImages()

      def saveThetaTxt(self):
            try:
                  self.alignFileName = QtGui.QFileDialog.getSaveFileName()
                  if string.rfind(str(self.alignFileName),".txt")==-1:
                        self.alignFileName=str(self.alignFileName)+".txt"
                  print str(self.alignFileName)
                  f=open(self.alignFileName,"w")
                  for i in arange(self.projections):
                        onlyfilename=self.selectedFiles[i].rfind("/")
                        print self.selectedFiles[i]
                        f.writelines(self.selectedFiles[i][onlyfilename+1:-3]+", "+str(self.theta[i])+"\n")
                        
                  f.close()
            except IOError:
                  print "choose file please"
            

#####!!!! just temp need to be fixed
      def saveImage(self):
            self.saveImageDir=QtGui.QFileDialog.getExistingDirectory()
            print self.saveImageDir
            for j in arange(self.data.shape[0]):
                  path=str(self.saveImageDir)+"/"+self.channelname[j]
                  try:
                        os.makedirs(path)
                  except OSError:
                        path=path
                  for i in arange(self.data.shape[1]):
                        temp_img=self.data[j,i,:,:]
                        temp=Image.fromarray(temp_img.astype(np.float32))

                        index=string.rfind(self.selectedFiles[i],"/")
                        temp.save(path+self.selectedFiles[i][index:-3]+"_"+self.channelname[j]+".tif")

      def multiplier(self):
            self.data =10*self.data
            self.updateImages()

      def divider(self):
            self.data = self.data*0.1
            self.updateImages()
            
      def normalize(self):
            self.normElement = self.projection.combo.currentIndex()
            normData = self.data[self.normElement,:,:,:]
            for i in arange(normData.shape[0]):
                  print normData.dtype
                  temp=normData[i,:,:][...]
                  tempMax=temp.max()
                  tempMin=temp.min()
                  temp=(temp-tempMin)/tempMax
                  self.data[self.normElement,i,:,:]=temp


      def showSinogram(self):
            #self.sino = QSelect2()
            #self.sino.setWindowTitle("Sinogram Window")
            self.tab_widget.removeTab(2)
            self.tab_widget.insertTab(2,self.createSinoWidget(),unicode("Sinogram"))

            self.sino.numb=len(self.channelname)
            for j in arange(self.sino.numb):
                  self.sino.combo.addItem(self.channelname[j])
            self.sino.show()
            self.sino.btn.clicked.connect(self.runCenterOfMass2)
            self.sino.btn.setText("center of mass")
            self.sino.btn2.clicked.connect(self.sinoShift)
            self.sino.sld.setRange(5,self.y-5)
            self.sino.lcd.display(5)
            self.sino.sld.valueChanged.connect(self.sino.lcd.display)
            self.sino.sld.valueChanged.connect(self.sinogram)
            
            self.sinoView.show()

      def sinogram(self):
            self.thickness=10
            self.sinoelement=self.sino.combo.currentIndex()
            sinodata=self.data[self.sinoelement,:,:,:]
            self.sinogramData=zeros([sinodata.shape[0]*self.thickness,sinodata.shape[2]],dtype=float32)

            for i in arange(self.projections):
                  self.sinogramData[i*self.thickness:(i+1)*self.thickness,:]=sinodata[i,self.sino.sld.value()-self.thickness/2:self.sino.sld.value()+self.thickness/2,:]
            global sinofig

            sinofig=self.sinogramData
            self.sinogramData[isinf(self.sinogramData)]=0.001
            self.sinoView.view.projView.setImage(self.sinogramData)
            self.sinoView.view.projView.setRect(QtCore.QRect(0,0,self.sinogramData.shape[0],self.sinogramData.shape[1]))
            self.sinoView.view.projData=self.sinogramData
            self.sinoView.view.getShape()
   
            
            #self.sinoView.getShape()
            #self.sinoView.setWindowTitle("Sinogram "+self.sino.combo.itemText(self.sinoelement)+" "+str(self.sino.sld.value()))

      def saveSinogram(self):
            j=Image.fromarray(self.sinogramData.astype(np.float32))
            j.save("sinogram.tif")

      def sinoShift(self):
            for i in arange(self.projections):
                  
                  self.data[:,i,:,:]=np.roll(self.data[:,i,:,:], self.sinoView.view.regShift[i], axis=2)

      def updateImages(self):
            self.projView.view.projView.update()
            self.imgProcess.view.projView.updateImage()
            self.sinoView.view.projView.updateImage()

      def openPositionPanel(self):
            self.positionPanel = PositionPanel()
            self.positionPanel.show()
            self.positionPanel.lbl1 = str(self.projView.view.projView.iniY)


            
            

####==============
      def readConfigFile(self):
            try:
                  fileName = QtGui.QFileDialog.getOpenFileName(self, "Select config file",
                        QtCore.QDir.currentPath(),"TXT (*.txt)")
                  ##### for future reference "All File (*);;CSV (*.csv *.CSV)"
                  print fileName
                  f=open(fileName,'r')
                  l = f.readlines()
                  for i in arange(len(l)):
                          ## image_tag
                        if string.find(l[i],"Image_Tag")!=-1:
                              self.ImageTag=l[i+1][:-1]
                          ## theta postion from extra pvs as csv
                        if string.find(l[i],"Theta position")!=-1:
                              print l[i+1][:-1]
                              self.thetaPos=int(l[i+1][:-1])
                  self.selectFiles()                                   
            except IOError:
                  print "choose file please"



            
####==============
class W2(QtGui.QDialog):
      def __init__(self, parent=None):
            super(W2, self).__init__(parent)
            self.textedit=QtGui.QLabel()
            self.textedit.setText("What")
            self.btn = QtGui.QPushButton('Click2')
            self.btn.setText("OK")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.textedit)
            vb.addWidget(self.btn)
            
            self.setLayout(vb)
            self.btn.clicked.connect(self.fireupWindows3)

      def fireupWindows3(self):
            
            self.accept()

####===============

class AlignWindow(QtGui.QWidget):
      def __init__(self):
            super(AlignWindow,self).__init__()
            self.initUI()

      def initUI(self):
            self.combo = QtGui.QComboBox(self)
            self.method = QtGui.QComboBox(self)

            self.btn = QtGui.QPushButton('Click2')
            self.btn2 = QtGui.QPushButton("Click3")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
            vb.addWidget(self.method)
            vb.addWidget(self.btn)
            vb.addWidget(self.btn2)
            self.setLayout(vb)
           

class QSelect(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect, self).__init__()
        
            self.initUI()
        
      def initUI(self):
            names=list()
            for i in arange(100):
                  names.append("")


            self.grid = QtGui.QGridLayout()
            self.lbl=QtGui.QLabel()
            self.lbl2=QtGui.QLabel()
            self.lbl.setText("closing this window won't affect your selection of the files")
            self.lbl2.setText("You should convert the files in order to generate sinogram or reconstructed data")
            self.btn=QtGui.QPushButton()
            self.btn.setText("Save Data in Memory")
            self.btn2=QtGui.QPushButton()
            self.btn2.setText("set Image Tag")
            self.btn3=QtGui.QPushButton("set Element")

            j = 0
            pos=list()
            for y in arange(10):
                  for x in arange(10):
                        pos.append((x,y))

            self.button=list()
            for i in names:
                  self.button.append(QtGui.QCheckBox(i))
                  self.grid.addWidget(self.button[j], pos[j][0], pos[j][1])
                  j = j + 1
            self.setLayout(self.grid)

            self.vb=QtGui.QVBoxLayout()
            self.vb2=QtGui.QVBoxLayout()
            
            self.vb.addWidget(self.lbl,11)
            self.vb.addWidget(self.lbl2,12)

            hb = QtGui.QHBoxLayout()
            hb.addWidget(self.btn2)
            hb.addWidget(self.btn3)
            self.vb2.addLayout(hb)
            self.vb2.addWidget(self.btn)

            self.grid.addLayout(self.vb,11,0,1,10)
            self.grid.addLayout(self.vb2,13,3,1,2)
        
            self.move(100, 100)
            self.setWindowTitle('Calculator')    
            self.show()
####==================
class QSelect2(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect2, self).__init__()
        
            self.initUI()
        
      def initUI(self):
 
            self.sld=QtGui.QSlider(QtCore.Qt.Horizontal, self)
            self.lcd=QtGui.QLCDNumber(self)
            self.combo = QtGui.QComboBox(self)
            self.btn = QtGui.QPushButton('Click2')
            self.btn.setText("Sinogram")
            self.btn2 = QtGui.QPushButton("shift data")
            self.btn3 = QtGui.QPushButton("X 10")
            self.btn4 = QtGui.QPushButton("/ 10")
            hb=QtGui.QHBoxLayout()
            hb.addWidget(self.btn3)
            hb.addWidget(self.btn4)
            self.btn3.setVisible(False)
            self.btn4.setVisible(False)
            self.lbl=QtGui.QLabel()
            self.lbl.setText("")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
            vb.addWidget(self.btn)
            vb.addWidget(self.btn2)
            vb.addWidget(self.lcd)
            vb.addWidget(self.sld)
            vb.addWidget(self.lbl)
            vb.addLayout(hb)
            self.setLayout(vb)
####===================
class QSelect3(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect3, self).__init__()
        
            self.initUI()
        
      def initUI(self):
 
                  
            self.combo = QtGui.QComboBox(self)
            self.method = QtGui.QComboBox(self)
            self.btn = QtGui.QPushButton('Click2')
            self.save = QtGui.QPushButton("Save tiff files")
            self.save.setHidden(True)
            self.btn.setText("Sinogram")
            self.sld=QtGui.QSlider(QtCore.Qt.Horizontal, self)
            self.lcd=QtGui.QLineEdit("0")
            self.lbl=QtGui.QLabel()
            self.lbl2=QtGui.QLabel("Center")
            self.lbl.setText("")

            centerBox = QtGui.QHBoxLayout()
            centerBox.addWidget(self.lbl2)
            centerBox.addWidget(self.lcd)
            self.methodname=["mlem", "gridrec", "art","pml_hybrid","pml_quad"]

            self.mulBtn = QtGui.QPushButton("x 10")
            self.divBtn = QtGui.QPushButton("/ 10")
            mdBox=QtGui.QHBoxLayout()
            mdBox.addWidget(self.mulBtn)
            mdBox.addWidget(self.divBtn)

            self.maxLbl = QtGui.QLabel("Max")
            self.minLbl = QtGui.QLabel("Min")
            self.maxText = QtGui.QLineEdit()
            self.minText = QtGui.QLineEdit()

            maxBox = QtGui.QHBoxLayout()
            minBox = QtGui.QHBoxLayout()
            maxBox.addWidget(self.maxLbl)
            maxBox.addWidget(self.maxText)
            minBox.addWidget(self.minLbl)
            minBox.addWidget(self.minText)
            
            self.betaName=QtGui.QLabel("Beta")
            self.deltaName = QtGui.QLabel("Delta")
            self.itersName=QtGui.QLabel("Iteration")
            self.beta = QtGui.QLineEdit("1")
            self.delta=QtGui.QLineEdit("0.01")
            self.iters = QtGui.QLineEdit("10")

            betaBox=QtGui.QHBoxLayout()
            deltaBox=QtGui.QHBoxLayout()
            itersBox= QtGui.QHBoxLayout()
            betaBox.addWidget(self.betaName)
            betaBox.addWidget(self.beta)
            deltaBox.addWidget(self.deltaName)
            deltaBox.addWidget(self.delta)
            itersBox.addWidget(self.itersName)
            itersBox.addWidget(self.iters)


            
            for k in arange(len(self.methodname)):
                  self.method.addItem(self.methodname[k])
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
            vb.addWidget(self.method)
            vb.addWidget(self.btn)
            vb.addWidget(self.save)
            vb.addLayout(centerBox)
            vb.addWidget(self.sld)
            vb.addWidget(self.lbl)
            vb.addLayout(mdBox)
            vb.addLayout(maxBox)
            vb.addLayout(minBox)
            vb.addLayout(itersBox)
            vb.addLayout(betaBox)
            vb.addLayout(deltaBox)
            self.setLayout(vb)

class QSelect4(QtGui.QWidget):
    
      def __init__(self):
            super(QSelect4, self).__init__()
        
            self.initUI()
        
      def initUI(self):
 


            self.sld=QtGui.QSlider(QtCore.Qt.Horizontal, self)
            self.lcd=QtGui.QLCDNumber(self)
            self.combo = QtGui.QComboBox(self)
            self.combo2=QtGui.QComboBox(self)
            self.combo3=QtGui.QComboBox(self)
            self.lbl1 = QtGui.QLabel("Set the size of the hotspot")
            self.lbl3 = QtGui.QLabel()
            self.lbl3.setText("Set a group number of the hot spot")
            for i in arange(5):
                  self.combo2.addItem(str(i+1))
            self.btn=QtGui.QPushButton("Hotspots to a line")
            self.btn2=QtGui.QPushButton("Hotspots to a sine curve")
            self.btn3=QtGui.QPushButton("set y")
            self.btn4=QtGui.QPushButton("Clear hotspot data")
##            self.btn = QtGui.QPushButton('Click2')
##            self.btn.setText("Sinogram")
##            self.btn2 = QtGui.QPushButton("shift data")
##            self.lbl=QtGui.QLabel()
##            self.lbl.setText("")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.combo)
##            vb.addWidget(self.btn)
##            vb.addWidget(self.btn2)

            vb.addWidget(self.lbl1)
            vb.addWidget(self.lcd)
            vb.addWidget(self.sld)
            vb.addWidget(self.combo3)

            hb1=QtGui.QVBoxLayout()
            hb1.addWidget(self.lbl3,0)
            hb1.addWidget(self.combo2)

            vb.addLayout(hb1)
            vb.addWidget(self.btn)
            vb.addWidget(self.btn2)
            vb.addWidget(self.btn3)
            vb.addWidget(self.btn4)
##            vb.addWidget(self.lbl)
            self.setLayout(vb)

class imageProcess(QtGui.QWidget):
      def __init__(self):
            super(imageProcess, self).__init__()
        
            self.initUI()

      def initUI(self):
            self.xSize = 20
            self.ySize = 20

##            self.mulBtn=QtGui.QPushButton("x10")
##            self.divBtn=QtGui.QPushButton("/10")
            self.bgBtn=QtGui.QPushButton("Bg Value")
            self.delHotspotBtn = QtGui.QPushButton("Delete HS")
            self.normalizeBtn = QtGui.QPushButton("Normalize")
            self.cutBtn = QtGui.QPushButton("Cut")
            self.gaussian33Btn = QtGui.QPushButton("3*3 gauss")
            self.gaussian55Btn = QtGui.QPushButton("5*5 gauss")
            self.xUpBtn=QtGui.QPushButton("x: +")
            self.xUpBtn.clicked.connect(self.xUp)
            self.xDownBtn = QtGui.QPushButton("x: -")
            self.xDownBtn.clicked.connect(self.xDown)
            self.yUpBtn = QtGui.QPushButton("y: +")
            self.yUpBtn.clicked.connect(self.yUp)
            self.yDownBtn = QtGui.QPushButton("y: -")
            self.yDownBtn.clicked.connect(self.yDown)

            self.xSizeLbl = QtGui.QLabel("x Size")
            self.ySizeLbl = QtGui.QLabel("y Size")

            self.xSizeTxt = QtGui.QLineEdit(str(self.xSize))
            self.ySizeTxt = QtGui.QLineEdit(str(self.ySize))



            self.combo1 = QtGui.QComboBox()
            self.combo2 = QtGui.QComboBox()
            
            hb1=QtGui.QHBoxLayout()
            hb1.addWidget(self.xUpBtn)
            hb1.addWidget(self.xDownBtn)
            hb2=QtGui.QHBoxLayout()
            hb2.addWidget(self.xSizeLbl)
            hb2.addWidget(self.xSizeTxt)
            hb3=QtGui.QHBoxLayout()
            hb3.addWidget(self.yUpBtn)
            hb3.addWidget(self.yDownBtn)
            hb4=QtGui.QHBoxLayout()
            hb4.addWidget(self.ySizeLbl)
            hb4.addWidget(self.ySizeTxt)

            vb1=QtGui.QVBoxLayout()
            vb1.addLayout(hb1)
            vb1.addLayout(hb2)
            vb2=QtGui.QVBoxLayout()
            vb2.addLayout(hb3)
            vb2.addLayout(hb4)
            xSG= QtGui.QGroupBox("x Size")
            xSG.setLayout(vb1)
            ySG = QtGui.QGroupBox("y Size")
            ySG.setLayout(vb2)
            vb3=QtGui.QVBoxLayout()
            vb3.addWidget(self.combo1)
            vb3.addWidget(self.combo2)
            vb3.addWidget(xSG)
            vb3.addWidget(ySG)



            hb6=QtGui.QHBoxLayout()
            hb6.addWidget(self.bgBtn,stretch =0 )
            hb6.addWidget(self.delHotspotBtn,stretch =0)
            
            hb7=QtGui.QHBoxLayout()
            hb7.addWidget(self.normalizeBtn,stretch =0)
            hb7.addWidget(self.cutBtn,stretch =0)
            
            hb8=QtGui.QHBoxLayout()
            hb8.addWidget(self.gaussian33Btn,stretch =0)
            hb8.addWidget(self.gaussian55Btn,stretch =0)
            
            vb3.addLayout(hb6)
            vb3.addLayout(hb7)
            vb3.addLayout(hb8)
 

            self.setLayout(vb3)
      def changeXSize(self):
            self.xSize=int(self.xSizeTxt.text())
      def changeYSize(self):
            self.ySize=int(self.ySizeTxt.text())
      def xUp(self):
            self.changeXSize()
            self.changeYSize()
            self.xSize+=1
            self.xSizeTxt.setText(str(self.xSize))
      def xDown(self):
            self.changeXSize()
            self.changeYSize()
            self.xSize-=1
            self.xSizeTxt.setText(str(self.xSize))
      def yUp(self):
            self.changeXSize()
            self.changeYSize()
            self.ySize+=1
            self.ySizeTxt.setText(str(self.ySize))
      def yDown(self):
            self.changeXSize()
            self.changeYSize()
            self.ySize-=1
            self.ySizeTxt.setText(str(self.ySize))

class PositionPanel(QtGui.QWidget):
      def __init__(self):
            super(PositionPanel, self).__init__()

            self.initUI()
  
      def initUI(self):
            self.lbl3=QtGui.QLabel("X position")
            self.lbl4=QtGui.QLabel("Y position")
            self.lbl1=QtGui.QLabel()
            self.lbl1.setText("")
            self.lbl2=QtGui.QLabel()
            self.lbl2.setText("")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.lbl3)
            vb.addWidget(self.lbl1)
            vb.addWidget(self.lbl4)
            vb.addWidget(self.lbl2)
            self.setLayout(vb)
            self.setGeometry(0,650,1000,150)

class IView(pg.GraphicsLayoutWidget):
    
      def __init__(self):
            super(IView, self).__init__()

            self.initUI()
            self.hotSpotNumb=0
        
      def initUI(self):
            self.show()
            self.p1=self.addPlot()
            
            self.projView = pg.ImageItem()
            self.projView.rotate(-90)
            self.p1.addItem(self.projView)

      def keyPressEvent(self, ev):

                  
            if ev.key() == QtCore.Qt.Key_Right:

                  self.getMousePos()
                  self.shiftnumb=1
                  self.shift()
                  self.projView.setImage(self.copy)
                  self.regShift[self.numb2]+=self.shiftnumb

            if ev.key() == QtCore.Qt.Key_Left:

                  self.getMousePos()
                  self.shiftnumb=-1
                  self.shift()
                  self.projView.setImage(self.copy)
                  self.regShift[self.numb2]+=self.shiftnumb

##            if ev.key() == QtCore.Qt.Key_N:
##                  print "n"
##                  print self.data.shape[0], self.hotSpotNumb
##                  print self.imageItem.iniX, self.imageItem.iniY
##                  self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,0]=self.imageItem.iniY
##                  self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,1]=self.imageItem.iniX
##                  if self.hotSpotNumb<self.data.shape[0]-1:
##                        self.hotSpotNumb+=1
##                        self.setImage(self.data[self.hotSpotNumb,:,:])
##
##            if ev.key() ==QtCore.Qt.Key_S:
##                  
##                  self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,0]=0
##                  self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,1]=0
##                  if self.hotSpotNumb<self.data.shape[0]-1:
##                        self.hotSpotNumb+=1
##                        self.setImage(self.data[self.hotSpotNumb,:,:])

                  
      def getMousePos(self):
            numb=self.projView.iniY
            self.numb2=int(numb/10)


      def shift(self):
            self.copy=self.projData
            self.copy[self.numb2*10:self.numb2*10+10,:]=np.roll(self.copy[self.numb2*10:self.numb2*10+10,:],self.shiftnumb,axis=1)
      def getShape(self):
            self.regShift=zeros(self.projData.shape[0],dtype=int)

class IView4(pg.QtGui.QWidget):
      def __init__(self):
            super(IView4, self).__init__()
        
            self.initUI()

      def initUI(self):
            self.show()

            hb1=QtGui.QHBoxLayout()

            self.view=IView()
            self.hist=pg.HistogramLUTWidget()
            self.hist.setImageItem(self.view.projView)


            hb1.addWidget(self.view)
            hb1.addWidget(self.hist,10)
            self.setLayout(hb1)


class IView2(pg.GraphicsLayoutWidget):
    
      def __init__(self):
            super(IView2, self).__init__()
        
            self.initUI()
            self.boxSize=20
            self.hotSpotNumb=0
            self.xSize=20
            self.ySize=20
        
      def initUI(self):
            self.show()

            self.p1=self.addPlot()
            
            self.projView=pg.ImageItem()
            self.projView.rotate(-90)
            self.projView.iniX=0
            self.projView.iniY=0
            self.ROI = pg.ROI([self.projView.iniX,self.projView.iniY],[20,20])
            self.p1.addItem(self.projView)
            self.p1.addItem(self.ROI)

##            self.hist = pg.HistogramLUTItem()
##            self.hist.setImageItem(self.projView)
##            self.addItem(self.hist)
      
      def mouseReleaseEvent(self,ev):

            self.ROI.setPos([self.projView.iniX-self.xSize/2,-self.projView.iniY-self.ySize/2])

      def keyPressEvent(self, ev):
##            if ev.key() == QtCore.Qt.Key_M:
##                  self.ROI.setPos([self.projView.iniX-10,self.projView.iniY-10])
            if ev.key() == QtCore.Qt.Key_N:
                  if self.hotSpotNumb<self.data.shape[0]:
                        print "n"
                        print "Total projections",self.data.shape[0], "current position",self.hotSpotNumb+1,"group number",  self.hotSpotSetNumb+1
                        self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,0]=self.projView.iniY
                        self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,1]=self.projView.iniX
                        print self.projView.iniX, self.projView.iniY

                        
                        self.hotSpotNumb+=1
                        if self.hotSpotNumb < self.data.shape[0]:
                              self.projView.setImage(self.data[self.hotSpotNumb,:,:])
                        else:
                              print "This is the last projection"


            if ev.key() ==QtCore.Qt.Key_S:
                  
                  self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,0]=0
                  self.posMat[self.hotSpotSetNumb,self.hotSpotNumb,1]=0
                  if self.hotSpotNumb<self.data.shape[0]-1:
                        self.hotSpotNumb+=1
                        self.projView.setImage(self.data[self.hotSpotNumb,:,:])


                  
##      def getMousePos(self):
##            numb=self.imageItem.iniX
##            self.numb2=int(numb/10)
##
##
##      def shift(self):
##            self.copy=self.getProcessedImage()
##            self.copy[self.numb2*10:self.numb2*10+10,:]=np.roll(self.copy[self.numb2*10:self.numb2*10+10,:],self.shiftnumb,axis=1)
##      def getShape(self):
##            self.regShift=zeros(self.getProcessedImage().shape[0],dtype=int)

class IView3(QtGui.QWidget):
      def __init__(self):
            super(IView3, self).__init__()
        
            self.initUI()

      def initUI(self):
            self.show()

            hb3=QtGui.QHBoxLayout()
            lbl1=QtGui.QLabel("x pos")
            self.lbl2=QtGui.QLabel("")
            lbl3=QtGui.QLabel("y pos")
            self.lbl4=QtGui.QLabel("")
            btn1=QtGui.QPushButton("position")
            hb3.addWidget(lbl1)
            hb3.addWidget(self.lbl2)
            hb3.addWidget(lbl3)
            hb3.addWidget(self.lbl4)
            hb3.addWidget(btn1)

            btn1.clicked.connect(self.updatePanel)
            
            hb2=QtGui.QHBoxLayout()
            hb1=QtGui.QHBoxLayout()
            vb1=QtGui.QVBoxLayout()
            self.view=IView2()
            self.sld=QtGui.QSlider(QtCore.Qt.Horizontal, self)
            self.lcd=QtGui.QLCDNumber(self)
            self.hist=pg.HistogramLUTWidget()
            self.hist.setImageItem(self.view.projView)

            hb2.addWidget(self.lcd)
            hb2.addWidget(self.sld)
            vb1.addLayout(hb3)
            vb1.addWidget(self.view)
            vb1.addLayout(hb2)
            hb1.addLayout(vb1)
            hb1.addWidget(self.hist,10)
            self.setLayout(hb1)

      def updatePanel(self):
            
            self.lbl2.setText(str(self.view.projView.iniY))
            self.lbl4.setText(str(self.view.projView.iniX))

            
class Manual(QtGui.QWidget):
      def __init__(self):
            super(Manual, self).__init__()

            self.initUI()
  
      def initUI(self):
            self.lbl1=QtGui.QLabel()
            self.lbl1.setText("Press N after clicking to the point where you are interested in")
            self.lbl2=QtGui.QLabel()
            self.lbl2.setText("Press S to skip to the next projection ")
            vb = QtGui.QVBoxLayout()
            vb.addWidget(self.lbl1)
            vb.addWidget(self.lbl2)
            self.setLayout(vb)
            self.setGeometry(0,650,1000,150)
            self.setFixedSize(1000,150)



#########################            
def main():
    
      app = QtGui.QApplication(sys.argv)
      ex = Example()
      sys.exit(app.exec_())

if __name__ == '__main__':
    main()
