# -*- coding: utf-8 -*-
#!/usr/bin/python


import sys
from sys import platform
import tkFileDialog
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

def openfile():
      try:  
            fileNametemp = QtGui.QFileDialog.getOpenFileNames(caption="Open File",
                                                                directory=QtCore.QDir.currentPath(),filter="h5 (*.h5)")
            fileNames=str(fileNametemp.join("\n")).split("\n")
            if fileNames == [""]:
                  raise IndexError
            RH2=fileNames
            #self.selectFiles()

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
