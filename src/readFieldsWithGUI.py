#! /usr/bin/env python

# *******************************************************************************
# This file is part of tools_piccante.
# 
# tools_piccante is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# tools_piccante is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with tools_piccante.  If not, see <http://www.gnu.org/licenses/>.
# 
# *******************************************************************************


# -*- coding: utf-8 -*-
#
import os, os.path, glob, sys, shutil, time
from PyQt4.QtCore import pyqtSlot
from PyQt4 import QtGui
import struct
#from scipy import *
import numpy as np



class Window(QtGui.QMainWindow):
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(100,100,600,500)
        self.setWindowTitle("Read Fields")
        
        self.addMyMenu()
        self.comp = 0
        self.textComp = 'x'
        self.outpath = os.path.join(os.getcwd(),"output.dat")
        self.inputpath     = os.getcwd()
        self.home()
        
    
    def home(self):
        #self.addMyButtons()
        self.addMyToolbar()
        
        
        self.Flag = False
        self.addLabels()
        self.addCompChoice()
        self.show()
    
    def addLabels(self):
        xpos = 10
        ypos = 200
        vsize = 20
        self.statusDisplay = QtGui.QLabel("Status: not ready!", self)
        self.statusDisplay.move(xpos,ypos)
        self.statusDisplay.resize(480,vsize)
        ypos+=vsize
        self.fileDisplay = QtGui.QLabel("No File Chosen", self)
        self.fileDisplay.move(xpos,ypos)
        self.fileDisplay.resize(480,vsize)
        ypos+=vsize
        self.fileInfo1 = QtGui.QLabel("File info:", self)
        self.fileInfo1.move(xpos,ypos)
        self.fileInfo1.resize(480,vsize)
        ypos+=vsize
        text = "Grid:"
        self.fileInfo2 = QtGui.QLabel(text, self)
        self.fileInfo2.move(xpos,ypos)
        self.fileInfo2.resize(480,vsize)
        ypos+=vsize
        self.fileInfo3 = QtGui.QLabel("", self)
        self.fileInfo3.move(xpos,ypos)
        self.fileInfo3.resize(480,vsize)
        
        self.progress = QtGui.QProgressBar(self)
        self.progress.setGeometry(10,80,250,20)
        
    def addMyToolbar(self):
        openAction = QtGui.QAction('Open File', self)
        openAction.triggered.connect(self.openFile)
        saveAction = QtGui.QAction('Write File', self)
        saveAction.triggered.connect(self.saveFile)
        
        self.toolbar = self.addToolBar("openAction")
        self.toolbar.addAction(openAction)
        self.toolbar.addAction(saveAction)
          
    def addMyMenu(self):
        openAction = QtGui.QAction("&Open", self)
        openAction.setShortcut("Ctrl+O")
        openAction.setStatusTip('Leave the App')
        openAction.triggered.connect(self.openFile)
        
        saveAction = QtGui.QAction("&Save", self)
        saveAction.setShortcut("Ctrl+S")
        saveAction.setStatusTip('output file')
        saveAction.triggered.connect(self.saveFile)
        
        quitAction = QtGui.QAction("&Quit", self)
        quitAction.setShortcut("Ctrl+Q")
        quitAction.setStatusTip('quit')
        quitAction.triggered.connect(self.close_app)
        
        mainMenu = self.menuBar()
        #mainMenu.setNativeMenuBar(False)
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(openAction)
        fileMenu.addAction(saveAction)
        fileMenu.addAction(quitAction)
        
        
        self.statusBar()
    
    def addCompChoice(self):
        xpos = 10
        ypos = 150
        xsize = 150 
        self.comp = 0
        self.componentLabel = QtGui.QLabel("Components = ?", self)
        self.componentLabel.move(xpos,ypos)
        self.componentLabel.resize(150,20)
        
        self.compChoice = QtGui.QComboBox(self)
        self.compChoice.move(xpos+xsize,ypos)
        self.compChoice.activated[str].connect(self.comp_choice)
        
    def addFileInfo(self):
        text = "Total points: %d = [%d : %d : %d] [Nx : Ny : Nz]" %(self.Ntot,self.Nx, self.Ny, self.Nz)
        self.fileInfo2.setText(text)
        
        
    def changeCompChoice(self):
        self.compChoice.clear()
        if self.Nc >=1:
            text = 'Components = %d' %self.Nc
            self.componentLabel.setText(text)
            self.compChoice.addItem("x")
        if self.Nc >=2 :
            self.compChoice.addItem("y")
        if self.Nc >=3 :
            self.compChoice.addItem("z")
        
    def addMyButtons(self):
        btn = QtGui.QPushButton("quit", self)
        btn.clicked.connect(self.close_app)
        btn.resize(70,20)
        btn.move(0,100)
        
        
        btn2 = QtGui.QPushButton("openFile", self)
        btn2.clicked.connect(self.openFile)
        btn2.resize(70,20)
        btn2.move(70,100)
        
        btn3 = QtGui.QPushButton("WRITE", self)
        btn3.clicked.connect(self.saveFile)
        btn3.resize(70,20)
        btn3.move(140,100)
            
    def comp_choice(self, text):
        label = "Chosen component: " + text
        self.componentLabel.setText(label)
        self.comp = 0
        self.textComp = text
        if text == "y":
            self.comp = 1
        elif text == "z":
            self.comp = 2
        else:
            pass
        print "changed component to " + text
        
    def close_app(self):
        mymessage = "Sei sicuro?"
        result = QtGui.QMessageBox.question(self, 'Message', mymessage, QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
 
        if result == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass
    
    def myDefaultOutput(self):
        if self.Nc == 1:
            self.outpath = os.path.join(os.getcwd(),"output.dat")
        else:
            self.outpath = os.path.join(os.getcwd(),"output_%s.dat"%self.textComp)
            
    def saveFile(self):
        if not(self.Flag):
            QtGui.QMessageBox.warning(self, "Attention", "No data file loaded: chose a file to open")
            return
        else:
            self.myDefaultOutput()
            self.outputname = QtGui.QFileDialog.getSaveFileName(self, 'Save File', self.outpath,)
            self.fileInfo3.setText("Output file: " + self.outputname)
            if self.comp >= self.Nc:
                self.comp = 0
                message = 'File has only %i components\nThe default (x) will be used instead' % self.Nc
                QtGui.QMessageBox.warning(self, "Attention", message)
            np.savetxt( str( self.outputname) ,self.myField[0,:,:,self.comp],fmt='%15.14e')
            print 'done'
            
    def openFile(self):
        self.filename = QtGui.QFileDialog.getOpenFileName(self, 'Open File', self.inputpath)
        self.Flag = self.test_file(self.filename)
        if self.Flag:
            self.myField = self.analize_field(self.filename)
            self.statusDisplay.setText("Status OK")
            self.fileDisplay.setText("Input file: " + self.filename)
            self.changeCompChoice()
            self.addFileInfo()
        else:
            self.fileDisplay.setText("ERROR")
            
    def printName(self):
        print self.filename
        
    def test_file(self, filename):
        f  = open(filename ,'rb')
        #- endianness 0:small - 1:big -#
    
        endianness = struct.unpack('i', f.read(4))[0]

        #- global grid dim -#
        Nx = struct.unpack('i', f.read(4))[0]
        Ny = struct.unpack('i', f.read(4))[0]
        Nz = struct.unpack('i', f.read(4))[0]
    
        if not(Nx > 0 and Ny > 0 and Nz > 0):
            return False
        if not(Nx < 1000000 and Ny < 10000 and Nz <10000):
            return False
    
        #- processor grid -#
        Npx = struct.unpack('i', f.read(4))[0]
        Npy = struct.unpack('i', f.read(4))[0]
        Npz = struct.unpack('i', f.read(4))[0]

        Nproc = Npx*Npy*Npz
    
    
        if not(Npx > 0 and Npy > 0 and Npz > 0):
            return False
        if not(Npx < 1000 and Npy < 1000 and Npz <1000):
            return False
    
        #- field components -#
        Nc = struct.unpack('i', f.read(4))[0]
        if not(Nc < 4 and Nc >0):
            return False
        
        return True
        f.close()
    
    def analize_field(self,filename):
        f  = open(filename ,'rb')
        #- endianness 0:small - 1:big -#
    
        endianness = struct.unpack('i', f.read(4))[0]

        #- global grid dim -#
        Nx = struct.unpack('i', f.read(4))[0]
        Ny = struct.unpack('i', f.read(4))[0]
        Nz = struct.unpack('i', f.read(4))[0]
        
        Ntot = Nx*Ny*Nz
        self.Ntot = Ntot
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        #- processor grid -#
        Npx = struct.unpack('i', f.read(4))[0]
        Npy = struct.unpack('i', f.read(4))[0]
        Npz = struct.unpack('i', f.read(4))[0]

        Nproc = Npx*Npy*Npz
    
        #- field components -#
        self.Nc = struct.unpack('i', f.read(4))[0]
    
        
        #- grid -> X -#
        #x = np.zeros((Nx))
        #for i in range(0,Nx):
        x = struct.unpack('f'*Nx, f.read(4*Nx))
        np.savetxt( 'x.dat' ,x[:],fmt='%15.14e')            

        #- grid -> Y -#
        y = np.zeros((Ny))
        for i in range(0,Ny):
            y[i] = struct.unpack('f', f.read(4))[0]

        #- grid -> Z -#
        z = np.zeros((Nz))
        for i in range(0,Nz):
            z[i] = struct.unpack('f', f.read(4))[0]

        #- loop on processors -#
        F = np.zeros((Nz,Ny,Nx,self.Nc))
        counter = 0
        prog = 0.0
        for nprocessor in range(0,Nproc):
        #-processor dims -#
            i0  = struct.unpack('i', f.read(4))[0]
            j0  = struct.unpack('i', f.read(4))[0]
            k0  = struct.unpack('i', f.read(4))[0]
            li0 = struct.unpack('i', f.read(4))[0]
            lj0 = struct.unpack('i', f.read(4))[0]
            lk0 = struct.unpack('i', f.read(4))[0]
        #print '>>> ',i0,j0,k0,li0,lj0,lk0
            NN=li0*lj0*lk0*self.Nc
            array=np.array(struct.unpack('f'*NN, f.read(4*NN))).reshape(lk0,lj0,li0,self.Nc)
    
            for k in range(0,lk0):
                for j in range(0,lj0):
                    for i in range(0,li0):
                        for c in range(0,self.Nc):
                            F[k+k0,j+j0,i+i0,c] = array[k,j,i,c]
                    counter += li0
                    prog = counter*(100.0/Ntot)
                    self.progress.setValue(prog)
        #np.savetxt( nameOutFile ,F[0,:,:,component],fmt='%15.14e')            

        f.close()
        print "done"
        return F
    
        
def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())

run()
