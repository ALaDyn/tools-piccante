#!/usr/bin/python

### loading shell commands
import os, os.path, glob, sys, shutil, time
from scipy.fftpack import fftn
#from PyQt4.QtGui import *
import struct
#from scipy import *
import numpy as np
import math

# from matplotlib import *
# from pylab import *
# import matplotlib as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
### --- ###

class field_analysis:
    
    def __init__(self,tmin, tmax,substeps, base, end):
        self.base = base
        self.end = end
        self.tmin = int(tmin)
        self.tmax = int(tmax)
        self.substeps = int(substeps)
        
        if 1000%substeps:
            print "ERROR: wrong number of timesteps!"
            sys.exit()
        self.Nt = substeps*(self.tmax-self.tmin)
        self.filenames = list()
        self.t = np.zeros((self.Nt))
        self.chose_file()
        
        self.init_values()
        self.set_frequency()
        self.basename = ("%s%03d"% (self.base,self.tmin))
        
        
    def analize_field(self, filename):
        f  = open(filename ,'rb')
        #- endianness 0:small - 1:big -#
    
        endianness = struct.unpack('i', f.read(4))[0]

        #- global grid dim -#
        Nx = struct.unpack('i', f.read(4))[0]
        Ny = struct.unpack('i', f.read(4))[0]
        Nz = struct.unpack('i', f.read(4))[0]
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        Ntot = Nx*Ny*Nz
        #- processor grid -#
        Npx = struct.unpack('i', f.read(4))[0]
        Npy = struct.unpack('i', f.read(4))[0]
        Npz = struct.unpack('i', f.read(4))[0]

        Nproc = Npx*Npy*Npz
    
        #- field components -#
        Nc = struct.unpack('i', f.read(4))[0]
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Nc = Nc
        
        
        #- grid -> X -#
        #x = np.zeros((Nx))
        #for i in range(0,Nx):
        x = struct.unpack('f'*Nx, f.read(4*Nx))
        
        #- grid -> Y -#
        y = struct.unpack('f'*Ny, f.read(4*Ny))
        self.x = x
        #- grid -> Z -#
        z = struct.unpack('f'*Nz, f.read(4*Nz))
        self.x = x
        self.y = y
        self.z = z
        #- loop on processors -#
        F = np.zeros((Nz,Ny,Nx,Nc))
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
            NN=li0*lj0*lk0*Nc
            array=np.array(struct.unpack('f'*NN, f.read(4*NN))).reshape(lk0,lj0,li0,Nc)
    
            for k in range(0,lk0):
                for j in range(0,lj0):
                    for i in range(0,li0):
                        for c in range(0,Nc):
                            F[k+k0,j+j0,i+i0,c] = array[k,j,i,c]
                            counter += li0
                            prog = counter*(100.0/Ntot)
                    
        #np.savetxt( nameOutFile ,F[0,:,:,component],fmt='%15.14e')            
        f.close()
        #print "done"
        return F

    def collect_file(self):
        name = "E_FIELD_000.000.bin.000"
        for i in os.listdir(os.getcwd()):
            if i.endswith(".bin.000") and i.startswith("E_FIELD_"): 
                #print i
                self.Nt +=1
                continue
            else:
                continue
    def chose_file(self):
        
        base = self.base
        end = self.end
        tindex = 0
        for i in range (self.tmin,self.tmax):
            for v in range (0,self.substeps):
                t = i + v*1.0/self.substeps
                
                self.t[tindex] = t
                name = base + ("%03d.%03d" % (i,v*1000/self.substeps) ) + end
                self.filenames.append(name)
                #print name
                tindex += 1
        self.Nt = len(self.filenames)
        
    def collect_data(self):
        self.alldata = np.zeros((self.Nt,self.Nz,self.Ny,self.Nx,self.Nc))
        for i in range(0,self.Nt):
            #print self.filenames[i]
            F = self.analize_field(self.filenames[i])
            self.alldata[i,:,:,:,:] = F[:,:,:,:]
            
            
    
    def init_values(self):
        self.analize_field(self.filenames[0])
        
        if(self.Nx>1):
            self.dx = self.x[1] - self.x[0]
        else:
            self.dx = 0
        if(self.Ny>1):
            self.dy = self.y[1] - self.y[0]
        else:
            self.dy = 0
        if(self.Nz>1):
            self.dz = self.z[1] - self.z[0]
        else:
            self.dz = 0
        if(self.Nt>1):
            self.dt = self.t[1] - self.t[0]
        else:
            self.dt = 0
        self.kx = np.empty_like(self.x)
        self.ky = np.empty_like(self.y)
        self.kz = np.empty_like(self.z)
        self.kt = np.empty_like(self.t)
        self.Lx = self.dx*self.Nx
        self.Ly = self.dy*self.Ny
        self.Lz = self.dz*self.Nz
        self.Lt = self.dt*self.Nt
        self.print_parameters()
    
    def print_parameters(self):
        print ("SIZE: [ Lx, Ly, Lz ] = [ %f, %f, %f ]" %(self.Lx,self.Ly,self.Lz))
        print (" Np : [ Nx, Ny, Nz ] = [ %d, %d, %d ]" %(self.Nx,self.Ny,self.Nz))
        print ("time span = [%f:%f]     Nt = %d" %(self.tmin, self.tmax,self.Nt))
        
    def set_frequency(self):
        if(self.Lx>0):
            self.dkx = 1/self.Lx
        else:
            self.dkx = 0
        for i in range(0,self.Nx):
            self.kx[i] = i*self.dkx

        if(self.Ly>0):
            self.dky = 1/self.Ly
        else:
            self.dky = 0
        for i in range(0,self.Ny):
            self.ky[i] = i*self.dky

        if(self.Lz>0):
            self.dkz = 1/self.Lz
        else:
            self.dkz = 0
        for i in range(0,self.Nz):
            self.kz[i] = i*self.dkz

        if(self.Lt>0):
            self.dkt = 1/self.Lt
        else:
            self.dkt = 0
        for i in range(0,self.Nt):
            self.kt[i] = i*self.dkt

        self.Lkx = self.dkx*(self.Nx/2)
        self.Lky = self.dky*(self.Ny/2)
        self.Lkz = self.dkz*(self.Nz/2)
        self.Lkt = self.dkt*(self.Nt/2)
        self.print_frequency()
        
    def print_frequency(self):
        print ("SIZE: [ Lkx, Lky, Lkz ] = [ %f, %f, %f ]" %(self.Lkx,self.Lky,self.Lkz))
        print (" dk : [ dkx, dky, dkz ] = [ %f, %f, %f ]" %(self.dkx,self.dky,self.dkz))
        print ("max_freq = %f     domega = %f" %(self.Lkt,self.dkt))
        
    def do_fft(self,zposition,comp):
        
        self.trasf = np.zeros((self.Nt,self.Ny,self.Nx))
        
        self.trasf = np.fft.fftn(self.alldata[:,zposition,:,:,comp])

    def select_omega(self,freq):
         ifreq = int(freq/self.dkt)
         selected=ifreq*self.dkt
         print ("selected ifreq = %d, freq=%4.3f" %(ifreq,selected))
         
         name = ("%s-kx-ky-omega%4.3f.txt"% (self.basename,selected))
         f1=open(name, 'w')
         for j in range(0,self.Ny/2):
             for i in range(0,self.Nx/2):
                 f1.write("%e %e %e\n" % (self.kx[i], self.ky[j], np.real(self.trasf[ifreq,i,j])))
         #np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
         f1.close()
         
    def kx_ky_oneTime(self,time,zposition,comp):
         itime = int((time-self.tmin)/self.dt)
         if itime>=self.Nt:
             itime = self.Nt-1
         selected=itime*self.dt + self.tmin
         
         print ("selected itime = %d  time = %4.3f" %(itime,selected))
         name = ("%s-kx-ky-time%3.3f.txt"% (self.basename,selected))
         
         onetimetrasf = np.zeros((self.Ny,self.Nx))
        
         onetimetrasf = np.fft.fftn(self.alldata[itime,zposition,:,:,comp])
         f1=open(name, 'w')
         
         for j in range(-self.Ny/2,self.Ny/2):
             ky = j*self.dky
             for i in range(-self.Nx/2,self.Nx/2):
                 kx = i*self.dkx
                 f1.write("%e %e %e\n" % (kx, ky, np.real(onetimetrasf[i,j])))
         f1.close()
         
    def avg_y(self):
         name = ("%s-kx-omega.txt"% (self.basename,))
        
         f1=open(name, 'w')
         for j in range(-self.Nt/2,self.Nt/2):
             kt = j*self.dkt
             for i in range(-self.Nx/2,self.Nx/2):
                 kx = i*self.dkx
                 f1.write("%e %e %e\n" % (kx, kt, np.real(np.sum(self.trasf[j,i,:]))))
         #np.savetxt( "kx-omega.txt" ,np.real(self.trasf[:,:,0]),fmt='%15.14e')
         f1.close()
    
def run():
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)
    path     = os.getcwd()
    #f        = open(os.path.join(path,'E_FIELD_000.000.bin.000'),'rb')
    component = 0
    basetime = 100
    deltaT = 10
    sub = 5
    endtime = basetime + deltaT
    #myname = "B_FIELD_"
    myname = "DENS_eleBulk_"
    myanalysis = field_analysis(basetime, endtime,substeps=sub, base=myname, end=".bin.000")
    myanalysis.collect_data()
    myanalysis.do_fft(zposition=0,comp=component)
    myanalysis.select_omega(freq=1.01)
    myanalysis.select_omega(freq=2.01)
    
    myanalysis.avg_y()
    myanalysis.kx_ky_oneTime(basetime,zposition=0,comp=component)
run()
    
def old():
    if len(sys.argv)<2:
        sys.exit('Usage: %s inputFile outputFile' % sys.argv[0])
    
    elif len(sys.argv)<3:
        outFile = 'default.dat'
    else:
        outFile = str(sys.argv[2])

    