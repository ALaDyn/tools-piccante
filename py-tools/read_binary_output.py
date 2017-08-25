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


######################################################################
# Name:         read_binary_output.py
# Author:       
__author__ = "Alberto Marocchino"
__credits__ = ["whomever", "you", "want"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alberto Marocchino"
__status__ = "Beta"
# Date:			2014-05-29
# Purpose:      reads binary frm output
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import struct
#from scipy import *
import numpy as np
# from matplotlib import *
# from pylab import *
# import matplotlib as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
### --- ###

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
path     = os.getcwd()
f        = open(os.path.join(path,'E_FIELD_000.000.bin'),'rb')

#- endianness 0:small - 1:big -#
endianness = struct.unpack('i', f.read(4))[0]

#- global grid dim -#
Nx = struct.unpack('i', f.read(4))[0]
Ny = struct.unpack('i', f.read(4))[0]
Nz = struct.unpack('i', f.read(4))[0]

#- processor grid -#
Npx = struct.unpack('i', f.read(4))[0]
Npy = struct.unpack('i', f.read(4))[0]
Npz = struct.unpack('i', f.read(4))[0]

Nproc = Npx*Npy*Npz
#- field components -#
Nc = struct.unpack('i', f.read(4))[0]

#- grid -> X -#
x = np.zeros((Nx))
for i in range(0,Nx):
	x[i] = struct.unpack('f', f.read(4))[0]

#- grid -> Y -#
y = np.zeros((Ny))
for i in range(0,Ny):
	y[i] = struct.unpack('f', f.read(4))[0]

#- grid -> Z -#
z = np.zeros((Nz))
for i in range(0,Nz):
	z[i] = struct.unpack('f', f.read(4))[0]

#- loop on processors -#
F = np.zeros((Nz,Ny,Nx,Nc))
for nprocessor in range(0,Nproc):
		#-processor dims -#
		i0  = struct.unpack('i', f.read(4))[0]
		j0  = struct.unpack('i', f.read(4))[0]
		k0  = struct.unpack('i', f.read(4))[0]
		li0 = struct.unpack('i', f.read(4))[0]
		lj0 = struct.unpack('i', f.read(4))[0]
		lk0 = struct.unpack('i', f.read(4))[0]
		print '>>> ',i0,j0,k0,li0,lj0,lk0
        locFields = np.zeros((li0,lj0,lk0,Nc))
        NN=li0*lj0*lk0*Nc
        array=numpy.array(struct.unpack('f'*NN, f.read(4*NN))).reshape(lk0,lj0,li0,Nc)
        
		for k in range(0,lk0):
			for j in range(0,lj0):
				for i in range(0,li0):
					for c in range(0,Nc):
						F[k+k0,j+j0,i+i0,c] = array[k,j,i,c]
                        
np.savetxt( 'F0.dat' ,F[0,:,:,1],fmt='%15.14e')			

f.close()
print path







# np.savetxt( 'F1.dat' ,F1[:,:,0],fmt='%15.14e')			
# np.savetxt( 'F2.dat' ,F2[:,:,0],fmt='%15.14e')			


#- processor -#
# for i in range(0,6):
# 	print struct.unpack('i', f.read(4))[0]
# 
# for i in range(0,360*64*1+1):
# 	print struct.unpack('f', f.read(4))[0]



#- endianess 0:small - 1:big -#
#struct.unpack('i', f.read(4))



# 
# # - #
# def read_ALaDyn_bin(dir_path,file_name,grid_no_grid):
# 
# 	# - #
# 	path     = os.path.join(os.path.join(dir_path,file_name))
# 	f        = open(path,'rb')
# 
# 	#- vector length -#
# 	struct.unpack('i', f.read(4))
# 	N_param = struct.unpack('i', f.read(4))[0]
# 	print N_param
# 	struct.unpack('i', f.read(4))
# 
# 
# 	struct.unpack('i', f.read(4))
# 	int_param=[]
# 	for i in range(0,N_param):
# 		int_param.append( struct.unpack('i', f.read(4)) )
# 	struct.unpack('i', f.read(4))
# 	nx= int_param[3][0]
# 	ny= int_param[4][0]
# 	nz= int_param[6][0]
# 	nproc_y = int_param[0][0]
# 	nproc_z = int_param[1][0]
# 
# 	struct.unpack('i', f.read(4))
# 	for i in range(0,N_param):
# 		struct.unpack('f', f.read(4))
# 	struct.unpack('i', f.read(4))
# 
# 
# 	#---***---#
# 	r = np.zeros((nx,ny,nz))
# 	print 'total grid size: n=(',nx,ny,nz,')'
# 	rr=[]
# 	print 'number of Np processors: Np_y=',nproc_y,'Np_z=', nproc_z
# 
# 	offsetz = 0
# 	for counter_z in range(0,nproc_z):
# 		offsety = 0;
# 
# 		for counter_y in range(0,nproc_y):
# 			struct.unpack('i', f.read(4))
# 			npx= struct.unpack('i', f.read(4))[0]
# 			npy= struct.unpack('i', f.read(4))[0]
# 			npz= struct.unpack('i', f.read(4))[0]
# 			struct.unpack('i', f.read(4))
# 	
# 			struct.unpack('i', f.read(4))
# 			for k in range(0,npz):
# 				for j in range(0,npy):
# 					for i in range(0,npx):
# 						r[i,j+offsety,k+offsetz] = struct.unpack('f', f.read(4))[0]
# 			struct.unpack('i', f.read(4))
# 			offsety += npy
# 		offsetz += npz;
# 		
# 	if grid_no_grid == 'nogrid':
# 		return r
# 
# 	#--- * --- * --- * --- * --- * ---#
# 	#- reading grid -#
# 	struct.unpack('i', f.read(4))
# 	X=[]; [X.append(struct.unpack('f', f.read(4))[0]) for i in range(0,nx)]	
# 	struct.unpack('i', f.read(4))
# 
# 	struct.unpack('i', f.read(4))
# 	Y=[];  [Y.append(struct.unpack('f', f.read(4))[0]) for i in range(0,ny)]
# 	struct.unpack('i', f.read(4))
# 
# 	struct.unpack('i', f.read(4))
# 	Z=[];  [Z.append(struct.unpack('f', f.read(4))[0]) for i in range(0,nz)]
# 	struct.unpack('i', f.read(4))
# 
# # 	x=np.zeros((nx,ny,nz)); y=z=x;
# # 	for k in range(0,nz):
# # 		for j in range(0,ny):
# # 			for i in range(0,nx):
# # 				x[i,j,k] = X[i]
# # 				y[i,j,k] = Y[j]
# # 				z[i,j,k] = Z[k]
# 	x=X; y=Y; z=Z;
# 				
# 	return (r,x,y,z)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
