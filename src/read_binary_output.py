#!/usr/bin/python
######################################################################
# Name:         read_binary_output.py
# Author:       
__author__ = "Alberto Marocchino"
__credits__ = ["whomever", "you", "want"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Alberto Marocchino"
__email__ = "albz_uk@gmail.com"
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
F = np.zeros((Nx,Ny,Nz,Nc))
for nprocessor in range(0,Npx*Npy*Npz):
		#-processor dims -#
		i0  = struct.unpack('i', f.read(4))[0]
		j0  = struct.unpack('i', f.read(4))[0]
		k0  = struct.unpack('i', f.read(4))[0]
		li0 = struct.unpack('i', f.read(4))[0]
		lj0 = struct.unpack('i', f.read(4))[0]
		lk0 = struct.unpack('i', f.read(4))[0]
		print '>>> ',i0,j0,k0,li0,lj0,lk0

		for k in range(k0,k0+lk0):
			for j in range(j0,j0+lj0):
				for i in range(i0,i0+li0):
					for c in range(0,Nc):
						F[i,j,k,c] = struct.unpack('f', f.read(4))[0]
# 					F1[i,j,k] = struct.unpack('f', f.read(4))[0]
# 					F2[i,j,k] = struct.unpack('f', f.read(4))[0]

np.savetxt( 'F0.dat' ,F[:,:,0,1],fmt='%15.14e')			

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
