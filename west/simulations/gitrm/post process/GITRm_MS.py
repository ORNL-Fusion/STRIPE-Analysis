import sys, os
import numpy as np
import io, libconf
import math
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import LogNorm
import netCDF4 as nc
from os import path


file = np.array(['../diiid-helicon/DIII-D_helicon_runs/gitrm-surface.nc'])
#The following parameters indicate the final no of bins from the refined data; should be divisible by refined binning dimesnsions
coarsen=1
rows=9          #Final no of rows i.e Angles
columns=10      #Final no of comlumns i.e Energies

#Number of Files and species
num_files=1
n_species=1

#Names in the figures
f_names=np.array(["Carbon", "Tungsten"])

E_show=np.zeros((num_files, n_species, rows, columns))
E_g=np.zeros((num_files, n_species))


def plotAll(string, max_e):
	
	print("Plotting  " + string )
	fig = plt.figure(figsize = (16,10))
	plt.rcParams.update({'font.size': 22})
	r=[4.721309631429217e+15 , 0.089862316980425*1e13]	
	for f in range(num_files):
		
		ds_all = nc.Dataset(file[f])	
		dimE_col=len(ds_all.dimensions['nEnergies'])
		dimA_row=len(ds_all.dimensions['nAngles'])
		
		for s in range (n_species):
			edist_per_species  = ds_all[string][0,s,:,:,:]
			E_dist=np.sum(edist_per_species,axis=0).T
			if (bool(coarsen)):
				array_coarse=np.zeros((rows,columns))
				a=int(dimA_row/rows)
				b=int(dimE_col/columns)
				for i in range(rows):
					for j in range(columns):
						for k in range(a):
							for l in range(b):
								array_coarse[i][j]=array_coarse[i][j]+E_dist[i*a+k][j*b+l];
				E_show[f,s,:,:] =array_coarse
			
			plt.subplot(num_files,n_species, f*n_species+s+1)
			E_plot=np.copy(E_show[f,s,:,:])
			E_plot[E_plot<=0]=np.nan
			cs0=plt.imshow(r[s]*E_plot, extent=[0, 1000, 0, 90], cmap='cool', origin='lower', 
				aspect='auto', interpolation='nearest')
			plt.xlabel('Energy [eV]')
			plt.ylabel('Angle [deg]')
			C = plt.colorbar(cs0, pad=0.01, format='%.3g')
			#plt.clim(0,max_e)
			given_str=f_names[s]
			plt.title(given_str, weight="bold")

	plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.9, wspace=0.4, hspace=0.6)
	fig.savefig(string+"A2_2M.jpg", dpi=500)
	
	if num_files>1:	
		for s in range(n_species):
			print("   Species: norm-diff", s, 100*LA.norm((E_show[0,s,:,:]-E_show[1,s,:,:]))/LA.norm(E_show[0,s,:,:]) )

def plotLine(string):
	
	print("Plottinng " + string)
	fig = plt.figure(figsize = (23,10))
	plt.rcParams.update({'font.size': 28})	
	a1=['-k', '-k']
	a2=['k', 'k']
	r=[4.721309631429217e+15 , 0.089862316980425*1e13]
        	
	for f in range(num_files):
		ds_all = nc.Dataset(file[f])
		dimE=len(ds_all.dimensions['nEnergies'])
		n_bins=60
		array_coarse=np.zeros((n_bins, 1))
	
		for s in range (n_species):
			edist_per_species  = ds_all[string][s,:,:,:]
			E_dist=np.sum(edist_per_species,axis=0)
			E_dist=np.sum(E_dist, axis=1)
			for i in range(n_bins):
				array_coarse[i]=np.sum(E_dist[int(i*dimE/n_bins):int((i+1)*dimE/n_bins)])
			print(np.sum(E_dist), np.sum(array_coarse))
			plt.subplot(1, 2, s+1)	
			plt.plot(np.linspace(0, 100, n_bins), r[s]*array_coarse, a1[s], markerfacecolor=a2[s],markeredgecolor=a2[s],markersize=6)
			plt.title(f_names[s], weight="bold")
			#plt.xlim([0,90])
			#plt.ylim([0, 2000])
			plt.xlabel("Energy (eV)")
			plt.ylabel(r"Particle rate $(s^{-1})$", fontsize=22)
		plt.savefig("E1_2M.jpg", dpi=500)
	
def plotLine2(string):
	
	print("Plottinng " + string)
	fig = plt.figure(figsize = (23,10))
	plt.rcParams.update({'font.size': 28})	
	a1=['-k', '-k']
	a2=['k', 'k']
	r=[4.721309631429217e+15 , 0.089862316980425*1e13]
        	
	for f in range(num_files):
		ds_all = nc.Dataset(file[f])
		dimE=len(ds_all.dimensions['nEnergies'])
		dimA=len(ds_all.dimensions['nAngles'])	
		for s in range (n_species):
			edist_per_species  = ds_all[string][s,:,:,:]
			E_dist=np.sum(edist_per_species,axis=0)
			E_dist=np.sum(E_dist, axis=0)
			print(np.shape(E_dist))
			plt.subplot(1, 2, s+1)
			plt.plot(np.linspace(0, 90, 90), r[s]*E_dist, a1[s], markerfacecolor=a2[s], markeredgecolor=a2[s],markersize=6)
			#plt.xlim([0,90])
			#plt.ylim([0, 2000])
			plt.xlabel("Angle (degrees)")
			plt.ylabel(r"Particle rate $(s^{-1})$", fontsize=22)
			plt.title(f_names[s], weight="bold")
		plt.savefig(string + "MS_.jpg", dpi=500)
	

def global_values(string):
	print("Calculating  " + string)
	for f in range(num_files):
		for s in range (n_species):
			ds_all = nc.Dataset(file[f])
			erosion = ds_all[string][s,:]
			E_g[f,s]=np.sum(erosion)
	
	if num_files>1:
		for s in range(n_species):
			print("   Species:  diff", s, E_g[0,s], E_g[1,s], 100*(E_g[0,s]-E_g[1,s])/E_g[0,s])

	
def main():
	plotAll('surfEDist', 5.5e5)
	#plotAll('surfReflDist', 3e4)
	plotLine2('surfReflDist')
	plotLine('surfSputtDist')
	#global_values('grossErosion')
	#global_values('grossDeposition')

if __name__=="__main__":
	main()
