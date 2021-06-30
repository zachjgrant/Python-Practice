#!/usr/bin/env python

#Python apadtation to capstone code with attempted MPI
from mpi4py import MPI
import numpy as np
import math
import statistics
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

start_time = time.time()
P = 22.5 #Nominal Power of Lamp at End of Life (W)
R = 3.81 #Chamber Radius (cm)
R_L = 2.54 #Lamp Radius (cm)
L = 50 #Length of Lamp (cm)
T_water = .95 #UV Transmittance
refl_SS = .30 #Max UV Reflectance for SS
T_Silica = .95 #Max UV Transmittance for Quartz Silica Sleeve
t = .3 #Thickness of quartz sleeve (cm)
alpha = math.log(T_water)
I_0 = T_Silica * (P / (2 * R_L * math.pi * L)) #Initial Intensity (W/cm^2)

x = np.linspace(R_L + (((R-R_L) / size)*rank), R_L + (((R-R_L) / size)*(rank+1)-1e-7/size),12700000/size)
#x = np.linspace(R_L, R, 12700001) #cm
I_xout = np.zeros(int(12700000/size), dtype = float)
I_xback = np.zeros(int(12700000/size), dtype = float)

I_xout = np.array(I_0 * (R_L / x) * np.exp(alpha * (x-R_L))) #Intensity for Outward Propogation (W/cm^2)
I_wall = I_0 * (R_L / R) * math.exp(alpha * (R-R_L)) #Intensity at the Wall (W/cm^2)
I_0back = refl_SS * I_wall #Initial Intensity for Inward Propogation (W/cm^2)
I_xback = np.array(I_0back * (R / x) * np.exp(alpha * (abs(x-R)))) #Intensity for Inward Propogation (W/cm^2)
I_tot = I_xout + I_xback
I_mean = np.mean(I_tot)

comm.Barrier()
I_mean = comm.gather(I_mean, root = 0)
#I_xback = comm.gather(I_xback, root = 0)
#I_xout1 = np.zeros(int(12700000), dtype = float)
#I_xback1 = np.zeros(int(12700000), dtype = float)
#num = np.zeros(size, dtype = float) + 12700000/size
#dis = np.linspace(0, 12700000-(12700000/size), size)

#comm.Gatherv(I_xout, [I_xout1, num, dis, MPI.FLOAT])
#comm.Gatherv(I_xback, [I_xback1, num, dis, MPI.FLOAT])

if rank == 0:
	#I_total = I_xout1 + I_xback1 #Total Intensity (W/cm^2)
	I_mean2 = np.mean(I_mean) #Average Intensity (W/cm^2) (Can assume this for turbulent flow)


	UVD = .016 #UV Dosage Standard (J/cm^2)
	t = UVD / I_mean2 #Time of treatment needed (s)

	A_c = (R**2 - R_L**2) * math.pi #Chamber Area (cm^2)
	v = L / t #Inlet Velocity (cm/s)
	Q = A_c * v * .015850323 #flow rate

	print("---%f gpm---" %Q)

	print("---%s seconds---" %(time.time()-start_time))
	

