import sys
import numpy
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import os
import os.path

path = os.getcwd()

# Parameters and number of subensembles 
kB = 1.0
T=0.8
Numeps = 8
NumIntermediates = 15
NumIterations = 145000

# Allocate storage for simulation data
N_k = numpy.zeros([Numeps], numpy.int32)
eps = numpy.zeros([Numeps], numpy.float64)
E_kn = numpy.zeros([Numeps,NumIterations], numpy.float64)
z_kn = numpy.zeros([Numeps,NumIterations], numpy.float64)

# Read in surface strength
infile = open(os.path.join(path,'surface.dat'), 'r')
lines = infile.readlines()
infile.close()
for k in range(Numeps):
    eps[k] = lines[k]

# Read the simulation data
for k in range(Numeps):
    filename = os.path.join(path,'energy%d.out' %k)
    print ("Reading %s..." % filename)
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    n = 0
    for line in lines:
         tokens = line.split()
         z = float(tokens[1])
         E = float(tokens[2])
         z_kn[k,n] = z
         E_kn[k,n] = E     
         n += 1
    N_k[k] = n

beta = 1 / (kB * T)

# Insert intermediate epsilon and corresponding blank energy
delta = (eps[Numeps - 1] - eps[0]) / (NumIntermediates - 1)

val_k = []
currentv = eps[0]

while(currentv <= eps[Numeps - 1]):
	if not any( currentv == a for a in eps):
		val_k = numpy.append(val_k,currentv)
	currentv = currentv + delta
eps = numpy.concatenate((eps,numpy.array(val_k)))

K = len(eps)
print(eps)
print(K)
Nall_k = numpy.zeros([K],numpy.int32)
Eall_kn = numpy.zeros([K,NumIterations],numpy.float64) 

for k in range (Numeps):
	Eall_kn[k,0:N_k[k]] = E_kn[k,0:N_k[k]]
	Nall_k[k] = N_k[k]

# Compute reduced potential energies
print "--Computing reduced energies..."

u_kln = numpy.zeros([K,K,NumIterations], numpy.float64)

for k in range(K):
	for l in range(K):
		u_kln[k,l,0:Nall_k[k]] = beta * ((eps[l]) / eps[k]) * Eall_kn[k,0:Nall_k[k]] 

# Initialize MBAR.
print ("Running MBAR...")
mbar = pymbar.MBAR(u_kln, Nall_k, verbose=True, relative_tolerance=1e-12)