#!/bin/bash
# Author: adpozuelo@gmail.com
# Version: 1.0
# Date: 03/2021
# Input files: in.mc and both log.mc.gpu.'units'.float & log.mc.gpu.'units' to process data

import matplotlib.pyplot as plt
import numpy as np
import csv

def readData(filename, a, b):
    with open(filename) as csvfile:
        myreader = csv.reader(csvfile, delimiter=' ')
        x = [] ; y = []
        for row in myreader:
            if len(row) > 1:
                if float(row[b]) > 0.0:
                    x.append((float(row[a])))
                    y.append((float(row[b])))
    return np.array(x), np.array(y)

def create_histogram_file(filename, esr):
    ehisto = np.zeros(nde * 2 + 1, dtype=int)
    for energy in esr:
        nei = int(np.rint(energy / delta_esr)) + nde;
        ehisto[nei] += 1
    htotal = np.sum(ehisto)
    fout = open(filename, 'w')
    for i in range(len(ehisto)):
        if ehisto[i] > 0:
            energy = (i - nde + 0.5) * delta_esr
            percentage = ehisto[i] / htotal
            fout.write('%g %g\n' % (energy, percentage))
    fout.close()

mc_file = open ('in.mc', 'r')
line_splitted = mc_file.readline().split()
nstep = line_splitted[0]
nequil = line_splitted[1]
units = line_splitted[3]
mc_file.readline() ; mc_file.readline()
line_splitted = mc_file.readline().split()
natoms = int(line_splitted[1])
line_splitted = mc_file.readline().split()
natoms += int(line_splitted[1])
mc_file.close()

delta_esr = 0.1
nde = 50000
ljToKelvin = 1.0
if units == 'K':
    delta_esr = 1000.0
    ljToKelvin = 120.0 * natoms

esr = []

is_simulation = False
with open('log.mc.gpu.' + units) as f:
    for line in f:
        line_splitted = line.split()
        if line_splitted[0] == nequil:
            is_simulation = True
        if is_simulation:
            esr.append(line_splitted[2])

esr = np.array(esr, dtype=float)
mc_ehisto_file_double = 'ehisto_mc_gpu_' + units + '.dat'
create_histogram_file(mc_ehisto_file_double, esr)

esr = []

is_simulation = False
with open('log.mc.gpu.' + units + '.float') as f:
    for line in f:
        line_splitted = line.split()
        if line_splitted[0] == nequil:
            is_simulation = True
        if is_simulation:
            esr.append(line_splitted[2])

esr = np.array(esr, dtype=float)
mc_ehisto_file_float = 'ehisto_mc_gpu_' + units + '.dat'
create_histogram_file(mc_ehisto_file_float, esr)

x1,y1 = readData(mc_ehisto_file_double, 0, 1)
x2,y2 = readData(mc_ehisto_file_float, 0, 1)
plt.scatter(x1,y1,label='Double', color='red')
plt.plot(x2,y2,label='Float')
plt.xlabel(units)
plt.ylabel('Percentage')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1, ncol=2, mode="expand",
           borderaxespad=0.)
plt.savefig('ehisto_float_vs_double_' + units + '.png')
plt.show()