#! /usr/bin/env python

import os
import sys
import re
import numpy as np

Nspins = 16
deph = 0
Omega = 1
Gamma = 0.01
typeint = 0
interactions = np.arange(1.0, 1.1, 1.0)
realisations = np.arange(2, 1001, 1)
detunings = np.arange (0.0, 0.1, 1.0)
gammabars = np.arange (0.0, 0.1, 1.0)
dephasingrates = np.arange (0.0, 0.1, 1.0)
#rads = np.arange (0.1, 1.01, 0.01) 

#def buildSubmit(filename, Nspins, gammabar, dephasingrate, detuning, interaction, deph, realisation):
def buildSubmit(filename, Nspins, gammabar, Omega, detuning, interaction, deph, Gamma, dephasingrate, typeint, realisation):

    inlines  = file('submit.sh').readlines()
    outlines = []
    for line in inlines:
        if 'PARAMETERS' in line:
            #line = line.replace('PARAMETERS=(0 0 0 0 0 0 0)', 'PARAMETERS=(%d %lf %lf %lf %lf %d %d)' % (Nspins, gammabar, dephasingrate, detuning, interaction, deph, realisation))
            line = line.replace('PARAMETERS=(0 0 0 0 0 0 0 0 0 0)', 'PARAMETERS=(%d %lf %lf %lf %lf %d %lf %lf %d %d)' % (Nspins, gammabar, Omega, detuning, interaction, deph, Gamma, dephasingrate, typeint, realisation))
 
        outlines.append(line)
    file(filename, 'w').writelines(outlines)   

for interaction in interactions:
    for dephasingrate in dephasingrates :
	for gammabar in gammabars :
	    for detuning in detunings:
	        for realisation in realisations:
			#print 'subjob run = no%d_g%.2lf_xi%.2lf_det%.2lf_R%.2lf_noise%d_run%d' % (Nspins, gammabar, dephasingrate, detuning, interaction, deph, realisation)
			print 'subjob run = no%d_f%.2lf_om%.2lf_det%.2lf_V%.2lf_n%d_gdec%.1lf_w%.1lf_v%d_run%d' % (Nspins, gammabar, Omega, detuning, interaction, deph, Gamma, dephasingrate, typeint, realisation)
			filename_submit = 'job_no%d_f%.2lf_om%.2lf_det%.2lf_V%.2lf_n%d_gdec%.1lf_w%.1lf_v%d_run%d.sh' % (Nspins, gammabar, Omega, detuning, interaction, deph, Gamma, dephasingrate, typeint, realisation)

			#buildSubmit(filename_submit, Nspins, gammabar, dephasingrate, detuning, interaction, deph, realisation)
			buildSubmit(filename_submit, Nspins, gammabar, Omega, detuning, interaction, deph, Gamma, dephasingrate, typeint, realisation)
			os.system('qsub ' + filename_submit)

