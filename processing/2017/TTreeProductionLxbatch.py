#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

queue = "2nd" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 

"""
command = [	#'./DarkMatterSearchLeadingMuons -1 dataMuMu_2016 B 2016',
		#'./DarkMatterSearchLeadingMuons -1 dataMuMu_2016 C 2016',
		#'./DarkMatterSearchLeadingMuons -1 dataMuMu_2016 G 2016',
		'./DarkMatterSearchLeadingMuons -1 dataMuMu_2017 B 2017',
		'./DarkMatterSearchLeadingMuons -1 dataMuMu_2017 C 2017',
		'./DarkMatterSearchLeadingMuons -1 dataMuMu_2017 D 2017',
		'./DarkMatterSearchLeadingMuons -1 dataMuMu_2017 E 2017',
		'./DarkMatterSearchLeadingMuons -1 dataMuMu_2017 F 2017',
		#'./DarkMatterSearchLeadingElectrons -1 dataEE_2016 B 2016',
		#'./DarkMatterSearchLeadingElectrons -1 dataEE_2016 C 2016',
		#'./DarkMatterSearchLeadingElectrons -1 dataEE_2016 G 2016',
		'./DarkMatterSearchLeadingElectrons -1 dataEE_2017 B 2017',
		'./DarkMatterSearchLeadingElectrons -1 dataEE_2017 C 2017',
		'./DarkMatterSearchLeadingElectrons -1 dataEE_2017 D 2017',
		'./DarkMatterSearchLeadingElectrons -1 dataEE_2017 E 2017',
		'./DarkMatterSearchLeadingElectrons -1 dataEE_2017 F 2017'
	  ]	
"""

commandLines = [ './Dilepton -1 MuMu-v1 B 2017' ]

#print '\nSending Lxbatch jobs to produce NTuples for Missing Mass CTPPS Analysis\n\n'
print '\nSending Lxbatch jobs to produce NTuples for PPS Analysis\n\n'
path = os.getcwd()

index = 1
for line in commandLines:  
    jobFileName = 'job_' + str(index) + '.sh'
    with open(jobFileName, 'w') as fout:
	fout.write("#!/bin/sh\n")
	fout.write("echo\n")
	fout.write("echo 'START---------------'\n")
	fout.write("\n")
	#fout.write("source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh\n")
	fout.write("source /cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/gcc/6.3.0/etc/profile.d/init.sh\n")
	#fout.write("source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.06/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n")
	fout.write("source /cvmfs/cms.cern.ch/slc6_amd64_gcc630/lcg/root/6.10.08/bin/thisroot.sh\n")
	fout.write("cd "+str(path)+"\n")
	fout.write(line+"\n")
	fout.write("echo 'STOP---------------'\n")
	fout.write("echo\n")
    os.system( "chmod 755 " + jobFileName )
    subCommand = "bsub -q " + queue + " -o logs -J job" + str(index) + " < " + jobFileName
    print subCommand
    os.system( subCommand )
    print "job nr " + str(index) + " submitted"
    index += 1

print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
