#!/bin/bash
# Submission script for my HMC master project to the schools compute cluster
# to submit this type qsub ./(name of script)
# (1) Run the job in the same directory as the qsub is issued 
#$ -cwd
#
# (2) Name the job 
#
#$ -N HMC-Fourier
#
#(3) send a notification email whent he job starts (b) and finishes (e) 
#
#$ -m be
#
# (4) set the queu which the job should run on. for undergrads at edinburgh we 
# can onlye acces the sopa.1.day queue
#
#$ -q sopa.1.day
#
# Estimate of how long the job should take...for this usually 15 mins should be
# enough apaprt from final runs hrs:mins:secs
#
#$ -l h_rt=0:15:0
#

dir -p Cluster_Data
./HMC 
 
