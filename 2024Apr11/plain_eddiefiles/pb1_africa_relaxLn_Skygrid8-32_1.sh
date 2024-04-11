#!/bin/sh 
################################################################################
################################################################################
# This script runs a simple BEAST run on EDDIE 
# 
# It requires the beauti XMLs to be present in the current working directory and 
# for all beast runs to require the same command line options.
#
# Works alongside arrayjobsubmit.sh as part of array job submission, but can
# easily be appropriated for submission of a single BEAST run
################################################################################
# Grid Engine options
#$ -N pb1_africa_relaxLn_Skygrid8-32_1.xml
#$ -cwd
#$ -pe sharedmem 2
#$ -l h_vmem=8G
#$ -l h_rt=160:00:00
#$ -M james.baxter@ed.ac.uk
#$ -P roslin_eeid_aiv
#$ -m baes
. /etc/profile.d/modules.sh
################################################################################
# load BEAGLE and BEAST
module load roslin/beast/1.10.4-beagle2
################################################################################
# Run the program
echo '=============================================='
echo '** Hello BEAST user !**'
echo "This job is running on $HOSTNAME"
echo 'Start BEAST with AIV data this is a 160 hour run'
beast pb1_africa_relaxLn_Skygrid8-32_1.xml 
echo '** Done **'
echo "============================================="
################################################################################
################################################################################
# END #
################################################################################
################################################################################
You can write whatever you want in it.
This is just an example.
