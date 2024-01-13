#!/bin/sh
################################################################################
################################################################################
# This script submits an array job to EDDIE
#
# It requires the beauti XMLs to be present in the current working directory and
# for all beast runs to require the same command line options.
#
# The options for beast and grid engine for each run are specified in beast.sh
#
# To submit: qsub -t 1-n arrayjob.sh, in which n = the total number of XMLs to 
# run.

################################################################################
################################################################################
# Grid Engine options (lines prefixed with #$)
#$ -N beast_array_submission
#$ -cwd                
#$ -l h_vmem=2G

################################################################################
xml_files=$(find /exports/eddie/scratch/jbaxter/ -type f -name "*.xml")

./beast.sh $xml_files 

################################################################################
################################################################################
# END #
################################################################################
################################################################################