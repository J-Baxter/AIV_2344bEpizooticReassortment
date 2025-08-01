#!/bin/bash

################################################################################
################################################################################
# Grid Engine options (lines prefixed with #$)
# Name job and set to use current working directory
#$ -N stagein 
#$ -cwd
#$ -q staging
#$ -l h_rt=12:00:00 
#$ -r yes
#$ -notify
trap 'exit 99' sigusr1 sigusr2 sigterm

################################################################################
# Source and destination directories
#
# Source path on DataStore in the staging environment
# Note: these paths are only available on the staging nodes
# It should start with one of /exports/csce/datastore, /exports/chss/datastore, /exports/cmvm/datastore or /exports/igmm/datastore
#
SOURCE=/exports/igmm/datastore/<SOURCE DIR>
#
# Destination path on Eddie. It should be on the fast Eddie HPC filesystem, starting with one of:
# /exports/csce/eddie, /exports/chss/eddie, /exports/cmvm/eddie, /exports/igmm/eddie or /exports/eddie/scratch, 
#
DESTINATION=/exports/eddie/scratch/jbaxter/<DESTINATION DIR>
  
# Perform copy with rsync
# Note: do not use -p or -a (implies -p) as this can break file ACLs at the destination
# rsync -rl ${SOURCE} ${DESTINATION}

# scp -r data/ jbaxter@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/jbaxter/
