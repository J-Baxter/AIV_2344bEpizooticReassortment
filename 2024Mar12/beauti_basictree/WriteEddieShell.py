# Writes a batch submission shell script to run BEAST (Rambaut et al. 2004) on a
# linux server running Oracle Grid Engine/Sun Grid Engine.
#
# Arguments: 1) path to beauti XML that you wish to use to run the BEAST analysis
#
# NB: Either run this script on an interactive node on the server, or make sure
#     you run in the same working directory as your beauti.xml. The path submitted 
#     to BEAST will be the path submitted in the first argument.
#
# Copyright (c) 2024 James Baxter

# import modules
import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(description="a script to write sun grid batch job submission for BEAST1")
    parser.add_argument("xml_path")
    args = parser.parse_args()
    return args


def get_filename(input_string):
    match = re.search(r'[^/]+$', input_string)
    if match:
        return match.group(0)
    else:
        return input_string


def main():
    inputs = parse_args()
    relative_path = inputs.xml_path
    shell_filename = re.sub('.xml$', '.sh', relative_path)
    xml_filename = get_filename(relative_path)

    
    with open(shell_filename, "w") as file:
        file.write('#!/bin/sh \n')
        file.write('################################################################################\n')
        file.write('################################################################################\n')
        file.write('# This script runs a simple BEAST run on EDDIE \n')
        file.write('# \n')
        file.write('# It requires the beauti XMLs to be present in the current working directory and \n')
        file.write('# for all beast runs to require the same command line options.\n')
        file.write('#\n')
        file.write('# Works alongside arrayjobsubmit.sh as part of array job submission, but can\n')
        file.write('# easily be appropriated for submission of a single BEAST run\n')
        file.write('################################################################################\n')
        file.write('# Grid Engine options\n')
        file.write('#$ -N ' + xml_filename + '\n')
        file.write('#$ -cwd\n')
        file.write('#$ -pe sharedmem 2\n')
        file.write('#$ -l h_vmem=8G\n')
        file.write('#$ -l h_rt=200:00:00\n')
        file.write('#$ -M james.baxter@ed.ac.uk\n')
        file.write('#$ -P roslin_eeid_aiv\n')
        file.write('#$ -m baes\n') 
        file.write('. /etc/profile.d/modules.sh\n')
        file.write('################################################################################\n')
        file.write('# load BEAGLE and BEAST\n')
        file.write('module load roslin/beast/1.10.4-beagle2\n')
        file.write('################################################################################\n')
        file.write('# Run the program\n')
        file.write("echo '=============================================='\n")
        file.write("echo '** Hello BEAST user !**'\n")
        file.write('echo "This job is running on $HOSTNAME"\n')
        file.write("echo 'Start BEAST with AIV data this is a 168 hour run'\n")
        file.write('beast '+ xml_filename + ' \n')
        file.write("echo '** Done **'\n")
        file.write('echo "============================================="\n')
        file.write('################################################################################\n')
        file.write('################################################################################\n')
        file.write('# END #\n')
        file.write('################################################################################\n')
        file.write('################################################################################\n')
        file.write("You can write whatever you want in it.\n")
        file.write("This is just an example.\n")
        file.close()
    
   
if __name__ == "__main__":
    main()
