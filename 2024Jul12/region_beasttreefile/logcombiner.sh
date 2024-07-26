#!/bin/bash

# Default values
burnin=25
prefix="input"
tree=false
resample=10000

# Function to display usage
function usage {
    echo "Usage: $0 -p prefix [-b burnin] [-t] [-r resample]"
    exit 1
}

# Parse command-line options
while getopts ":p:b:tr:" opt; do
  case ${opt} in
    p )
      prefix=$OPTARG
      ;;
    b )
      burnin=$OPTARG
      ;;
    t )
      tree=true
      ;;
    r )
      resample=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      usage
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      usage
      ;;
  esac
done
shift $((OPTIND -1))

# Check if prefix is provided
if [ -z "$prefix" ]; then
    echo "Prefix is required"
    usage
fi

# Get list of input files with the given prefix
if [ "$tree" = true ]; then
    input_files=("${prefix}"*.trees)
else
    input_files=("${prefix}"*.log)
fi

# Check if any input files were found
if [ ${#input_files[@]} -eq 0 ]; then
    echo "Error: No input files found with prefix '$prefix'"
    exit 1
fi

# Output filename
if [ "$tree" = true ]; then
    output_file="${prefix}_combined.trees"
else
    output_file="${prefix}_combined.log"
fi

# Run logcombiner
if [ "$tree" = true ]; then
    logcombiner -burnin "$burnin" -trees -resample "$resample" 
"${input_files[@]}" "$output_file"
else
    logcombiner -burnin "$burnin" -resample "$resample" 
"${input_files[@]}" "$output_file"
fi

# Optionally, check if the output file exists and inform the user
if [ -f "$output_file" ]; then
    echo "Output file '$output_file' generated successfully."
else
    echo "Error: Output file '$output_file' was not generated."
fi
