#!/bin/bash
# Fraser King - 2017

# extract.sh - This program runs in tandem with read_hdf5.py to
# extract relevant data from the CloudSat HDF4 files within
# a range around superstations of interest

# Variable setup
declare input_path=$1
declare output_path=$2
declare total=$(find $input_path -type f -name "*.hdf" | wc -l)
declare count=0

echo -e "Running CSA CloudSat Extraction Script.\n\nSee extraction_log.txt in this directory to follow the extraction progress.\n"

find $input_path -type f -name "*.hdf" | while read filename;
 do 
 	# Conversion runloop
    echo -e $(date -u) "\nConverting ${filename} to h5";
    name=$(basename "$filename" .hdf);
    ./h4toh5 $filename "${output_path}/${name}.h5";
    ((count++));

    echo "Extracting data via Python...";
    h5path="${output_path}/${name}.h5";
    h5name="${name}.h5";
    python3 extract_ecmwf_temps.py $h5path $h5name $output_path;
    echo "Extraction complete";

    echo -e "Deleting temporary .h5 file..."
    rm $h5path;

    echo -e "File complete! Progress: ${count}      /${total}\n";
 done > extraction_log.txt

 echo -e "\nExtraction complete!"