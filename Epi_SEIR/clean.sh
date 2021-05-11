#!/bin/sh

if [ "$1" == "" ]; then
	echo "Error: Must specify configuration file path"
	exit 1
fi

CONFIG_FILE_PATH=$1
LOGS_DIR=$(grep "LOGFILE_PATH" $CONFIG_FILE_PATH | tr -s " " | cut -d " " -f 2 | cut -d "'" -f 2)
OUTPUT_DIR=$(grep "OUTPUT_DIR" $CONFIG_FILE_PATH | tr -s " " | cut -d " " -f 2 | cut -d "'" -f 2)

mv ${LOGS_DIR%%/}/*.log ${LOGS_DIR%%/}/logfile_archive
rm ${OUTPUT_DIR%%/}/*.csv
rm ${OUTPUT_DIR%%/}/*.parquet
