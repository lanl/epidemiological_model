#!/bin/sh

LOGS_DIR="/Users/jkeithley/Documents/CIMMID/human/dengue_model/epi_seir/logs"
HUMAN_OUTPUT_DIR="/Users/jkeithley/Documents/CIMMID/human/dengue_model/epi_seir/human_model_output"

mv $LOGS_DIR/*.log $LOGS_DIR/logfile_archive
rm $HUMAN_OUTPUT_DIR/*.csv
rm $HUMAN_OUTPUT_DIR/*.parquet
