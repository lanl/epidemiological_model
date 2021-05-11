#!/bin/sh

LOGS_DIR="logs"
HUMAN_OUTPUT_DIR="human_model_output"

mv $LOGS_DIR/*.log $LOGS_DIR/logfile_archive
rm $HUMAN_OUTPUT_DIR/*.csv
rm $HUMAN_OUTPUT_DIR/*.parquet
