#!/bin/bash

# Ermittelt den Pfad des Shell-Skripts und gibt ihn aus
SCRIPT_DIR=$(dirname "$(realpath "$0")")
cd $SCRIPT_DIR

docker build -t automatedpertools/programm:1.0 .
docker run --rm -ti --name AutomatedPERTools automatedpertools/programm:1.0