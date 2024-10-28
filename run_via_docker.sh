#!/bin/bash

# Ermittelt den Pfad des Shell-Skripts und gibt ihn aus
SCRIPT_DIR=$(dirname "$(realpath "$0")")
cd $SCRIPT_DIR

docker build -t automatedpertools/programm:1.0 .
docker run -ti --name AutomatedPERTools automatedpertools/programm:1.0

docker cp AutomatedPERTools:./usr/local/bin/server_run_collection .
docker cp AutomatedPERTools:./usr/local/bin/experiment.log .

docker rm AutomatedPERTools