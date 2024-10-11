@echo off

REM Get the path of the script and change directory to it
SET SCRIPT_DIR=%~dp0
cd /d %SCRIPT_DIR%

REM Build the Docker image
docker build -t automatedpertools/programm:1.0 .

REM Run the Docker container
docker run -ti --name AutomatedPERTools automatedpertools/programm:1.0

docker cp AutomatedPERTools:./usr/local/bin/data.pickle .

docker rm AutomatedPERTools

