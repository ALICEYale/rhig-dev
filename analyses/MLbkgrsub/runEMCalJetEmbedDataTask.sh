#!/bin/bash

if [ -z "$1" ]; then
	DATATYPE="AOD"
else
	DATATYPE="${1}"
fi

if [ -z "$2" ]; then
	PERIOD="LHC15o"
else
	PERIOD="${2}"
fi

if [ -z "$3" ]; then
    FILELIST="files_test.txt"
else
	FILELIST="${3}"
fi

if [ -z "$4" ]; then
	NEVENTS="5000"
else
	NEVENTS="${4}"
fi

root -b -q -x runEMCalJetEmbedDataTask.C\(\""${DATATYPE}"\",\""${PERIOD}"\",\""${FILELIST}"\","${NEVENTS}"\)
