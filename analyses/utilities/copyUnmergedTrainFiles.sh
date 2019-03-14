#!/bin/bash

echo ===========================
BASEDIR=$PWD
BASEFOLDER="ManualMerging"
BASE=$BASEDIR/$BASEFOLDER
PASS=pass1

#—————————————————
YEAR=2015
PERIOD="LHC15o"
#TRAIN=2195
#TIME=_20170615-2337
#TRAIN=2267
#TIME=_20170713-1715
#TRAIN=2272
#TIME=_20170721-1802 
TRAIN=2675
TIME=_20171031-2234
#—————————————————

SUBDIR="Train_$TRAIN"
mkdir $BASEFOLDER
mkdir $BASEFOLDER/$SUBDIR
cd $BASE/$SUBDIR

#- - - - - - - Loop - - - - - - - - - - - -
#####LHC15n#####
#FolderNumber=000245952
#FolderNumber=000246113
#FolderNumber=000246809
#FolderNumber=000246087
FolderNumber=000245952

i=10
while [  $i -lt 16 ]; do

        ##The file name must be /00$i/ for i<10 and /0$i/ for i<100 and /$i/ for everything else 
	inFILE="alien:///alice/data/$YEAR/$PERIOD/$FolderNumber/$PASS/PWGJE/Jets_EMC_PbPb/$TRAIN$TIME/Stage_2/0$i/AnalysisResults.root"
        #inFILE="alien:///alice/data/$YEAR/$PERIOD/$FolderNumber/pass1/PWGJE/Jets_EMC_PbPb/$TRAIN$TIME/Stage_1/00${i}/AnalysisResults.root"
	echo  - - - - - - - - - - - - - - 
	echo  start copying File $i
	echo  $inFILE

	outputFILE="${FolderNumber}_AnalysisResults_${i}.root"
	echo Copy to local file: ${outputFILE}

	alien_cp $inFILE $outputFILE

	let i=i+1 
done
#- - - - - - - Loop - - - - - - - - - - - -

cd $BASEDIR
echo ===========================

