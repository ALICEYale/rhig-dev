#!/bin/bash

echo ===========================
BASEDIR=$PWD
BASEFOLDER="ManualMerging"
BASE=$BASEDIR/$BASEFOLDER
PASS=pass1
STAGE=1
MERGESTAGE=Merge${STAGE}
#—————————————————
YEAR=2015
PERIOD="LHC15o"
#TRAIN=2195
#TIME=_20170615-2337
#TRAIN=2267
#TIME=_20170713-1715
TRAIN=2675
TIME=_20171031-2234
#—————————————————

SUBDIR="Train_$TRAIN"
cd $BASE/$SUBDIR
mkdir $BASE/$SUBDIR/$MERGESTAGE
cd $BASE/$SUBDIR/
#cd $BASE/$SUBDIR/Merge1

#- - - - - - - Loop - - - - - - - - - - - -
#####LHC15o#####
#FolderNumber =000246809
#FolderNumber =000245952
#FolderNumber =000246087
FolderNumber=000245952

i=4
while [  $i -lt 16 ]; do

	echo The counter is $i
	inFILE1="${FolderNumber}_AnalysisResults_${i}.root"
	let i=i+$STAGE
	inFILE2="${FolderNumber}_AnalysisResults_${i}.root"
	let i=i+$STAGE 
	inFILE3="${FolderNumber}_AnalysisResults_${i}.root"
	#let i=i+$STAGE 
	#inFILE4="${FolderNumber}_AnalysisResults_${i}.root"
	
	echo Add files $inFILE1 , $inFILE2 , $inFILE3 ,  #$inFILE4        
	hadd $BASE/$SUBDIR/$MERGESTAGE/${FolderNumber}_AnalysisResults_${i}.root $inFILE1 $inFILE2 $inFILE3 #$inFILE4

	let i=i+$STAGE 
done
#- - - - - - - Loop - - - - - - - - - - - -

cd $BASEDIR
echo ===========================

