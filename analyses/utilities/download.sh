#!/bin/bash

# Script to download files from the grid, and write them into a similar directory structure on your machine.
# Edit the config fields below to set the appropriate path of the files you want to download, and the destination path.
# It will also automatically write the file list if desired
# Author: Salvatore
# Edits by: Raymond

# Select system and period
# Possible systems are "PbPb", "pPb", and "pp"
system="PbPb"
period="LHC15o"
# Select options:
# Max number of files to have locally
maxFiles=2
# Max number of files to download
maxFilesDownloaded=10
# Where to write the files
# Remember to include the trailing slash!
destPrefix="$HOME/code/alice/data/${period}/"
# Name of the folder containing the file list
# It will go in the destPrefix/../$fileListFolderName folder
# An empty filename will prevent the file list from being generated
fileListFolderName="fileLists"

# Other defaults
# They override the values set below
# data can be either "sim" or "data"
data=""
year=""
run=""
pass=""
# NOTE: The filetype may need to be followed by the production number, depending on the data set
fileType=""
productionNumber=""
# pT Hard Bins (when applicable)
ptHardBin="4"
# Filename is often "root_archive.zip" or "aod_archive.zip", depending on the data set
filename=""

# Below are known systems, which should work out of the box.
# LHC12a15a: Jet-jet pp 2.76 TeV MC in pt hard bins anchored to LHC11a.
# WARNING: Does _NOT_ contain EMCal clusters!!
if [[ "$system" == "pp" && "$period" == "LHC12a15a" ]];
then
    data=${data:-"sim"}
    year=${year:-"2012"}
    run=${run:-"146805"}
    if [[ -z "${ptHardBin}" ]]; then
        echo "ERROR: pt hard bin is not selected! Please select one and try again!"
        exit 1
    fi
    fileType=${fileType:-"AOD"}
    filename=${filename:-"root_archive.zip"}

    searchpath="/alice/${data}/${year}/${period}/${run}/${ptHardBin}"
fi
# LHC15g1a: Jet-jet pp 2.76 TeV MC in pt hard bins anchored to LHC11a.
# This production _does_ contain EMCal clusters!
# https://alice.its.cern.ch/jira/browse/ALIROOT-6125
# https://twiki.cern.ch/twiki/bin/view/ALICE/JetMCProductionsCrossSections#LHC15g1a_Jet_Jet_events_pp_colli
if [[ "$system" == "pp" && "$period" == "LHC15g1a" ]];
then
    data=${data:-"sim"}
    year=${year:-"2015"}
    run=${run:-"146805"}
    if [[ -z "${ptHardBin}" ]]; then
        echo "ERROR: pt hard bin is not selected! Please select one and try again!"
        exit 1
    fi
    fileType=${fileType:-"AOD"}
    filename=${filename:-"aod_archive.zip"}

    searchpath="/alice/${data}/${year}/${period}/${ptHardBin}/${run}/${fileType}"
fi
# LHC12a15e_fix: Jet-jet pp 2.76 TeV MC in pt hard bins anchored to LHC11h.
if [[ "$system" == "pp" && "$period" == "LHC12a15e_fix" ]];
then
    data=${data:-"sim"}
    year=${year:-"2012"}
    run=${run:-"169838"}
    if [[ -z "${ptHardBin}" ]]; then
        echo "ERROR: pt hard bin is not selected! Please select one and try again!"
        exit 1
    fi
    fileType=${fileType:-"AOD"}
    productionNumber=${productionNumber:-"149"}
    filename=${filename:-"root_archive.zip"}

    searchpath="/alice/${data}/${year}/${period}/${run}/${ptHardBin}/${fileType}${productionNumber}"
fi
# LHC16e1: pp 5.02 TeV MC anchored to LHC15n
# https://alice.its.cern.ch/jira/browse/ALIROOT-6671
if [[ "$system" == "pp" && "$period" == "LHC16e1" ]];
then
    data=${data:-"sim"}
    year=${year:-"2016"}
    run=${run:-"244351"}
    if [[ -z "${ptHardBin}" ]]; then
        echo "ERROR: pt hard bin is not selected! Please select one and try again!"
        exit 1
    fi
    fileType=${fileType:-"AOD"}
    filename=${filename:-"root_archive.zip"}

    searchpath="/alice/${data}/${year}/${period}/${run}/${ptHardBin}/"
    if [[ "${fileType}" == "AOD" ]];
    then
        # Add "AOD" to searchpath and change filename
        searchpath="${searchpath}${fileType}/"
        filename="aod_archive.zip"
    fi
fi
# LHC16j5: pp 5.02 TeV jet-jet PYTHIA anchored to LHC15o
# https://alice.its.cern.ch/jira/browse/ALIROOT-6905
if [[ "$system" == "pp" && "$period" == "LHC16j5" ]];
then
    data=${data:-"sim"}
    year=${year:-"2016"}
    run=${run:-"246945"}
    if [[ -z "${ptHardBin}" ]]; then
        echo "ERROR: pt hard bin is not selected! Please select one and try again!"
        exit 1
    fi
    fileType=${fileType:-"AOD"}
    filename=${filename:-"root_archive.zip"}
    productionNumber=${productionNumber:-"200"}

    # AOD200: /alice/sim/2016/LHC16j5/20/246945/AOD200/0004/AliAOD.root
    searchpath="/alice/${data}/${year}/${period}/${ptHardBin}/${run}/"
    if [[ "${fileType}" == "AOD" ]];
    then
        # Add "AOD" to searchpath and change filename
        searchpath="${searchpath}${fileType}${productionNumber}/"
        filename="AliAOD.root"
    fi
fi
# LHC15i2c: pp 7 TeV D mseon MC
if [[ "$system" == "pp" && "$period" == "LHC15i2c" ]];
then
    data=${data:-"sim"}
    year=${year:-"2015"}
    run=${run:-"120822"}
    pass=${pass:-"4"}
    fileType=${fileType:-"AOD"}
    filename="aod_archive.zip"

    searchpath="/alice/${data}/${year}/${period}/${run}/${pass}/${fileType}"
fi
# pPb 2013 data
if [[ "$system" == "pPb" && "$period" == *"LHC13"* ]];
then
    data=${data:-"data"}
    year=${year:-"2013"}
    # Corresponds to a run in LHC13b
    run=${run:-"000195344"}
    pass=${pass:-"4"}
    fileType=${fileType:-"AOD"}
    filename=${filename:-"aod_archive.zip"}

    if [[ "${period}" == "LHC13d" ]];
    then
        run="000195724"
    fi

    searchpath="/alice/${data}/${year}/${period}/${run}/pass${pass}/${fileType}"
fi
# LHC11a: pp 2.76 TeV data
if [[ "$system" == "pp" && "$period" == "LHC11a" ]];
then
    data=${data:-"data"}
    year=${year:-"2011"}
    run=${run:-"000146858"}
    pass=${pass:-"4"}
    fileType=${fileType:-"AOD"}
    filename=${filename:-"AliAOD.root"}
    productionNumber=${productionNumber:-"113"}

    # /alice/data/2011/LHC11a/000146860/ESDs/pass4_with_SDD/AOD113/0002
    searchpath="/alice/${data}/${year}/${period}/${run}/ESDs/pass${pass}_with_SDD/${fileType}${productionNumber}"
fi
# LHC11h: 2011 PbPb 2.76 data
if [[ "$system" == "PbPb" && "$period" == "LHC11h" ]];
then
    data=${data:-"data"}
    year=${year:-"2011"}
    run=${run:-"000168464"}
    pass=${pass:-"2"}
    fileType=${fileType:-"AOD"}
    filename=${filename:-"root_archive.zip"}
    productionNumber=${productionNumber:-"145"}

    # /alice/data/2011/LHC11h_2/000167693/ESDs/pass2/AOD145
    searchpath="/alice/${data}/${year}/${period}/${run}/ESDs/pass${pass}/${fileType}${productionNumber}"
fi
# LHC15o: 2015 PbPb 5.02 data
# Defaults to using AOD194, but can modify back to the original pass1 AOD by unsetting the production number
# and specifying the filename as "aod_archive.zip".
if [[ "$system" == "PbPb" && "$period" == "LHC15o" ]];
then
    data=${data:-"data"}
    year=${year:-"2015"}
    # The first run from the Calo run list
    run=${run:-"000246928"}
    pass=${pass:-"1"}
    fileType=${fileType:-"AOD"}
    filename=${filename:-"root_archive.zip"}
    productionNumber=${productionNumber:-"194"}

    # AOD:
    # /alice/data/2015/LHC15o/000246928/pass1/AOD/001
    # AOD194:
    # /alice/data/2015/LHC15o/000245683/pass1/AOD194
    searchpath="/alice/${data}/${year}/${period}/${run}/pass${pass}/${fileType}${productionNumber}"
fi
# 2017 5 TeV pp reference
if [[ "$system" == "pp" && "$period" == "LHC17p" || "$period" == "LHC17q" ]];
then
    data=${data:-"data"}
    year=${year:-"2017"}
    # The first run from the Central Barrel Tracking list
    run=${run:-"000282343"}
    pass=${pass:-"1_FAST"}
    fileType=${fileType:-"AOD"}
    filename=${filename:-"AliAOD.root"}
    productionNumber=${productionNumber:-"208"}

    # AOD208:
    # /alice/data/2017/LHC17p/000282343/pass1_FAST/AOD208/0001/
    searchpath="/alice/${data}/${year}/${period}/${run}/pass${pass}/${fileType}${productionNumber}"
fi

# Fall back to default if it is an unknown period
if [[ -z "$searchpath" ]];
then
    echo "Unkonwn period and system selected. Using default path and values."
    searchpath="/alice/${data}/${year}/${period}/${run}/${pass}/${fileType}"
fi

echo "Looking for file ${filename} in path ${searchpath}"

read -p "Search the grid? Press enter to continue or ctrl-c to cancel!"

files=$(alien_find  ${searchpath} ${filename})

echo "Files found on the grid:"
echo "${files}"

echo "Downloading at most $maxFilesDownloaded files to have at most $maxFiles on the system!"
read -p "Ready to download the files? Press enter to continue or ctrl-c to cancel!"

nFiles=0
nFilesDownlaoded=0

for file in $files
do
    if [[ "${nFilesDownlaoded}" -ge "${maxFilesDownloaded}" ]]; then
        break
    fi
    if [[ "${nFiles}" -ge "${maxFiles}" ]]; then
        break
    fi

    if [[ ! $file == /alice* ]]; then
        continue
    fi
    localpath=${file/\/alice\/${data}\/${year}\/${period}\//}
    localpath=${localpath/\/${filename}/}

    dest="${destPrefix}${localpath}"
    fileDest="${dest}/${filename}"
    if [ ! -f ${fileDest} ]; then
        mkdir -p ${dest}
        echo "Copying alien://${file} to ${fileDest}"
        alien_cp alien://${file} ${fileDest}
        ((nFilesDownlaoded++))
    else
        echo "File already exists: ${fileDest}"
    fi
    ((nFiles++))
done

# Generate file list
if [[ -n "$fileListFolderName" ]];
then
    fileListDest="${destPrefix/\/${period}\//}/${fileListFolderName}"
    echo "Generating filelist ${period}.txt in ${fileListDest}"
    fileListFilename="${fileListDest}/${period}${fileType}.txt"
    # "%/" removes trailing slash
    find ${destPrefix%/} -name ${filename} > ${fileListFilename}
    # Add the proper extensions if necessary
    # ie ../root_archive.zip#AliAOD.root
    fileExtension="${filename##*.}"
    if [[ "${fileExtension}" == "zip" ]];
    then
        if [[ "${fileType}" == "AOD" ]];
        then
            addAfterHash="AOD"
        else
            addAfterHash="ESDs"
        fi

        # Append to each line
        awk -v addAfterHash="${addAfterHash}" '{print $0 "#Ali" addAfterHash ".root"}' ${fileListFilename} > ".tmp" && mv ".tmp" ${fileListFilename}
        echo "Created file list at \"${fileListFilename}\""
        echo "NOTE: The filename has had the file type (AOD or ESD) appended to it!"
    fi
else
    echo "Not generating file list!"
fi
