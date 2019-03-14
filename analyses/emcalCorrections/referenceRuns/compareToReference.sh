#!/usr/bin/env bash

# New Filename is defined from the emcalCorrections directory
newFilename=${1:-AnalysisResults.root}
# Reference filename is defined from the referenceRuns directory
referenceFilename=${2:-AnalysisResults.root}

referenceRunsPath="../emcalCorrections/referenceRuns"

python prepareReferenceForComparison.py -f "${referenceFilename}"
echo -e "INFO: Finished reference file prep!"

# Move to the utilities directory because it is easiest to the run the comparison script when localized from there
cd ../../utilities
echo -e "\nINFO: Starting comparison!"
python compareHistos.py -f "${referenceRunsPath}/../${newFilename}" "${referenceRunsPath}/${referenceFilename}" -y "${referenceRunsPath}/referenceMap.yaml" --overrideCorrectionsMap
