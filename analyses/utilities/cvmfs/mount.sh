#!/usr/bin/env bash

echo "Requesting sudo to mount CVMFS repos!"

# Mount main alice repo
sudo mount -t cvmfs alice.cern.ch /Users/Shared/cvmfs/alice.cern.ch
# OCDB
sudo mount -t cvmfs alice-ocdb.cern.ch /Users/Shared/cvmfs/alice-ocdb.cern.ch
