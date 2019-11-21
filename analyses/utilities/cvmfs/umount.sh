#!/usr/bin/env bash

echo "Requesting sudo to unmount CVMFS repos!"

# Mount main alice repo
sudo umount /Users/Shared/cvmfs/alice.cern.ch
# OCDB
sudo umount /Users/Shared/cvmfs/alice-ocdb.cern.ch
