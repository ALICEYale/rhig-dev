#!/usr/bin/env bash

echo "Requesting sudo to unmount CVMFS repos!"

# Mount main alice repo
sudo umount /cvmfs/alice.cern.ch
# OCDB
sudo umount /cvmfs/alice-ocdb.cern.ch
