#!/bin/bash -e

# This script *MUST* be sourced, otherwise export PATH does
# not work for future shells. To source, use:

#	source SCRIPTNAME
# or
#	. SCRIPTNAME

# There seems to be no other way to update the PATH for future shell 
# invocations, if I missed something let me know robert@drive5.com.

sudo apt update
sudo apt install -y awscli

git clone https://github.com/rcedgar/palm_annot.git
chmod +x palm_annot/bin/* palm_annot/py/* palm_annot/bash/*
sudo cp -v palm_annot/bin/* /usr/bin

export PATH=$PATH:~/palm_annot/py
