#!/bin/bash


# DOWNLOADING GENE EXPRESSION DATA
## installing gdc-client tool
wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip
unzip gdc-client_v1.6.1_Ubuntu_x64.zip
sudo mv gdc-client /usr/bin && rm gdc-client_v1.6.1_Ubuntu_x64.zip
gdc-client --help

##downloadi files using gdc-client tool
mkdir tnbc
gdc-client download -m manifest.txt -d tnbc/


# DRUGS CONVERSION AND ENERGY MINIMIZATION
## split the multi compound sdf files using open babel
obabel -isdf fda.sdf -osdf -O *.sdf --split
obabel -isdf repurposed.sdf -osdf -O *.sdf --split

## enery minimization using open babel
for file in *.sdf;
do obminimize -o hin -sd -ff MMFF94 -h "$file" > "${file%.sdf}-minimized.hin";
done