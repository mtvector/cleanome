#!/bin/bash
#download mitofinder sif (had to do directly from web)
#https://cloud.sylabs.io/library/remiallio/default/mitofinder
#singularity pull library://remiallio/default/mitofinder

singularity run ../remiallio_default_mitofinder.sif -j rhemac10 -a ./NC_005943.1.fa -r ./homo_sapiens.gb -o 2 --blast-identity-prot 30 --blast-size 20 --blast-identity-nucl 30 

