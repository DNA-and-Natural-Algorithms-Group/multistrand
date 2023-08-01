# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

## Run this file to re-create the figures from the MS-2.0 paper.
## Run this as ./overnight.sh > overnight-DATE.txt

## This is the sample trajectory, figure 2
## single process, <1s

python sample_trace_fig1.py > sample_trace_fig1_output.txt;

## These are the first passage figure
## 8 processes, < 15 minutes

python barplots.py 8 2400 bonnet > barplots_output_bonnet.txt;
python barplots.py 8 2400 flamm  > barplots_output_flamm.txt;
python barplots.py 8 2400 yurke  > barplots_output_yurke.txt;
python barplots.py 8 2400 yurke2 > barplots_output_yurke2.txt;

## These are the comparison between first-step and trajectory mode
## 8 processes < 1 hr

python case1.py plots 3200			> case1_plots.txt;
python case1.py slowDownStudy 3200		> case1_slowdown.txt;

## These are the trajectory studies of P0, P3, P4.
## 8 processes, <1 hr

python case3.py 8 3200 P4			> case3_P4.txt	
python case3.py 8 3200 P3;			> case3_P3.txt
python case3.py 8 3200 P0;			> case3_p0.txt

## Rickettsia.
## 8 processes, < 4 hrs

python barplots.py 8 9600 rickettsia > barplots_output_rickettsia.txt;


## These are the iso-random studies 
## 8 processes, < 24 hrs

python case2.py generate iso-random 15 120 ir15-120 	> case2_iso15.txt;
python case2.py generate iso-random 25 120 ir25-120 	> case2_iso25.txt;				



