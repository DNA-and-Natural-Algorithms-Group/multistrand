## Run this file to re-create the figures from the MS-2.0 paper.
## Run this as ./overnight.sh > overnight-DATE.txt

## This is the sample trajectory, figure 2
## single process, <1s

python sample_trace_fig1.py > sample_trace_fig1_output.txt;

## These are the first passage figure
## 8 processes, < 15 minutes

python barplots.py 8 2400 > barplots_output.txt;

## These are the comparison between first-step and trajectory mode
## 8 processes < 1 hr

python case1.py plots 1600;
python case1.py slowDownStudy 1600;

## These are the trajectory studies of P0, P3, P4.
## 8 processes, <24 hr

python case3.py 8 1600 P4;		
python case3.py 8 1600 P3;
python case3.py 8 1600 P0;

## These are the iso-random studies 
## 8 processes, < 24 hrs


#python hybridization_F2.py generate iso-random 15 200 ir15-200;
#python hybridization_F2.py generate iso-random 25 50 ir25-50;				

python case2.py generate iso-random 15 120 ir15-120;
python case2.py generate iso-random 25 120 ir25-120;				


