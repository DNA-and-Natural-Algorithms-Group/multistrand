## Run this file to re-create the figures from the MS-2.0 paper.
## Run this as ./overnight.sh > overnight-DATE.txt

## This is the sample trajectory, figure 2
## single process, <1s

python sample_trace_fig1.py > sample_trace_fig1_output.txt;


## These are the first passage figure
## 8 processes, < 1 hr

python barplots 8 800 > barplots_output.txt;


## These are the comparison between first-step and trajectory mode
## 8 processes < 1 hr


python hybridization_F1.py plots 1600;
python hybridization_F1.py slowDownStudy 1600;

## These are the trajectory studies of P0, P3, P4.
## 8 processes, <1 hr

python hybridization_F3.py 8 1600 P4;
python hybridization_F3.py 8 1600 P3;
python hybridization_F3.py 8 1600 P0;

## These are the iso-random studies 
## 8 processes, < 24 hrs


#python hybridization_F2.py generate iso-random 15 200 ir15-200;
#python hybridization_F2.py generate iso-random 25 50 ir25-50;				

python hybridization_F2.py generate iso-random 15 500 ir15-500;
python hybridization_F2.py generate iso-random 25 500 ir25-500;				#<-- we got here.


