
Frits Dannenberg May 17th 2017  

This folder contains the file required to generate the plots and tables of the hybridization case study in the Multistrand 2.0 journal version  

barplots: The figure 3 barplots of first passage times  
Case 1: Table of Gao 2006 and the concentration-dependent rates and First Step mode comparison  
Case 2: Scatterplots of association and dissociation rates for randomized and Gao 2006 sequences  
Case 3: Heathmap for hybridization states  
sample_trace_fig1 : The sample trace in figure 2 



To test these files, try 

barplots 2 30
case1.py plots 160  
case1.py slowDownStudy 1200  
case2.py generate iso-random 15 10 ir15-10  
case3.py 2 24 test2  
sample_trace_fig1

Warning! These scripts will create several compute threads  
