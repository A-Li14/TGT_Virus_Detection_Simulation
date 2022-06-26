# TGT_Virus_Detection_Simulation

This repository contains code for simulating virus detection using ONT sequencing and BLAST. In addition, two web applications are included which describe the results of simulated ONT sequencing and a hypergeometric model of virus detection using ONT sequencing and BLAST. 

-----------------------------------------------------------------------------------------
Simulating virus detection with ONT sequencing and BLAST consists of the following steps:
1. Select a virus genome and a host genome. 
2. Choose a ratio of virus and host genomes. 
3. Fragment copies of the virus and host genomes and add errors to the fragments to generate reads.
4. Choose a reference database to BLAST the simulated reads against. 
5. Run BLAST. 

The following scripts in the "Simulation" folder implement the steps described above: 
1. virus_read_detection_simulation.py
2. loop_virus_read_detection_simulation.py
3. virus_read_detection_blastn_search.py

virus_read_detection_simulation.py is the main script. It takes arguments that describe the parameters of the simulation including the selected virus and host genomes, the quantity of each, the error rate, and the reference database. The script generates the reads and then calls a function described in virus_read_detection_blastn_search.py in a loop to BLAST the reads. 

loop_virus_read_detection_simulation.py is an alternate main script that extends virus_read_detection_simulation.py by taking an additional argument specifying the number of times to repeat the simulation. This script is only for use in simulations where a subset of the reads are BLASTed. For example, simulations can be run to obtain a distribution for the number of reads sequenced before obtaining the first virus sequence. 

------------------------------------------------------------------------------------------

The following scripts in the "Dash Apps" folder implement web applications that summarize different steps in the simulation of virus detection.
1. dashApp_fragmentation.py
2. Hypergeometric app.py

dashApp_fragmentation.py simulates ONT sequencing of a chosen reference genome and summarizes the results in terms of the number of reads, the range of read lengths, and a histogram of read lengths. 

Hypergeometric app.py models ONT sequencing using the hypergeometric probability distribution. Given a prespecified number of virus and background reads, the app plots the likelihood of sequencing a virus read as a function of the total number of reads sequenced. In addition, the app calculates the number of reads that need to be sequenced to achieve threshold levels of sensitivity. 


