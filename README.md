# TGT_Virus_Detection_Simulation

This repository contains code for simulating virus detection using ONT sequencing and BLAST. In addition, a web application is included which describes the application of the hypergeometric probability distribution to the problem of virus detection using ONT sequencing and BLAST. 


-----------------------------------------------------------------------------------------
Simulating virus detection with ONT sequencing and BLAST consists of the following steps:
1. Select a virus genome and a host genome. 
2. Choose a ratio of virus and host genomes. 
3. Fragment copies of the virus and host genomes and add errors to the fragments to generate reads.
4. Choose a reference database to BLAST the simulated reads against. 
5. Run BLAST. 

The following scripts implement the steps described above: 
1. virus_read_detection_simulation.py
2. loop_virus_read_detection_simulation.py
3. virus_read_detection_blastn_search.py

virus_read_detection_simulation.py is the main script. It takes arguments that describe the parameters of the simulation including the selected virus and host genomes, the quantity of each, the error rate, and the reference database. The script generates the reads and then calls virus_read_detection_blastn_search.py in a loop to BLAST the reads. 

loop_virus_read_detection_simulation.py is an alternate main script that extends virus_read_detection_simulation.py by taking an additional argument specifying the number of times to repeat the simulation. This script is only for use in simulations where a subset of the reads are BLASTed. For example, simulations can be run to obtain a distribution for the number of reads sequenced before obtaining the first virus sequence. 

------------------------------------------------------------------------------------------




