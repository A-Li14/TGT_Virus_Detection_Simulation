###This file analyzes the results of looped virus detection simulations


import glob
import pandas as pd
import re


def get_runtime(results_txt_file) :
    '''Return the runtime of a simulation. Also include the number of lines'''
    with open(results_txt_file) as file :
        line = file.readline()
        count = 1
        nextline = file.readline()
        while nextline != "" :
            line = nextline
            count = count+1
            nextline = file.readline()
    return [line[16:(line.find("min")-1)],count]



####Figure 2###########################################

output_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/Analyzing_Results/Runtime.csv"

###Find all txt files in a results directory
results_txt = glob.glob("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/*/*.txt",recursive=True)

sim_data = pd.DataFrame()


###Extract data from each text file; each text file summarizes a simulation
runtime = []
lines = []
sim_num = []
for path in results_txt :
    temp = get_runtime(path)
    runtime.append(temp[0])
    lines.append(temp[1])
    ind = path.find("Simulation")
    ind = path.find("Simulation",ind+1)
    sim_num.append(path[ind:path.find("/",ind)])


###Store simulation data together
sim_data ["sim_num"] = sim_num
sim_data["runtime"] = runtime
sim_data["lines"] = lines


###Exclude simulation data from a prior experiment
###Different versions of simulations were labeled "Simulation1", "Simulation2", etc. 
###Similar simulations were labeled "Simulation1.1", "Simulation1.2", etc. 

# sim_data.drop(sim_num)
drop_ind = []
drop_ind = [i for (i,j) in enumerate(sim_data["sim_num"]) if j.find("Simulation1.")!=-1]
drop_ind = drop_ind + [i for (i,j) in enumerate(sim_data["sim_num"]) if j == "Simulation1"]
sim_data = sim_data.drop(drop_ind)
sim_data["runtime_num"] = sim_data["runtime"].astype("float")
sim_data.sort_values("lines")
sim_data.sort_values("runtime_num")

sim_data.drop("runtime",1)



###Select and analyze data from Simulation4
inds = [i for i in range(len(sim_data)) if sim_data.iloc[i]["sim_num"].find("Simulation4")!=-1]

runtime_length = sim_data.iloc[inds]
runtime_length = runtime_length.sort_values("sim_num")

runtime_length["Read_Length"] = ["200-300","500-600","1000-1100"]

plot = runtime_length.plot.scatter(x="Read_Length",y="runtime_num",xlabel="Read Lengths",ylabel="Minutes",title="Read Mapping Duration for 5000 MVM Genomes",legend=None)
plot.get_figure().savefig("Runtime by Average Read Length.png")


###Add total number of reads to the data set from result csv files
results_csv = glob.glob("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/*/*result.csv",recursive=True)

###Load simulation results
total_reads = []
sim_num = []
for path in results_csv :
    temp = pd.read_csv(path)
    read_ids = temp.query_id.str.strip("read").astype("int")
    total_reads.append(max(read_ids))
    ind = path.find("Simulation")
    ind = path.find("Simulation",ind+1)
    sim_num.append(path[ind:path.find("/",ind)])

sim_data2 = pd.DataFrame()
sim_data2["sim_num"] = sim_num
sim_data2["Total_Reads"] = total_reads


###Exclude data from Simulation1
drop_ind = []
drop_ind = [i for (i,j) in enumerate(sim_data2["sim_num"]) if j.find("Simulation1.")!=-1]
drop_ind = drop_ind + [i for (i,j) in enumerate(sim_data2["sim_num"]) if j == "Simulation1"]
sim_data2 = sim_data2.drop(drop_ind)

sim_data.sim_num==sim_data.sim_num

sim_data["Total_Reads"] = sim_data2.Total_Reads



###Output data file
sim_data.drop("runtime",1).to_csv(output_path,index=None)




###Table 2; different virus to host ratios#######################
tab2_path = glob.glob("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/*6.[45]/*/*.txt",recursive=True)

tab2_data = pd.DataFrame()


###Load data from text result files
runtime = []
lines = []
sim_num = []
rep_num = []
for path in tab2_path :
    temp = get_runtime(path)
    runtime.append(temp[0])
    lines.append(temp[1])
    ind = path.find("Simulation")
    ind = path.find("Simulation",ind+1)
    sim_num.append(path[ind:path.find("/",ind)])
    ind = path.find("Rep")
    rep_num.append(path[ind:path.find("/",ind)])


tab2_data ["sim_num"] = sim_num
tab2_data["rep_num"] = rep_num
tab2_data["runtime"] = runtime
tab2_data["lines"] = lines
len(tab2_data)


###Load data from csv result files
results_csv = glob.glob("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/*6.[45]/*/*result.csv",recursive=True)

total_reads = []
sim_num = []
rep_num = []
for path in results_csv :

    temp = pd.read_csv(path)
    read_ids = temp.query_id.str.strip("read").astype("int")
    total_reads.append(max(read_ids))
    ind = path.find("Simulation")
    ind = path.find("Simulation",ind+1)
    sim_num.append(path[ind:path.find("/",ind)])
    ind = path.find("Rep")
    rep_num.append(path[ind:path.find("/",ind)])

sim_tab2 = pd.DataFrame()
sim_tab2["sim_num"] = sim_num
sim_tab2["rep_num"] = rep_num
sim_tab2["Total_Reads"] = total_reads

len(sim_tab2)
sum(sim_tab2.sim_num==tab2_data.sim_num)==len(sim_tab2)
sum(sim_tab2.rep_num==tab2_data.rep_num)==len(sim_tab2)

tab2_data['Total_Reads'] = total_reads


###Output data
tab2_data.to_csv("Simulation runtime med + low concentrations.csv")






####Figure 3###########################################
###1000 repeated simulations to map against the hypergeometric distribution predictions
results_txt6_3 = glob.glob("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES10%ReadErrorRate_10power2VirusHostRatio_Simulation6.3/*/*.txt",recursive=True)
len(results_txt6_3)

sim_data6_3 = pd.DataFrame()


###Load data
runtime6_3 = []
lines6_3 = []
sim_num6_3 = []
rep_num6_3 = []
for path in results_txt6_3 :
    temp = get_runtime(path)
    runtime6_3.append(temp[0])
    lines6_3.append(temp[1])
    ind = path.find("Simulation")
    ind = path.find("Simulation",ind+1)
    sim_num6_3.append(path[ind:path.find("/",ind)])
    ind = path.find("Rep")
    rep_num6_3.append(path[ind:path.find("/",ind)])


sim_data6_3["sim_num"] = sim_num6_3
sim_data6_3["runtime"] = runtime6_3
sim_data6_3["lines"] = lines6_3
sim_data6_3["Rep"] = rep_num6_3

sim_data6_3

sim_data6_3["runtime"] = sim_data6_3.runtime.astype("float")
sim_data6_3.runtime.mean()


###Load data
###Get total reads aligned in each simulation; since mapping stopped when a virus was detected, the highest query id will be the
###number of reads
paths = glob.glob("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES10%ReadErrorRate_10power2VirusHostRatio_Simulation6.3/*/*RVDB_result.csv",recursive=True)
total_reads = []
rep_num = []
avg_read_length = []
for path in paths :
    temp = pd.read_csv(path)
    read_ids = temp.query_id.str.strip("read").astype("int")
    total_reads.append(max(read_ids))
    ind = path.find("Rep")
    rep_num.append(path[ind:path.find("/",ind)])
    avg_read_length.append(temp.drop_duplicates("query_id").query_length.mean())



read_id_df = pd.DataFrame()
read_id_df["Rep"] = rep_num
read_id_df["Reads_Mapped"] = total_reads
read_id_df["avg_read_length"] = avg_read_length

sum(read_id_df.Rep==sim_data6_3.Rep) == 1000

sim_data6_3["Reads_Mapped"] = total_reads
sim_data6_3["avg_read_length"] = avg_read_length

sim_data6_3.plot.scatter(x="Reads_Mapped",y="runtime")

# sim_data6_3.to_csv("Sim6_Hypergeometric Comparison/Figure3_data.csv")
sim_data6_3.to_csv("Sim6_Hypergeometric Comparison/Simulation_1000x_runtime.csv")




###Compare simulation results to hypergeometric simulation
from scipy.stats import hypergeom
import math

###Calculate approximate parameters to the hypergeometric distribution based on simulation data
host_bases = 2412249049 ###CHO
viral_bases = 5149*1538 ###MVM
average_read_length = (500+5000)/2

host_reads = math.ceil(host_bases/average_read_length)
viral_reads = math.ceil(viral_bases/average_read_length)


###Calculate likelihood of detecting virus deads
p_detect = 1-hypergeom.pmf(M=host_reads+viral_reads,n=viral_reads,N=range(2000),k=0)



###Plot likelihood of detecting virus reads from the hypergeometric distribution and simulation results
plt.plot(p_detect)

sim_data6_3.Reads_Mapped.max()
sim_p_detect = [(sim_data6_3.Reads_Mapped<i).sum()/1000 for i in range(sim_data6_3.Reads_Mapped.max())]


plt.plot(sim_p_detect)
plt.plot(p_detect)

sim_df = pd.DataFrame()
sim_df["Reads_Sequenced"] = range(sim_data6_3.Reads_Mapped.max())
sim_df["p_detect"] = sim_p_detect

sim_df.to_csv("Sim6_Hypergeometric Comparison/Simulation 1000x Sensitivity by Reads Mapped.csv")


