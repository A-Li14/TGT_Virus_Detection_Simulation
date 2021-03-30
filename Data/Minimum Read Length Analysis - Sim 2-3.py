
# files = os.listdir("..")
import pandas as pd
from matplotlib import pyplot as plt

# file_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES10%ReadErrorRate_10power2VirusHostRatio_Simulation1/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation1_RVDB_result.csv"

##Simulation 2
# file_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES10%ReadErrorRate_10power2VirusHostRatio_Simulation2/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation2_RVDB_result.csv"
# output_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/Analyzing_Results/Sim2_Accuracy_By_Length/Sim2_10error_Alignment_Rate.csv"

##Simulation 2.1
# file_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES15%ReadErrorRate_10power2VirusHostRatio_Simulation2.1/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation2.1_RVDB_result.csv"
# output_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/Analyzing_Results/Sim2_Accuracy_By_Length/Sim2.1_15error_Alignment_Rate.csv"

##Simulation 2.2
# file_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES20%ReadErrorRate_10power2VirusHostRatio_Simulation2.2/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation2.2_RVDB_result.csv"
# output_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/Analyzing_Results/Sim2_Accuracy_By_Length/Sim2.2_20error_Alignment_Rate.csv"

##Simulation 3
file_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES10%ReadErrorRate_10power2VirusHostRatio_Simulation3/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation3_RVDB_result.csv"
output_path = "/home/gridsan/alexsli/Alexander_Notebook/Simulation2/Analyzing_Results/Sim3_Accuracy_By_Length/Sim3_10error_Alignment_Rate.csv"


df = pd.read_csv(file_path)    
df.head()
df.shape

###two criteria for finding MVM alignments: mvm or minute in the hit definition
df['sbjct_title'].str.lower().str.find("mvm")
sum(df['sbjct_title'].str.lower().str.find("mvm")!=-1) ###Count alignments to MVM variants
max(df['query_id'].str.strip("read").astype("int"))  ###Count the number of reads aligned

df['sbjct_title'].str.lower().str.find("minute")
sum(df['sbjct_title'].str.lower().str.find("minute")!=-1) ###Count alignments to MVM variants
max(df['query_id'].str.strip("read").astype("int"))  ###Count the number of reads aligned

###histogram of the alignment lengths
df.columns
min(df.align_length)
plt.hist(df.align_length,bins=30)



#####using read length data

###read the csv file containing read lengths; compare to the query lengths for alignments
df_lengths = pd.read_csv(file_path.strip("_RVDB_result.csv")+".csv",header=None)
df_lengths.head()

###subset df_lengths to the reads that aligned against mvm or minute
matched = (df['sbjct_title'].str.lower().str.find("minute")!=-1) | (df['sbjct_title'].str.lower().str.find("mvm")!=-1)
df[matched].shape
df[matched]

read_lengths = df_lengths[1].value_counts()
df[df["query_id"].str.contains("read30933")]


# max(read_lengths)
matched_lengths = df[matched].drop_duplicates(subset=["query_id"])['query_length'].value_counts()

# length_bins = range(0,4000,100) ###Specify how to bin read lengths for measuring alignment accuracy
# length_bins = [0,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1250,1500,2000,4000]
length_bins = list(range(0,4000,50))
# temp = read_lengths.groupby(pd.cut(read_lengths.index, length_bins)).sum()
match_df = pd.DataFrame(data = {"Total_Reads":read_lengths.groupby(pd.cut(read_lengths.index, length_bins)).sum(),
    "Aligned_Reads":matched_lengths.groupby(pd.cut(matched_lengths.index, length_bins)).sum()})

match_df["Sensitivity"] = match_df["Aligned_Reads"]/match_df["Total_Reads"]
match_df


plot = match_df.plot(y="Sensitivity",xlabel="Read Length",ylabel="Proportion",title="Sensitivity by Read Length")
plot.get_figure().savefig("Sim3_Accuracy_By_Length/Alignment Rate by Read Length.png")
# plt.plot([i for i in match_df.index.categories.to_series()],match_df["Alignment_Rate"])
# plt.xlabel("Read Length")
# plt.xticks(match_df.index)
# plt.show()
match_df.to_csv(output_path)


##########----------------------------------------------------------------

###Plot classification accuracy for different error rates

ds15 = pd.read_csv("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES15%ReadErrorRate_10power2VirusHostRatio_Simulation3.1/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation3.1_RVDB_result.csv")
ds15_lengths = pd.read_csv("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES15%ReadErrorRate_10power2VirusHostRatio_Simulation3.1/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation3.1.csv",header=None)

ds20 = pd.read_csv("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES20%ReadErrorRate_10power2VirusHostRatio_Simulation3.2/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation3.2_RVDB_result.csv")
ds20_lengths = pd.read_csv("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/MVM_100Threads_YES20%ReadErrorRate_10power2VirusHostRatio_Simulation3.2/MVM_FragmentsCountLength_100Threads_10power2VirusHostRatio_Simulation3.2.csv",header=None)


###subset ds15_lengths to the reads that aligned against mvm or minute
matched = (ds15['sbjct_title'].str.lower().str.find("minute")!=-1) | (ds15['sbjct_title'].str.lower().str.find("mvm")!=-1)
ds15[matched].shape
ds15[matched]

read_lengths = ds15_lengths[1].value_counts()

# max(read_lengths)
matched_lengths = ds15[matched].drop_duplicates(subset=["query_id"])['query_length'].value_counts()

# length_bins = range(0,4000,100) ###Specify how to bin read lengths for measuring alignment accuracy
# length_bins = [0,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1250,1500,2000,4000]
length_bins = list(range(0,4000,50))
# temp = read_lengths.groupby(pd.cut(read_lengths.index, length_bins)).sum()
match_ds15 = pd.DataFrame(data = {"Total_Reads":read_lengths.groupby(pd.cut(read_lengths.index, length_bins)).sum(),
    "Aligned_Reads":matched_lengths.groupby(pd.cut(matched_lengths.index, length_bins)).sum()})

match_ds15["Sensitivity"] = match_ds15["Aligned_Reads"]/match_ds15["Total_Reads"]
match_ds15

####-------------------------
###subset ds20_lengths to the reads that aligned against mvm or minute
matched = (ds20['sbjct_title'].str.lower().str.find("minute")!=-1) | (ds20['sbjct_title'].str.lower().str.find("mvm")!=-1)
ds20[matched].shape
ds20[matched]

read_lengths = ds20_lengths[1].value_counts()

# max(read_lengths)
matched_lengths = ds20[matched].drop_duplicates(subset=["query_id"])['query_length'].value_counts()

# length_bins = range(0,4000,100) ###Specify how to bin read lengths for measuring alignment accuracy
# length_bins = [0,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1250,1500,2000,4000]
length_bins = list(range(0,4000,50))
# temp = read_lengths.groupby(pd.cut(read_lengths.index, length_bins)).sum()
match_ds20 = pd.DataFrame(data = {"Total_Reads":read_lengths.groupby(pd.cut(read_lengths.index, length_bins)).sum(),
    "Aligned_Reads":matched_lengths.groupby(pd.cut(matched_lengths.index, length_bins)).sum()})

match_ds20["Sensitivity"] = match_ds20["Aligned_Reads"]/match_ds20["Total_Reads"]
match_ds20



match_df.plot(y="Sensitivity",xlabel="Read Length",ylabel="Proportion",title="Sensitivity by Read Length")
match_ds15.plot(y="Sensitivity",xlabel="Read Length",ylabel="Proportion",title="Sensitivity by Read Length")
match_ds20.plot(y="Sensitivity",xlabel="Read Length",ylabel="Proportion",title="Sensitivity by Read Length")


match_df["Error_Rate"] = "10%"
match_ds15["Error_Rate"] = "15%"
match_ds20["Error_Rate"] = "20%"


sensitivity_df = pd.concat([match_df,match_ds15,match_ds20])
sensitivity_df.plot(y="Sensitivity",xlabel="Read Length",ylabel="Proportion",title="Sensitivity by Read Length")

sensitivity_df.to_csv("/home/gridsan/alexsli/Alexander_Notebook/Simulation2/Analyzing_Results/Sim3_Accuracy_By_Length/Classification Accuracy by Read Length.csv")


# fig, ax = plt.subplots(figsize=(3,1))
# # use unstack()
# sensitivity_df.groupby(['Error_Rate'])['Sensitivity'].plot(ax=ax)

# # fig, axes = plt.subplots(ncols=2)
# sensitivity_df[sensitivity_df.Error_Rate=="10%"].Sensitivity.plot(ax=axes[0])
# sensitivity_df[sensitivity_df.Error_Rate=="10%"].Sensitivity.plot(ax=axes[0])
# sensitivity_df[sensitivity_df.Error_Rate=="10%"].Sensitivity.plot(ax=axes[0])





# reads = set(df['query_id'].str.strip("read").astype("int"))
# check = [i in reads for i in range(max(reads))]
# missing = [i for i in range(max(reads)) if check[i]==False]

# missing
# 1 in reads
# sum(df['query_id'].str.find("read50")!=-1)
# df['query_id'=="read50"]


# def extract_species(blast_df,start_ind=15) :
#     ind1 = blast_df['sbjct_title'].str.find("|",start_ind)
#     output = []
#     for i in range(len(df)) :
#         ind2 = blast_df['sbjct_title'].str.find("|",ind1[i])
#         output.append(blast_df['sbjct_title'])
#     ind2 = blast_df['sbjct_title'].str.find("|",max(ind1))




# df['sbjct_title'].str.find("|",15)
