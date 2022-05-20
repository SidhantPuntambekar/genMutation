#This script will annotate the synonymous mutations data from the curated data file fake_transcript_variants_sorted_v3_ad_syn_processed.gz
#Sidhant Puntambekar
#Currently takes 31 minutes to run with 80000 sites in synon and all site ranges from genome in a bottle 

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import scipy.stats as stats
#%matplotlib inline

#Load and check synonymous mutations
df = pd.read_csv('../data/fake_transcript_variants_sorted_v3_ad_syn_processed', delimiter="\t")
#print("Initial Load: ")
#print(df.head())
#print(len(df))
#print(df.tail())

#Load and check genome in a bottle file
dfGRCh = pd.read_csv('../data/GRCh38/GRCh38_alldifficultregions.bed', delimiter="\t")
#print("Initial Load: ")
#print(dfGRCh.head())
#print(len(dfGRCh))
#print(dfGRCh.tail())

#dfGRCh['all_sites'] = dfGRCh['left_sites'].astype(str) + dfGRCh['right_sites'].astype(str)
#print("New Col: ")
#print(dfGRCh.head())
#print(len(dfGRCh))
#print(dfGRCh.tail())

#Filter to just chromosome 1
dfGRChChrom1 = dfGRCh[dfGRCh['chrom'] == 'chr1']
#print("New Col: ")
#print(dfGRChChrom1.head())
#print(len(dfGRChChrom1))
#print(dfGRChChrom1.tail())

#Filter based on curated chromosome 1 sites
dfGRChChrom1FilteredLeftSeries = dfGRChChrom1[dfGRChChrom1['left_sites'] >= 43947376] #43947606 lower bound of synon mutation chrom1 sites, 68232294 upper bound of synon mutation chrom1 sites
dfGRChChrom1FilteredLeftSeries = dfGRChChrom1FilteredLeftSeries[dfGRChChrom1FilteredLeftSeries['left_sites'] <= 68232294 + 5000]
print(dfGRChChrom1FilteredLeftSeries.head())
print(len(dfGRChChrom1FilteredLeftSeries))
print(dfGRChChrom1FilteredLeftSeries.tail())

#Insert column of trues for col representing if site is in difficult sequencing sites in GRCh38 reference genome
dfLength = len(df)
df.insert(7, 'site_in_genome_bottle', [False for x in range(dfLength)])
print(df.head())
print(df.tail())

#Compare sites in df and filtered dfGRChChrom1FilteredLeftSeries
siteCol = df.loc[:,'pos']
siteArray = siteCol.values

leftSites = dfGRChChrom1FilteredLeftSeries.loc[:,'left_sites']
leftSitesArray = leftSites.values

rightSites = dfGRChChrom1FilteredLeftSeries.loc[:,'right_sites']
rightSitesArray = rightSites.values

inSiteBool = np.array([False]*202176)

print("Synonymous Mutation Sites: ", siteArray)
print("Difficult Sequence Left Bounds: ", leftSitesArray)
print("Difficult Sequence Right Bounds: ", rightSitesArray)

#df = df.reset_index(drop=True)
#dfGRChChrom1FilteredLeftSeries = dfGRChChrom1FilteredLeftSeries.reset_index(drop=True)

#comparison_col = np.where(df['pos'].values >= dfGRChChrom1FilteredLeftSeries['left_sites'].values, True, False)

#df["site_is_in_genome_bottle"] = comparison_col
#print(df.head())

indexInSite = np.array([0]*202176)
j = 0

#print(len(leftSitesArray))
for i in range(len(df)):
	#for j in range(len(dfGRChChrom1FilteredLeftSeries)):
		#print(siteArray[i])
	while (leftSitesArray[j] <= siteArray[i] and j < len(leftSitesArray)):
		j = j + 1
		#print(j)	
	print("LeftSite", leftSitesArray[j-1])
	print("Synon Site", siteArray[i])
	print(j-1)	
	if (siteArray[i] >= leftSitesArray[j-1] and siteArray[i] <= rightSitesArray[j-1]):
		inSiteBool[i] = True
		#print(i) #Print index of position in synon mutation file
		indexInSite[i] = i
	else:
		inSiteBool[i] = False

for i in range(20001):
	print(indexInSite[i])
print("InSiteBool Array length: ", len(inSiteBool))
#print(len(leftSitesArray))

df.loc[indexInSite, 'site_in_genome_bottle'] = True
print(df.head(65))
print(df.tail(65))
#print(i) #Print location in array where inSiteBool array contains a true
#df.to_csv('synonAnnotated.csv', index=False)
