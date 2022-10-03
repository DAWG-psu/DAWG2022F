# DAWG workshop - Oct 7th
## Processing sequence data into abundance matrix with taxonomy identification

## Background

### From sample to sequence to data

#### Type of sequencing ####
![image](https://user-images.githubusercontent.com/77017866/193500525-9bb0e0c1-0473-4c4b-8544-ee5e54297c1a.png )

**Amplicon sequencing**
Amplicon sequencing of marker-genes (e.g. 16S, 18S, ITS) involves using specific primers that target a specific gene or gene fragment. It is used for taxonomic classification and community structure changes in response to a treatment, environmental gradients and time. Also, this technology can focus on microbial DNAs extracted from sample (avoid host DNA contamination)

**Metagenomics**
Shotgun metagenomic sequencing aims to amplify all the accessible DNA of a mixed community. It uses random primers and therefore suffers much less from pcr bias. Metagenomics provides a window into the taxonomy and functional potential of a sample. However, it takes more data to complete, and if your sample contains complicated host DNA, one would end up discarding too much data


#### Sequencing platform ####
![image](https://user-images.githubusercontent.com/77017866/193501230-0921699c-b36a-4fea-88b5-10d316db19a8.png)

#### Sampling strategy and collection metadata ####

1. What is your comparison (continuous value? groups? treatment? spatial and/or temporal variations?)
2. How many samples will be collected?
3. Including controls (i.e., DNA extraction control, sampling control, library control, etc..)

![image](https://user-images.githubusercontent.com/77017866/193501804-eccb70c6-cbab-4b80-b009-3e274f10b336.png)

##### General pipeline of microbiome study

1. Sampling
2. DNA extraction/Library preparation
3. Sequencing
4. Data generation
5. Downstream analysis

#### More information about Microbiome center KickStart Workshop at https://github.com/Penn-State-Microbiome-Center/KickStart-Workshop-2022

## Background

## DADA2 pipeline
#### What is DADA2?
