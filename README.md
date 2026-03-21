# Taxonomic classification of shotgun metagnomics data from gut microbiomes of vegans and omnivores

# 1 | Introduction

# 2 | Methods

## 2.1 | Data Description
```
prefetch SRR8146970 SRR8146975 SRR8146971 SRR8146989 SRR8146974 SRR8146986
fasterq-dump SRR*
```

## 2.2 | Quality Control
```
fastqc -t 8 genomes/*.fastq -o fastqc_results
multiqc fastqc_results/ -o multiqc_results
```

## 2.3 | Taxonomic Classfication
```
kraken-biom bracken_results/*.report -o table.biom 
```

## 2.4 | Diversity Measures

## 2.5 | Differential Abundance

## 2.6 | Visualization

# 3 | Results

# 4 | Discussion

# References
De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., ... & Ercolini, D. (2019). Distinct genetic and functional traits of human intestinal Prevotella copri strains are associated with different habitual diets. Cell host & microbe, 25(3), 444-453.
