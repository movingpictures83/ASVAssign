# ASVAssign
# Language: R
# Input: TXT (keyword, value pairs)
# Output: prefix
# Tested with: PluMA 1.1, R 4.0.0
# dada2_1.18.0

PluMA plugin to take a collection of reads and map them to Amplicon
Sequence Variants (ASVs) (Callahan et al, 2016).  The reads are assumed
to have been processed using the DADA2 R package (Callahan et al, 2017)
 
The plugin accepts as input a TXT file with (keyword, value) pairs.  Three keywords are accepted:

nochim: Sequences with chimeras removed
DB: Reference database
specDB: Supplemental reference database for species.

The output prefix is used to produce three output files:

<prefix>.fa: Sequence corresponding to each ASV
<prefix>.counts.tsv: Amount of each ASV in each sample
<prefix>.taxonomy.tsv: Taxonomic classification of each ASV
