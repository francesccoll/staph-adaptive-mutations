# The mutational landscape of *Staphylococcus aureus* during colonisation

This GitHub project contains the data and code necessary to reproduce the findings of the study 'The mutational landscape of *Staphylococcus aureus* during colonisation', and includes the following directories:
* adaptive_mutations_analysis: contains multiple directories with scripts executing genomic analyses:
  * applied to all isolates:
    - amrfinder: Files and scripts used to run AMRFinderPlus from *S. aureus* assemblies.
    - assembly: Files and script to run *de novo* assembly from Illumina short reads.
    - mapping_ref: scripts used to map Illumina short reads to *S. aureus* reference genomes.
    - mlst: scripts to run *in silico* MLST typing.
  * applied to multiple isolates per CC
    - pairwise_snp_distances: scripts to create pairwise SNP matrix for all isolates (across CCs)
    - phylogeny_per_cc: scripts to generate core-genome phylogenetic per clonal complex.
  * applied to isolates of the same host
    - identify_ancestral: pipeline to re-construct the MRCA sequence of each host clonal strain.
    - mapping_mrca: pipeline used to identify *de novo* mutations arising in each host clonal strain.
  * applied to all hosts
    - association_mrca: R scripts used to filter *de novo* mutations per host clonal strain and to identify functional units with an excess of protein-altering mutations.
* agr_mutants: R script is used to calculate agr-mutant frequencies in different subsets of patients and to test the effect of several variables on the emergence of agr mutants.
* data: directory with input data files used by scripts in this GitHub project.
* genetic_diversity: scripts to summarise within-host genetic diversity.
* growth_curves: growth curves data of knockout strains and natural mutants grown under sub-inhibitory concentrations of daptomycin, and R script used to analyse this data.
* growth_curves_nitrogen: growth curves data of transposon knockout strains and natural mutants grown under multiple nitrogen sources in Biolog PM3 plates, and R script used to analyse this data.

# Citation
Coll F, Blane B, Bellis K, *et al.* The mutational landscape of *Staphylococcus aureus* during colonisation. *Nature Communications* (accepted). 2024. 
Previous BioRxiv version available here: https://www.biorxiv.org/content/10.1101/2023.12.08.570284v1
