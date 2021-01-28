# RNAseq_analysis_Cauris

RNAseq analysis pipeline for analysis of Illumina short-read sequencing data from Candida auris.

C. auris genome sequence and annotation file (obtained from NCBI) are provided. For current version see NCBI.

Clone the repository and add the read files in the base directory. To start the analysis use: `bash analysis_script.sh`

The `blast_analysis.sh` script creates a BLAST database from C. albicans protein fasta sequences and performs BLASTp with all C. auris proteins supplied in fasta format on this database to find homologues.

