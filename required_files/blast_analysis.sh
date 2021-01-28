#!/bin/bash

# Makes BLAST database from all C. albcians protein sequences provided in FASTA format in a single file
makeblastdb -in C_albicans_A22_proteins.fa -out C_albicans_db -dbtype prot

# BLAST search with the C. albicans database and all C. auris proteins provided in FASTA format in a single file
blastp -db C_albicans_db -query C_auris_proteins.fasta -out C_auris_blast_results.txt



