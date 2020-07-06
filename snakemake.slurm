#!/bin/bash -l
#SBATCH -J snakemake
#SBATCH -N 1 #number of nodes to be used
#SBATCH --ntasks-per-node=20
#SBATCH -p ulimited2 #name of the job queue

module load shared
module load shared bwa/0.7.17-r1194-dirty
module load shared samtools/1.9.0
module load shared gatk/4.0.12.0

snakemake -j 20 --cluster-config cluster_slurm.yaml --cluster "sbatch -A {cluster.account} -p {cluster.partition} --output {cluster.output} --error {cluster.error} --nodes {1} --ntasks {20}"
 

exit
