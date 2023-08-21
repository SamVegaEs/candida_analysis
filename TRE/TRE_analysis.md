#Analysis of assemblies re-arragements by TRE. 

Promer alignment of Assemblies. Assemblies constructed by Jordan at NIAB.

1. Alignment of our assemblies against the reference A22 C. albicans

```bash
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
for Query in $(ls niab_assemblies/c.albicans/UoK_July2023/assemblies/*/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f5 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_A22
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```

2. Coverage of Nanopore reads over Assembly. For the alignment of ONT reads versus Nanopore assembly use the program minimap:


```bash
for Assembly in $(ls niab_assemblies/c.albicans/UoK_July2023/assemblies/*/*.fasta); do
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
Strain=$(echo $Assembly | rev | cut -f2 -d '/'| rev)
Organism=$(echo $Reference | rev | cut -f2 -d '/' | rev)
Reads=$(ls niab_assemblies/c.albicans/UoK_July2023/porechop/$Strain/*.fastq.gz)
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/minimap/$Organism/vs_${Strain}
mkdir $OutDir
ProgDir=git_repos/tools/seq_tools/genome_alignment
sbatch $ProgDir/minimap/slurm_minimap2_2.sh $Reference $Reads $OutDir
done
```


