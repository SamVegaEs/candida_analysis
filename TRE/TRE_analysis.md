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

3. SNP Calling of strains.

3.1 First create  a custom SnpEff genome database as SC5314 has still not been created as a reference. 

```bash
SnpEff=/home/sv264/local/bin/snpEff
nano $SnpEff/snpEff.config
```
```
Add the following lines to the section with databases:
#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
# Bc16 genome
Bc16v1.0.genome: BC-16
# P414 genome
P414v1.0.genome: 414
# 62471 genome
62471v1.0.genome: 62471
# R36_14 genome
R36_14v1.0.genome: R36_14
# SCRP371 genome
SCRP371v1.0.genome: SCRP371
# P. stipis
Ps589v1.0.genome: 589
PsY-11545v1.0.genome: Y-11545
PsCBS6054v1.0.genome: CBS6054
PsAB108.genome: ab1082
PsAB918.genome: ab918
# C. albicans 
Calbicans.genome: SC5314
```

Collect input files

```bash
Organism="C.albicans"
Strain="SC5314"
DbName="Calbicans"
ProjDir=$PWD
Reference=$(ls $ProjDir/assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
Gff=$(ls $ProjDir/assembly/misc_publications/c.albicans/*.gff)
SnpEff=/home/sv264/local/bin/snpEff
mkdir $SnpEff/data/${DbName}
cp $Reference $SnpEff/data/${DbName}/sequences.fa
cp $Gff $SnpEff/data/${DbName}/genes.gff

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v ${DbName}
```



