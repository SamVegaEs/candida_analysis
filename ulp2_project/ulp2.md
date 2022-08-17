Alignment of Assemblies (sofmasked) for ULP2 project.

Assemblies done by Jordan R. Price from NIAB.

1. Promer alignment of Assemblies against the reference A22 C. albicans

For ab791 and ab758

```bash
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
for Query in $(ls niab_assemblies/*/*/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_A22
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
For ab863

```bash
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
for Query in $(ls niab_assemblies/*/ab863/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_A22
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
2. BWA alignment of Illumina raw reads versus A22 for read coverage plots.  


Reads for A22 Illumina coverage were obtained from: https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&page=4&page_size=10&acc=SRR850113&display=download
https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r97#Sec10

```bash
for Reference in $(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta); do
for StrainPath in $(ls -d raw_dna/paired/c.albicans/SC5314); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
Ref=$(echo $Reference | rev | cut -f1 -d '/' | rev)
F_Read=$(ls $StrainPath/*.fastq.gz)
R_Read=$(ls $StrainPath/*.fastq.gz)
echo $F_Read
echo $R_Read
echo $Ref
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Ref}
ProgDir=~/git_repos/tools/seq_tools/genome_alignment/bwa
sbatch $ProgDir/sub_bwa_slurm.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

Generation of the .tsv files required for circos and for IGV.

```bash
for Sam in $(ls analysis/genome_alignment/minimap/*/SC5314/*_aligned_sorted.bam); do
  Target=$(echo $Sam | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
  echo "$Strain-$Target"
  OutDir=$(dirname $Sam)
  samtools depth -aa $Sam > $OutDir/${Strain}_${Target}_depth.tsv
done

for Strain in ab55 ab791; do
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_*/*_depth.tsv); do
    echo ${Cov} | cut -f4,5,6 -d '/' --output-delimiter " - "
    cat $Cov | cut -f3 | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
  done
done > analysis/genome_alignment/minimap/read_coverage.txt
```
```bash
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_ab*/*_depth.tsv); do
    Strain=$(echo $Cov | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Cov | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Cov)
    ProgDir=git_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_depth.tsv > $OutDir/${Organism}_${Strain}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_depth_10kb.tsv
  done
```
(For circos scripts check specific folder).

3. Minimap2 alignment for read coverage plot agains reference A22, using ONT data (Unmasked assembly).

For ab791 and ab758

```bash
for Assembly in $(ls niab_assemblies/*/*/*.fasta); do
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
Strain=$(echo $Assembly | rev | cut -f2 -d '/'| rev)
Organism=$(echo $Reference | rev | cut -f3 -d '/' | rev)
Reads=$(ls projects/niab_assemblies/*/$Strain/*.fastq.gz)
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/minimap/$Organism/vs_${Strain}
ProgDir=git_repos/tools/seq_tools/genome_alignment
sbatch $ProgDir/minimap/slurm_minimap2.sh $Reference $Reads $OutDir
done
```
For ab863

```bash
for Assembly in $(ls niab_assemblies/*/ab863/*.fasta); do
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
Strain=$(echo $Assembly | rev | cut -f2 -d '/'| rev)
Organism=$(echo $Reference | rev | cut -f3 -d '/' | rev)
Reads=$(ls projects/niab_assemblies/*/$Strain/*.fastq.gz)
Prefix="${Organism}_${Strain}"
OutDir=~/analysis/genome_alignment/minimap/$Organism/vs_${Strain}
ProgDir=git_repos/tools/seq_tools/genome_alignment
sbatch $ProgDir/minimap/slurm_minimap2.sh $Reference $Reads $OutDir
done
```

Generation of the .tsv files required for circos and for IGV.

```bash
for Sam in $(ls analysis/genome_alignment/minimap/*/vs_ab*/*_aligned_sorted.bam); do
  Target=$(echo $Sam | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
  echo "$Strain-$Target"
  OutDir=$(dirname $Sam)
  samtools depth -aa $Sam > $OutDir/${Strain}_${Target}_depth.tsv
done

for Strain in ab55 ab791; do
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_*/*_depth.tsv); do
    echo ${Cov} | cut -f4,5,6 -d '/' --output-delimiter " - "
    cat $Cov | cut -f3 | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
  done
done > analysis/genome_alignment/minimap/read_coverage.txt
```
```bash
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_ab*/*_depth.tsv); do
    Strain=$(echo $Cov | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Cov | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Cov)
    ProgDir=git_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_depth.tsv > $OutDir/${Organism}_${Strain}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_depth_10kb.tsv
  done
```
(For circos scripts refer to specific folder)

To visualise in IGV the new bam files we need the index of that bam file. Minimap did not generate the index files of the alignments, so I will try to index them using samtools.


```bash
for File in $(ls analysis/genome_alignment/minimap/*/vs_ab791/*.fasta_aligned_sorted.bam); do
Strain=$(echo $File | rev | cut -f2 -d '/'| rev)
File=analysis/genome_alignment/minimap/*/$Strain/*.fasta_aligned_sorted.bam
samtools index <$File> -o $File
done
```

In order to extract the reads from the bam files we will use samtools. Example:

```bash
#Run at cd /analysis/genome_alignment/minimap/c.albicans/vs_XX. The -h option was added, because when sorting a message error appeared saying that the header was missing.

samtools view A22_Chromosomes_29.fasta_aligned_sorted.bam "Ca22chr1A_C_albicans_SC5314:1293000-1293500" -h > ab863_C1_2.bam

#Sorting of the bam file generated.

samtools sort -n ab863_C1_2.bam > ab863_C1_2_sorted.bam

#Now that we have the reads in a bam file, we can extract the fastq reads using bedtools.

bedtools bamtofastq -i ab863_C1_2_sorted.bam -fq ab863_C1_2_sorted.bam.fq

samtools view A22_Chromosomes_29.fasta_aligned_sorted.bam "Ca22chrRA_C_albicans_SC5314" | awk '$10 ~/ACTTCTTGGTGTACGGATGTCTA/' | awk '{print $1, $10}' > telseq.fa
```

4. Plot the coverage of Nanopore reads in the minichromosome region using R. Example:

```
#Extract the reads for ab863 and ab791 respect ab55 in Chromosome R and in Chromosome 1:

cat analysis/genome_alignment/minimap/c.albicans/vs_ab55/c.albicans_vs_ab55_depth.tsv  | grep 'Ca22chrRA_C_albicans_SC5314' > tmp_CR_ab55.tsv

cat analysis/genome_alignment/minimap/c.albicans/vs_ab791/c.albicans_vs_ab791_depth.tsv | grep 'Ca22chrRA_C_albicans_SC5314' > tmp_CR_ab791.tsv

cat analysis/genome_alignment/minimap/c.albicans/vs_ab863/c.albicans_vs_ab863_depth.tsv  | grep 'Ca22chrRA_C_albicans_SC5314' > tmp_CR_ab863.tsv


/Users/sv264/OneDrive - University of Kent/Desktop/tmp_CR_ab55.ts
#Now we open the folders in R.

tmp_ab55 <- read.delim("C:/Users/sv264/OneDrive - University of Kent/Desktop/tmp_CR_ab55.tsv", header=FALSE)
tmp_ab791 <- read.delim("C:/Users/sv264/OneDrive - University of Kent/Desktop/tmp_CR_ab791.tsv", header=FALSE)
tmp_ab863 <- read.delim("C:/Users/sv264/OneDrive - University of Kent/Desktop/tmp_CR_ab863.tsv", header=FALSE)

library(ggplot2)

p <- ggplot(data=tmp_ab55, aes(x=tmp_ab55$V2, y=tmp_ab55$V3)) + geom_line() + labs(x = "Position (bp)", y = "Coverage")+ geom_vline(xintercept = 1300000, colour = 'red', linetype = "dashed") + geom_vline(xintercept = 2550000, colour = 'red', linetype = "dashed") 
outfile= paste("ab55", "vs_ab55", "chromosome1", "minion.jpg", sep = "_")
ggsave(outfile , plot = p, width = 20, height = 5, units = 'in', limitsize = TRUE)

p2 <- ggplot(data=tmp_ab791, aes(x=tmp_ab791$V2, y=tmp_ab791$V3)) + geom_line() + labs(x = "Position (bp)", y = "Coverage") 
outfile= paste("ab791", "vs_ab55", "chromosomeR", "minion.jpg", sep = "_")
ggsave(outfile , plot = p2, width = 20, height = 5, units = 'in', limitsize = TRUE)

Change the name of the columns before merging data

names(tmp_ab791)[names(tmp_ab791) == "V1"] <- "V4"
names(tmp_ab791)[names(tmp_ab791) == "V2"] <- "V5"
names(tmp_ab791)[names(tmp_ab791) == "V3"] <- "V6"
names(tmp_ab863)[names(tmp_ab863) == "V1"] <- "V7"
names(tmp_ab863)[names(tmp_ab863) == "V2"] <- "V8"
names(tmp_ab863)[names(tmp_ab863) == "V3"] <- "V9"

Merge all the data

tmp_all <- cbind(tmp_ab55, tmp_ab791, tmp_ab863)

All data in the same plot:

plot <- plot (tmp_ab55$V2, tmp_ab55$V3, ylim=c(0,500), xlim=c(1870000,2200000), type="l", xlab="Position", ylab = "Coverage")+lines(tmp_ab55$V2, tmp_ab791$V3, xlim=c(1870000,2200000), ylim=c(0,500), type="l", col=2)+lines(tmp_ab55$V2, tmp_ab863$V3, type="l", xlim=c(1870000,2200000), ylim=c(0,500), col=3)
```


