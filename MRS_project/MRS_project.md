Alignments of assemblies. Assemblies done by Jordan R. Price in NIAB.

1. Mummer alignment of MRS strains versus reference assembly A22. 

```bash
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
for Query in $(ls niab_assemblies/*/ab1*/*.fasta); do
Strain=$(echo $Query | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_A22
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```
(For generation of circos plots check specific folder).

2. Minimap2 alignment for read depth plots.

Matt ran minimap2 himself and gave me the alignments generated to process them for circos. 

For generation of the .tsv files required for circos and IGV for ab1095, ab1907 and ab1098 versus reference A22.

```bash
for Sam in $(ls analysis/genome_alignment/minimap/*/vs_ab109*/*_aligned_sorted.bam); do
  Target=$(echo $Sam | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
  echo "$Strain-$Target"
  OutDir=$(dirname $Sam)
  samtools depth -aa $Sam > $OutDir/${Strain}_${Target}_depth.tsv
done

for Strain in ab1095, ab1097, ab1098; do
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_ab109*/*_depth.tsv); do
    echo ${Cov} | cut -f4,5,6 -d '/' --output-delimiter " - "
    cat $Cov | cut -f3 | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
  done
done > analysis/genome_alignment/minimap/read_coverage.txt
```
```bash
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_ab109*/*_depth.tsv); do
    Strain=$(echo $Cov | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Cov | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Cov)
    ProgDir=git_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_depth.tsv > $OutDir/${Organism}_${Strain}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_depth_10kb.tsv
  done
```
(For circos plots generation check specific folder)


