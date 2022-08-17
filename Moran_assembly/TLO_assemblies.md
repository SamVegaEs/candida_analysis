Project started for Minion assembly for Candida tlo project (Moran collaboration). 


#Identify sequencing coverage

Allow perl to read the content of the programs: sudo chmod a=rx /git_repos/tools/seq_tools/dna_qc

For Minion data:


```bash
for RawData in $(ls projects/candida_tlo/raw_data_minion/*/*/cc*.gz); do
echo $RawData;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
GenomeSz=14
OutDir=$(dirname $RawData)
mkdir -p $OutDir
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done


  for StrainDir in $(ls -d projects/candida_tlo/raw_data_minion/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```bash
for RawData in $(ls projects/candida_tlo/raw_data_minion/*/*/BC*.gz); do
echo $RawData;
ProgDir=~/git_repos/tools/seq_tools/dna_qc;
GenomeSz=14
OutDir=$(dirname $RawData)
mkdir -p $OutDir
sbatch $ProgDir/sub_count_nuc_slurm_2.sh $GenomeSz $RawData $OutDir
done


  for StrainDir in $(ls -d projects/candida_tlo/raw_data_minion/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```




```
Coverage.

cc03    126.04
cc10    144
cc12    180.43
cc16    183.98

BC04    15.36
BC05    14.94
BC06    13.97
```

#Read correction using Canu

Run independently worked, but not all together with the same code. 

```bash
for TrimReads in $(ls projects/candida_tlo/raw_data_minion/*/*/cc*q.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.8/$Organism/"$Strain"
ProgDir=~/git_repos/tools/seq_tools/assemblers/canu
sbatch $ProgDir/sub_canu_correction_slurm.sh $TrimReads 14m $Strain $OutDir
done
```

```bash
for TrimReads in $(ls projects/candida_tlo/raw_data_minion/*/*/BC*.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.8/$Organism/"$Strain"
ProgDir=~/git_repos/tools/seq_tools/assemblers/canu
sbatch $ProgDir/sub_canu_correction_slurm_2.sh $TrimReads 14m $Strain $OutDir
done
```

#Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.8/*/*/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=~/git_repos/tools/seq_tools/assemblers/SMARTdenovo
sbatch $ProgDir/sub_SMARTdenovo_slurm.sh $CorrectedReads $Prefix $OutDir
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=~/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  sbatch $ProgDir/sub_quast_slurm.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=~/git_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d ~/git_repos/tools/gene_prediction/busco/saccharomycetes_odb10)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
sbatch $ProgDir/sub_busco3_slurm.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```bash	
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo.dmo.lay.txt      1168     11      569     400     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo.dmo.lay.txt      1141     21      594     402     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo.dmo.lay.txt      1198     34      546     393     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo.dmo.lay.txt      1167     15      579     391     2137
```

#Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls projects/candida_tlo/raw_data_minion/*/$Strain/*q.gz)
# ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon2_$Iterations"
ProgDir=~/git_repos/tools/seq_tools/assemblers/racon
sbatch $ProgDir/sub_racon_slurm_2.sh $Assembly $ReadsFq $Iterations $OutDir
done
``` 

I skipped the remove contaminants bits, at it seems specific for the requirements NCBI has when the genomes are submitted.
How neccessary is quast, if you use BUSCO? 

```bash
ProgDir=~/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*16/racon2_10/*.fasta | grep 'round_10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
sbatch $ProgDir/sub_quast_slurm.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/cc16/racon2_10/*_round_9.fasta); do
echo $Assembly
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=~/git_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d ~/git_repos/tools/gene_prediction/busco/saccharomycetes_odb10)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
sbatch $ProgDir/sub_busco3_slurm.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/cc16/*/assembly/*/short_summary*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
cc03

short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_1.txt         1466    13      466     205     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_2.txt         1491    16      429     217     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_3.txt         1489    11      439     209     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_4.txt         1503    12      428     206     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_5.txt         1451    8       463     223     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_6.txt         1483    9       446     208     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_7.txt         1451    10      468     218     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_8.txt         1460    16      470     207     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_9.txt         1489    11      427     221     2137
short_summary.specific.saccharomycetes_odb10.cc03_smartdenovo_racon_round_10.txt        1474    12      431     232     2137

cc10

short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_1.txt         1438    54      503     196     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_2.txt         1466    39      474     197     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_3.txt         1450    51      491     196     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_4.txt         1463    46      475     199     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_5.txt         1493    45      457     187     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_6.txt         1486    45      459     192     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_7.txt         1475    44      448     214     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_8.txt         1476    43      471     190     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_9.txt         1462    40      486     189     2137
short_summary.specific.saccharomycetes_odb10.cc10_smartdenovo_racon_round_10.txt        1490    48      459     188     2137

cc12

short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_1.txt         1448     57      477     212     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_2.txt         1507     61      451     179     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_3.txt         1491     53      455     191     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_4.txt         1488     62      457     192     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_5.txt         1486     51      467     184     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_6.txt         1464     58      484     189     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_7.txt         1477     52      451     209     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_8.txt         1468     55      463     206     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_9.txt         1510     47      429     198     2137
short_summary.specific.saccharomycetes_odb10.cc12_smartdenovo_racon_round_10.txt        1510     52      434     193     2137

cc16

short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_1.txt         1467     14      462     208     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_2.txt         1484     21      452     201     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_3.txt         1470     28      461     206     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_4.txt         1477     16      459     201     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_5.txt         1478     19      452     207     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_6.txt         1489     21      456     192     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_7.txt         1477     23      447     213     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_8.txt         1467     26      463     207     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_9.txt         1491     20      448     198     2137
short_summary.specific.saccharomycetes_odb10.cc16_smartdenovo_racon_round_10.txt        1479     21      453     205     2137
```

codes to rewrite:

#Remove mitochondrial DNA

Create the database for Candida mitochondrial DNA.

Splitting sequences by long repeats of ambiguous base N

```bash
 cat calb_mtDNA.fa | perl -p -e 's/N\n/N/' | perl -p -e 's/^N+//;s/N+$//;s/N{200,}/\n>split\n/' >calb_mtDNA_split.fa
```

```bash
./bwa64 index -p calb_mtDNA -a is calb_mtDNA_split.fa >bwa.log 2>&1 &
```

Using an exclusion database with deconseq:

```bash
 for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/*.fasta | grep 'round_10'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "calb_mtDNA"; do
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      sbatch $ProgDir/sub_deconseq_no_retain_slurm.sh $Assembly $Exclude_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
for Exclude_db in "calb_mtDNA"; do
echo $Exclude_db
for File in $(ls assembly/*/*/*/*/log.txt | grep "$Exclude_db"); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Bad\n";
done
done
```

```bash
calb_mtDNA
cc03    12      0
cc10    16      1
cc12    15      1
cc16    17      3
```
#Repeat Masking.

```bash
 # for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta); do
  for Assembly in $(ls assembly/SMARTdenovo/*/*/deconseq_calb_mtDNA/contigs_min_500bp_filtered_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
    ProgDir=git_repos/tools/seq_tools/repeat_masking
    sbatch $ProgDir/rep_modeling.sh $Assembly $OutDir
  done
```

sbatch $ProgDir/transposonPSI.sh $Assembly $OutDir

The TransposonPSI masked bases is normally used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files. I did not use it for candida assemblies as there was an error in the contigs naming. 

```bash
 # for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta); do
  for Assembly in $(ls assembly/SMARTdenovo/*/*/deconseq_calb_mtDNA/contigs_min_500bp_filtered_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
    ProgDir=git_repos/tools/seq_tools/repeat_masking
    sbatch $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```



#Alignment of Assemblies (hardmasked)


Promer alignment of Assemblies.

Alignment of our assemblies against the reference A22 C. albicans

```bash
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_A22
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```

Generating circos plots with Mummer alignment.

```bash
ProgDir=/~/git_repos/scripts/circos_plots/genome_alignment/C.albicans/TLO_project/vs_A22/cc03
circos -conf $ProgDir/A22_vs_cc03_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/589_vs_Y-11545_v2_circos.png
mv $OutDir/circos.svg $OutDir/589_vs_Y-11545_v2_circos.svg
ls $PWD/$OutDir/589_vs_Y-11545_v2_circos.png
```

#Alignment of Assemblies (sofmasked)


Promer alignment of Assemblies.

Alignment of our assemblies against the reference A22 C. albicans

```bash
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_A22
OutDir=analysis/genome_alignment/mummer/sofmasked/$Organism/$Strain/$Prefix
ProgDir=git_repos/tools/seq_tools/genome_alignment/MUMmer
sbatch $ProgDir/sub_nucmer_slurm.sh $Reference $Query $Prefix $OutDir
done
```



#Read coverage. Unmasked assembly.

Coverage of Nanopore reads over Assembly. For the alignment of ONT reads versus Nanopore assembly use the program minimap:


```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Reference=$(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/'| rev)
Organism=$(echo $Reference | rev | cut -f2 -d '/' | rev)
Reads=$(ls projects/candida_tlo/*/*/$Strain/*.fastq.gz)
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/minimap/$Organism/vs_${Strain}
ProgDir=git_repos/tools/seq_tools/genome_alignment
sbatch $ProgDir/minimap/slurm_minimap2.sh $Reference $Reads $OutDir
done
```

To generate the coverage plots with circos we need to generate the .tsv files

```bash
for Sam in $(ls analysis/genome_alignment/minimap/*/vs_*/*_aligned_sorted.bam); do
  Target=$(echo $Sam | rev | cut -f3 -d '/' | rev)
  Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
  echo "$Strain-$Target"
  OutDir=$(dirname $Sam)
  samtools depth -aa $Sam > $OutDir/${Strain}_${Target}_depth.tsv
done

for Strain in cc03 cc10 cc12 cc16; do
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_*/*_depth.tsv); do
    echo ${Cov} | cut -f4,5,6 -d '/' --output-delimiter " - "
    cat $Cov | cut -f3 | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
  done
done > analysis/genome_alignment/minimap/read_coverage.txt
```
```bash
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_*/*_depth.tsv); do
    Strain=$(echo $Cov | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Cov | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Cov)
    ProgDir=git_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_depth.tsv > $OutDir/${Organism}_${Strain}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_depth_10kb.tsv
  done
```


To visualise in IGV the new bam files we need the index of that bam file. Minimap did not generate the index files of the alignments, so I will try to index them using samtools.


```bash
for File in $(ls analysis/genome_alignment/minimap/*/*/*.fasta_aligned_sorted.bam); do
Strain=$(echo $File | rev | cut -f2 -d '/'| rev)
File=analysis/genome_alignment/minimap/*/$Strain/*.fasta_aligned_sorted.bam
samtools index <$File> -o $File
done
```

In order to extract the reads from the bam files we will use samtools:

```bash
#Run at cd analysis/genome_alignment/minimap/c.albicans/vs_ccXX. The -h option was added, because when sorting a message error appeared saying that the header was missing.

samtools view A22_Chromosomes_29.fasta_aligned_sorted.bam "Ca22chr1A_C_albicans_SC5314:8000-10000" -h > cc10_C1_region_8.bam

#Sorting of the bam file generated.

samtools sort -n cc10_C1_region_8.bam > cc10_C1_region_8_sorted.bam

#Now that we have the reads in a bam file, we can extract the fastq reads using bedtools.

bedtools bamtofastq -i cc10_C1_region_8_sorted.bam -fq cc10_C1_region_8_sorted.bam.fq

samtools view A22_Chromosomes_29.fasta_aligned_sorted.bam "Ca22chr4A_C_albicans_SC5314" | awk '$10 ~/ACTTCTTGGTGTACGGATGTCTA/' | awk '{print $1, $10}' > telseq.fa
```
