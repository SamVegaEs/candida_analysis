#Analysis of assemblies re-arragements by TRE. 

Promer alignment of Assemblies. Assemblies constructed by Jordan at NIAB.

# 1. Alignment of our assemblies against the reference A22 C. albicans

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

# 2. Coverage of Nanopore reads over Assembly. For the alignment of ONT reads versus Nanopore assembly use the program minimap:

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

# 3. SNP Calling of strains.

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
This database of Candida albicans has both A and B genomes. 

3.2 Prepare the indexed genomes required for GATK. If not present the program will run but not generate the .vcf files.

Prepare genome reference indexes required by GATK if needed:

For SC5314_A: 

```bash
for Reference in $(ls assembly/misc_publications/c.albicans/A22_Chromosomes_29.fasta); do
OutName=$(echo $Reference | sed 's/.fasta/.dict/g')
OutDir=$(dirname $Reference)
ProgDir=~/local/bin/picard-2.27.3
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
done
```
For SC5314_B:

```bash
for Reference in $(ls assembly/misc_publications/c.albicans/C_albicans_A22/C_albicans_A22_B.fa); do
OutName=$(echo $Reference | sed 's/.fa/.dict/g')
OutDir=$(dirname $Reference)
ProgDir=~/local/bin/picard-2.27.3
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
done
```

3.3 SNP Calling of TRE strains.

Done as: AB55 vs C3, C5 and F2. 

First all A22_A and then all A22_B

3.3.1 Run bwa from scripts. Modify them within the script.

3.3.2 Run samtools from scripts. Modify them within the script.

3.3.3 Run gatk from scripts. Modify them within the script.

3.3.4 Run picard from scripts. Modify them within the script.

3.4 SNPs filtering

Run with the script in sbtacth_run. After vcf is generated filter SNPs:

--> SC5314_A vs C3

```bash
Vcf=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C3/bwagatk_C3.vcf)
vcftools=/home/sv264/local/bin/vcftools_0.1.13/bin
vcflib=/home/sv264/miniconda3/pkgs/vcflib-1.0.0_rc2-h56106d0_2/bin
mq=40
qual=30
dp=10
gq=30
na=1.00
removeindel=N
echo "count prefilter"
cat ${Vcf} | grep -v '#' | wc -l

export LD_LIBRARY_PATH=/home/sv264/miniconda3/pkgs/bzip2-1.0.8-h7b6447c_0/lib
export LD_LIBRARY_PATH=/home/sv264/miniconda3/lib
export PYTHONPATH=/usr/lib64/python2.7
$vcflib/vcffilter -f "QUAL > $qual & MQ > $mq" $Vcf \
| $vcflib/vcffilter -g "DP > $dp & GQ > $gq" > ${Vcf%.vcf}_qfiltered.vcf

echo "count qfilter"
cat ${Vcf%.vcf}_qfiltered.vcf | grep -v '#' | wc -l
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --remove-indels --recode --out ${Vcf%.vcf}_qfiltered_presence
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --keep-only-indels --recode --out ${Vcf%.vcf}_indels_qfiltered_presence
```
```
count prefilter
92679

count qfilter
90699

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 76957 out of a possible 90699 Sites
Run Time = 1.00 seconds

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 13692 out of a possible 90699 Sites
Run Time = 0.00 seconds
```
Collect VCF stats
General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="C3"
  VcfTools=/home/sv264/local/bin/vcftools_0.1.13/perl
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C3/*_qfiltered_presence*.vcf | grep -v 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C3/*_qfiltered_presence*.vcf | grep 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/_indels.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
Isolate="C3"
  for Vcf in $(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C3/*_qfiltered_presence*.vcf | grep -v 'indels'); do
      ProgDir=~/git_repos/scripts/popgen/snp
      python2.7 $ProgDir/similarity_percentage.py $Vcf
  done
```

Annotate VCF files

```bash
Organism="Calbicans"
Isolate="SC5314"
Strain="C3"
DbName="Calbicans" 
CurDir=/home/sv264
cd $CurDir
  for Vcf in $(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_${Strain}/*_qfiltered_presence.recode.vcf | grep -v 'indels'); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $Vcf)
    SnpEff=/home/sv264/local/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 ${DbName} $Vcf > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'intron_variant') || (ANN[*].EFFECT has 'splice_region_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sv264/local/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
    printf "$filename\t$AllSnps\t$GeneSnps\t$CdsSnps\t$NonsynSnps\t$SynSnps\n"
done
```
```bash
bwagatk_C3_qfiltered_presence.recode.vcf        76957   36462   35967   13280   22687

Warnings detected:
WARNING_TRANSCRIPT_INCOMPLETE   118
WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS 306
WARNING_TRANSCRIPT_NO_START_CODON       56
WARNING_TRANSCRIPT_NO_STOP_CODON        206
```

--> SC5314_A vs C5

```bash
Vcf=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C5/bwagatk_C5.vcf)
vcftools=/home/sv264/local/bin/vcftools_0.1.13/bin
vcflib=/home/sv264/miniconda3/pkgs/vcflib-1.0.0_rc2-h56106d0_2/bin
mq=40
qual=30
dp=10
gq=30
na=1.00
removeindel=N
echo "count prefilter"
cat ${Vcf} | grep -v '#' | wc -l

export LD_LIBRARY_PATH=/home/sv264/miniconda3/pkgs/bzip2-1.0.8-h7b6447c_0/lib
export LD_LIBRARY_PATH=/home/sv264/miniconda3/lib
export PYTHONPATH=/usr/lib64/python2.7
$vcflib/vcffilter -f "QUAL > $qual & MQ > $mq" $Vcf \
| $vcflib/vcffilter -g "DP > $dp & GQ > $gq" > ${Vcf%.vcf}_qfiltered.vcf

echo "count qfilter"
cat ${Vcf%.vcf}_qfiltered.vcf | grep -v '#' | wc -l
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --remove-indels --recode --out ${Vcf%.vcf}_qfiltered_presence
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --keep-only-indels --recode --out ${Vcf%.vcf}_indels_qfiltered_presence
```
```
count prefilter
82913

count qfilter
81003

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 68344 out of a possible 81003 Sites
Run Time = 2.00 seconds

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 12557 out of a possible 81003 Sites
Run Time = 0.00 seconds

```
Collect VCF stats
General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="C5"
  VcfTools=/home/sv264/local/bin/vcftools_0.1.13/perl
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C5/*_qfiltered_presence*.vcf | grep -v 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C5/*_qfiltered_presence*.vcf | grep 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/_indels.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
Isolate="C5"
  for Vcf in $(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_C5/*_qfiltered_presence*.vcf | grep -v 'indels'); do
      ProgDir=~/git_repos/scripts/popgen/snp
      python2.7 $ProgDir/similarity_percentage.py $Vcf
  done
```
Annotate VCF files

```bash
Organism="Calbicans"
Isolate="SC5314"
Strain="C5"
DbName="Calbicans" 
CurDir=/home/sv264
cd $CurDir
  for Vcf in $(ls ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_${Strain}/*_qfiltered_presence.recode.vcf | grep -v 'indels'); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $Vcf)
    SnpEff=/home/sv264/local/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 ${DbName} $Vcf > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'intron_variant') || (ANN[*].EFFECT has 'splice_region_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sv264/local/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
    printf "$filename\t$AllSnps\t$GeneSnps\t$CdsSnps\t$NonsynSnps\t$SynSnps\n"
done
```
```bash
bwagatk_C5_qfiltered_presence.recode.vcf        68344   31940   31485   11749   19736
Warnings detected:
WARNING_TRANSCRIPT_INCOMPLETE   105
WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS 283
WARNING_TRANSCRIPT_NO_START_CODON       42
WARNING_TRANSCRIPT_NO_STOP_CODON        187
```
--> SC5314_A vs F2

```bash
Vcf=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_F2/bwagatk_F2.vcf)
vcftools=/home/sv264/local/bin/vcftools_0.1.13/bin
vcflib=/home/sv264/miniconda3/pkgs/vcflib-1.0.0_rc2-h56106d0_2/bin
mq=40
qual=30
dp=10
gq=30
na=1.00
removeindel=N
echo "count prefilter"
cat ${Vcf} | grep -v '#' | wc -l

export LD_LIBRARY_PATH=/home/sv264/miniconda3/pkgs/bzip2-1.0.8-h7b6447c_0/lib
export LD_LIBRARY_PATH=/home/sv264/miniconda3/lib
export PYTHONPATH=/usr/lib64/python2.7
$vcflib/vcffilter -f "QUAL > $qual & MQ > $mq" $Vcf \
| $vcflib/vcffilter -g "DP > $dp & GQ > $gq" > ${Vcf%.vcf}_qfiltered.vcf

echo "count qfilter"
cat ${Vcf%.vcf}_qfiltered.vcf | grep -v '#' | wc -l
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --remove-indels --recode --out ${Vcf%.vcf}_qfiltered_presence
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --keep-only-indels --recode --out ${Vcf%.vcf}_indels_qfiltered_presence
```
```
count prefilter
92473

count qfilter
90463

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 76679 out of a possible 90463 Sites
Run Time = 2.00 seconds

After filtering, kept 1 out of 1 Individuals
Outputting VCF file...
After filtering, kept 13642 out of a possible 90463 Sites
Run Time = 1.00 seconds 
```
Collect VCF stats
General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="F2"
  VcfTools=/home/sv264/local/bin/vcftools_0.1.13/perl
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_F2/*_qfiltered_presence*.vcf | grep -v 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
  VcfFiltered=$(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_F2/*_qfiltered_presence*.vcf | grep 'indels')
  Stats=$(echo $VcfFiltered | sed 's/.vcf/_indels.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
Isolate="F2"
  for Vcf in $(ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_F2/*_qfiltered_presence*.vcf | grep -v 'indels'); do
      ProgDir=~/git_repos/scripts/popgen/snp
      python2.7 $ProgDir/similarity_percentage.py $Vcf
  done
```
Annotate VCF files

```bash
Organism="Calbicans"
Isolate="SC5314"
Strain="F2"
DbName="Calbicans" 
CurDir=/home/sv264
cd $CurDir
  for Vcf in $(ls ls analysis/popgen/SNP_calling/c.albicans/SC5314/vs_${Strain}/*_qfiltered_presence.recode.vcf | grep -v 'indels'); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $Vcf)
    SnpEff=/home/sv264/local/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 ${DbName} $Vcf > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant') || (ANN[0].EFFECT has 'intron_variant') || (ANN[*].EFFECT has 'splice_region_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'frameshift_variant') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'disruptive_inframe_deletion') || (ANN[0].EFFECT has 'disruptive_inframe_insertion') || (ANN[0].EFFECT has 'exon_loss_variant') || (ANN[0].EFFECT has 'splice_donor_variant') || (ANN[0].EFFECT has 'splice_acceptor_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/home/sv264/local/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
    printf "$filename\t$AllSnps\t$GeneSnps\t$CdsSnps\t$NonsynSnps\t$SynSnps\n"
done
```
```bash
bwagatk_F2_qfiltered_presence.recode.vcf        76679   36323   35827   13227   22600

Warnings detected:
WARNING_TRANSCRIPT_INCOMPLETE   118
WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS 295
WARNING_TRANSCRIPT_NO_START_CODON       55
WARNING_TRANSCRIPT_NO_STOP_CODON        207
```







