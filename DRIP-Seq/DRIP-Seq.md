# DRIP-SEQ ANALYSIS

# Data extraction and QC 

1. Extract data

```bash
tar -xvf *.tar
```

2. Decompress data

```   
gzip -d *.gz
```

3. Check Data integrity: 
```
md5sum -c MD5.txt
```

4. Trim reads with trimgalore!

```bash
trim_galore --paired *_1.fq *_2.fq
```

5. QC with fastqc

```bash
fastqc *_val_*.fq
```

6. Reads can then be re-compressed. Bowtie can use .gz files so trimmed reads can also be zipped. 

# Peak Calling

For this analysis I followed to different tutorials:

https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/09_data_visualization.md
https://github.com/macs3-project/MACS/issues/356

First one is good to learn and to follow programs for pipeline and see what to do, but second is better as it has spike in correction methods. 

1. Align reads to reference assembly (A22) with bowtie2 

Run with slurm bash code. Modify the bowtie file accordingly. To submit the job:

```bash
sbatch ~/git_repos/tools/seq_tools/genome_alignment/bowtie/bowtie2.sh
```
Running options: 

```bash
bowtie2 -x REF -1 READS_1.fq.gz (or READS_1.fq) -2 READS_2.fq.gz (or READS_2.fq) -S alignment.sam 
```
example:

```bash
bowtie2 -x ~/assembly/misc_publications/c.albicans/A22_Chromosomes_29 -1 ~/raw_dna/novogene/r_loops/batch_2/X204SC23101275-Z01-F001/01.RawData/C55_1/*_val_1.fq.gz -2 ~/raw_dna/novogene/r_loops/batch_2/X204SC23101275-Z01-F001/01.RawData/C55_1/*_val_2.fq.gz -S ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab55/ab55_c_albicans.sam 
```
Bowtie gives the results in sam format. Next step is transforming from sam to bam. 
IMPORTANT: Copy the output file generated with the alignment rates for each sample. Needed later on for spike in normalization. 

2. sam to bam

Run the program with run_samtools_chip code the folder where the original sams are. Modify accordingly.

```bash
sbatch ~/git_repos/scripts/sbatch_scripts/run_samtools_chip.sh
```

Running options:

```bash
samtools view -S -b ALIGNMENT.sam > ALIGNMENT.bam
samtools sort ALIGNMENT.bam -o ALIGNMENT_sorted.bam
samtools index -b ALIGNMENT_sorted.bam
```
example: 
```bash
samtools view -S -b ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep4/ab55_rh_nacl/ab55_rh_nacl_c_albicans.sam > ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep4/ab55_rh_nacl/ab55_rh_nacl_c_albicans.bam
samtools sort ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep4/ab55_rh_nacl/ab55_nacl_c_albicans.bam -o ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep4/ab55_rh_nacl/ab55_rh_nacl_c_albicans_sorted.bam
samtools index -b ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep4/ab55_rh_nacl/ab55_rh_nacl_c_albicans_sorted.bam
```

3. Remove duplicated reads
  
This can be done with different programs: Sambamba or Macs2. I used Macs2, following the one of tutorial 2 to keep the same pipeline than spike in correction. Also the output generated is .bed, which takes less memory: https://github.com/macs3-project/MACS/issues/356

Tutorial one uses sambamba: 

```bash
sbatch ~/git_repos/scripts/sbatch_scripts/run_sambamba_chip.sh
```

Running options:

```bash
sambamba view -h -t 2 -f bam \
-F "[XS] == null and not unmapped and not duplicate" \
ALIGNMENT_sorted.bam > ALIGNMENT_sorted_filtered.bam
```
example:

```bash
sambamba view -h -t 2 -f bam \
-F "[XS] == null and not unmapped and not duplicate" \
~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab55/ab55_c_albicans_sorted.bam > ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab55/ab55_c_albicans_sorted_filtered.bam
```

Tutorial two uses macs2: 

```bash
sbatch ~/git_repos/scripts/tools/seq_tools/macs2/macs2_filterdup.sh
```

Running options

```bash
macs2 filterdup -f BAMPE -i ALIGNMENT_sorted.bam --keep-dup=1 -o ALIGNMENT_sorted_filterdup.bed
```
example:

```bash
macs2 filterdup -f BAMPE -i ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab1048_in/ab1048_in_c_albicans_ecoli_sorted.bam --keep-dup=1 -o ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab1048_in/ab1048_in_c_albicans_ecoli_sorted_filterdup.bed
```

4. Run Macs2

The program can be run in default options but it will not correct the spike-in, to do so it needs to be run in different steps. All the steps are summarised in the github tutorial

```bash
sbatch ~/git_repos/scripts/tools/seq_tools/macs2/macs2_covertracks.sh
```

#Step3: Extend ChIP sample to get ChIP coverage track

#Step4: Build local bias track from control (input)
        Build the slocal background
        Build llocal background
        Combine and generate the maximum background noise
        Then, take the maximum then by comparing with d background
        Finally, combine with the genome wide background using bdgopt subcommand

#Step 5: scaling chip and control based on spike in. Normalization will be done as: Take the sample with the lowest number of mapped reads for the spike-in (minMap). Take the minMap and divide by the total number of mapped reads in the sample to compute a normliaztion factor. Example:
Number of reads of spike in in each of the samples:
SampleA: wc -l SampleA_sorted_filterdup.bed: 4578/38464 = 0.11
INPUT: wc -l INPUT_sorted_filterdup.bed: 4578/4578= 1

#STEP6: Compare ChIP and local lambda to get the scores in pvalue or qvalue
#STEP7: Call peaks on score track using a cutoff
#-c for qval: 1.301. The scores in the output from bdgcmp are in -log10 form, so if you need the cutoff as 0.05, the -log10 value is about 1.3
#-c for pval: 2
#-l is same -d as STEP2: 147
#-g: Maximum gap between stronger peaks, better to set it as the tag size. DEFAULT: 30

Example:

```bash
Step3: macs2 pileup -f BEDPE -B -i ALIGNMENT_sorted_filterdup.bed -o sampleA_sorted_filterdup.pileup.bdg
Step4: macs2 pileup -f BEDPE -i ALIGNMENT_sorted_filterdup.bed -B --extsize 73 -o INPUT_d_bg.bdg
       macs2 pileup -f BEDPE -i ALIGNMENT_sorted_filterdup.bed -B --extsize 500 -o INPUT_1k_bg.bdg
       macs2 bdgopt -i INPUT_1k_bg.bdg -m multiply -p 0.147 -o INPUT_1k_bg_norm.bdg
       macs2 pileup -f BEDPE -i ALIGNMENT_sorted_filterdup.bed -B --extsize 5000 -o INPUT_10k_bg.bdg
       macs2 bdgopt -i INPUT_10k_bg.bdg -m multiply -p 0.0147 -o INPUT_10k_bg_norm.bdg
       macs2 bdgcmp -m max -t INPUT_1k_bg_norm.bdg -c INPUT_10k_bg_norm.bdg -o INPUT_d_1k_10k_bg_norm.bdg
       macs2 bdgopt -i INPUT_d_1k_10k_bg_norm.bdg -m max -p 65.94 -o INPUT_local_bias_raw.bdg
Step5: macs2 bdgopt -i SampleA_sorted_filterdup.pileup.bdg -p 0.11 -o SampleA_sorted_filterdup_scale.pileup.bdg
       macs2 bdgopt -i INPUT_local_bias_raw.bdg -m multiply -p 1 -o INPUT_local_lambda.bdg
Step6: macs2 bdgcmp -t SampleA_sorted_filterdup_scale.pileup.bdg -c INPUT_local_lambda.bdg -m qpois -o SampleA_qvalue.bdg
       macs2 bdgcmp -t SampleA_sorted_filterdup_scale.pileup.bdg -c INPUT_local_lambda.bdg -m ppois -o SampleA_pvalue.bdg
Step7: macs2 bdgbroadcall -i SampleA_qvalue.bdg -c 1.301 -l 147 -g 30 -o SampleA_qval_peaks.bed
       macs2 bdgbroadcall -i SampleA_pvalue.bdg -c 2 -l 147 -g 30 -o SampleA_pval_peaks.bed
```
PEAKS CALLING: DONE


