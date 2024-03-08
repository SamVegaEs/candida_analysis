DATA EXTRACTION AND QC FROM NOVOGENE. 

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


CHIP-Seq DATA analysis.

For this analysis I followed to different tutorials:

https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/09_data_visualization.md
https://github.com/macs3-project/MACS/issues/356

First one is good to learn and to follow programs for pipeline and see what to do, but second is better as it has spike in correction methods. 

1. Align reads to reference assembly (A22) with bowtie2 

#Run with slurm bash code. Modify the bowtie file accordingly. To submit the job:

sbatch ~/git_repos/tools/seq_tools/genome_alignment/bowtie/bowtie2.sh

Running options: 

```bash
bowtie2 -x REF -1 READS_1.fq.gz (or READS_1.fq) -2 READS_2.fq.gz (or READS_2.fq) -S alignment.sam 
```
example:

```bash
bowtie2 -x ~/assembly/misc_publications/c.albicans/A22_Chromosomes_29 -1 ~/raw_dna/novogene/r_loops/batch_2/X204SC23101275-Z01-F001/01.RawData/C55_1/*_val_1.fq.gz -2 ~/raw_dna/novogene/r_loops/batch_2/X204SC23101275-Z01-F001/01.RawData/C55_1/*_val_2.fq.gz -S ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab55/ab55_c_albicans.sam 
```
#bowtie gives the results in sam format. Next step is transforming from sam to bam. 

Copy the output file generated with the alignment rates for each sample. Needed later on for spike in normalization. 

2. sam to bam

#Run the program with run_samtools_chip code the folder where the original sams are. Modify accordingly.

sbatch ~/git_repos/scripts/sbatch_scripts/run_samtools_chip.sh

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
  
This can be done with different programs:

Tutorial one uses sambamba: 

sbatch ~/git_repos/scripts/sbatch_scripts/run_sambamba_chip.sh

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

Followed the one of tutorial 2 to keep the same pipeline than spike in correction. Also the output generated is .bed, which takes less memory: https://github.com/macs3-project/MACS/issues/356

Running options

```bash
macs2 filterdup -f BAMPE -i ALIGNMENT_sorted.bam --keep-dup=1 -o ALIGNMENT_sorted_filterdup.bed
```
example:

```bash
macs2 filterdup -f BAMPE -i ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab1048_in/ab1048_in_c_albicans_ecoli_sorted.bam --keep-dup=1 -o ~/analysis/genome_alignment/bowtie2/chip_seq/r_loops/rep2/ab1048_in/ab1048_in_c_albicans_ecoli_sorted_filterdup.bed
```

4. Run Macs2


