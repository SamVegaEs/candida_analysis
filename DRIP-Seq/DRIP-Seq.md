DATA EXTRACTION AND QC FROM NOVOGENE. 

1. Extract data

tar -xvf X204SC23101273-Z01-F001.tar

gzip -d *.gz

#Check Data integrity: 
md5sum -c MD5.txt

fastqc *.fq




2. Trim reads with trimgalore!

trim_galore --paired *_1.fq *_2.fq


1. First we check the quality of the data with fastqc 

fastqc *_val_*.fq


CHIP-Seq DATA analysis.

3. bowtie 


#Modify the bowtie file accordingly. To submit the job:

sbatch ~/git_repos/tools/seq_tools/genome_alignment/bowtie/bowtie2.sh



#bowtie gives the results in sam format. Next step is transforming from sam to bam. 


4. sam to bam

#Run the program with run_samtools_chip code the folder where the original sams are. Modify accordingly.

sbatch ~/git_repos/scripts/sbatch_scripts/run_samtools_chip.sh
