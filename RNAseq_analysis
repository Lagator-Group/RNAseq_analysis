## Quality control with trim-galore

# Install trim-galore to your conda environment
conda install -c bioconda trim-galore

# For paired end reads, give both R1 and R2 files and specify with option -paired
# -q 20 - quality score cut off 
# --phred33 - starts after 33 characters
# --fastqc - does quality check after trimming
# --length 50 - removes reads under 50bp long
# -o - output directory 
# -j 8 - uses 8 CPUs (multi threading)
trim_galore --paired yourfilename_R1.fastq.gz yourfilename_R2.fastq.gz -q 20 --phred33 --fastqc --length 50 -o fastq_trimmed/ -j 8 


## Remove ribosomal RNA that survived depletion with BBDuk

# Install bbmap to conda environment
conda install -c bioconda bbmap

# Download ribokmers.fa.gz file from this post - https://www.biostars.org/p/159959/
# Add to directory and decompress with gunzip

# Make directories to store removed rRNA reads and the reads you will use without rRNA reads
mkdir fastQ_trimmed_norRNA fastQ_trimmed_rRNA

# Use bbduk.sh to remove rRNA reads and spit out final reads that you will use for your alignment 
# If you specify in= and in2= bbduk will assume paired ends
# k=31 - will use kmers of length 31
# ref= - pass the file downloaded above here
bbduk.sh in=fastq_trimmed/yourfilename_R1.fq.gz in2=fastq_trimmed/yourfilename_R2.fq.gz out=fastQ_trimmed_norRNA/yourfilename_R1.fastq.gz out2=fastQ_trimmed_norRNA/yourfilename_R2.fastq.gz outm=fastQ_trimmed_rRNA/yourfilename_R1.fastq.gz outm2=fastQ_trimmed_rRNA/yourfilename_R2.fastq.gz k=31 ref=ribokmers.fa

# To QC, can see how many rRNA reads were removed and retained using seqkit 
conda install -c bioconda seqkit
seqkit stats fastQ_trimmed_norRNA/*.fastq.gz
seqkit stats fastQ_trimmed_rRNA/*.fastq.gz


## Align mRNA reads against reference genome

# Install bowtie2 to your conda environment
conda install -c bioconda bowtie2

# Download your reference genome from NCBI (FASTA and GTF files) and add to your directory
# (GTF is an annotation table that we will use to transform the SAM alignment info into the table which counts how many times a read aligned to each gene)

# Build the Bowtie2 database from FASTA file (in this example, BW25113)
# --threads 10 - specify number of threads for multithreading
bowtie2-build --threads 10 reference_genome.fna reference_genome

# Create directory to store SAM files (Sequence Alignment Map)
mkdir Bowtie2_SAM

# Align reads (NB this is the most time intensive process)
 bowtie2 -x reference_genome -1 fastQ_trimmed_norRNA/yourfilename_R1.fastq.gz -2 fastQ_trimmed_norRNA/yourfilename_R2.fastq.gz -S Bowtie2_SAM/yourfilename.sam --no-mixed --threads 10


## Convert SAM to BAM files

# Make directory to store BAM files 
mkdir BAM

# Install samtools to your conda environment
conda install -c bioconda samtools

# Next convert the SAM files in to BAM (same thing, just binary encoded so saves space)
# -bS means BAM output and SAM input
# -@ 10 - number of threads
samtools view -bS Bowtie2_SAM/yourfilename.sam -@ 10 > BAM/yourfilename.bam

## Sort and index BAM files
# -@ 10 - number of threads
# -o - output directory
samtools sort BAM/yourfilename.bam -@ 10 -o BAM/yourfilename_sorted.bam 

# The indexing should create a bunch of files with the suffix .bai - no need to specify outputs or delete anything
samtools index BAM/yourfilename_sorted.bam -@ 10

# You can delete the unsorted BAM files if the sorting worked 
# To confirm, look at filesize with 
ls -lh
# Should be able to tell if it worked if the filesize isn't zero


## Generate count table

# Convert the alignment info into the count table with featureCounts (in the conda package "subreads")
conda install -c bioconda subread

# -p - specifies that fragments (or templates) will be counted instead of reads. This is only applicable for paired-end reads.
# -T 5 - specifies number of threads
# -t CDS - specifies the feature type
# -o - output file
featureCounts -a reference_genome.gtf -p -T 5 -t CDS -o featureCounts_table.txt yourfilename.bam
