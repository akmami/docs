# High-Quality HiFi Read Generation with PBSIM3

This guide outlines the process for simulating high-fidelity (HiFi) reads using PBSIM3, converting the output to BAM format, sorting and merging the files, and generating FASTQ files for further analysis. 
This process simulates long-read sequencing using the `pbsim3` tool, processes SAM/BAM files, and generates HiFi reads for downstream applications.

## Prerequisites

Ensure the following tools are installed and available:

1. **PBSIM3**: A tool for simulating PacBio HiFi reads. [PBSIM3 GitHub Repository](https://github.com/yukiteruono/pbsim3)

2. **Samtools**: A suite of tools for handling SAM/BAM files. [Samtools GitHub Repository](https://github.com/samtools/samtools)

3. **CCS (Circular Consensus Sequencing)**: Tool for generating consensus reads from BAM files. Installation instructions are provided below. [CCS Tool](https://github.com/PacificBiosciences/ccs)

4. **pigz**: A parallel implementation of gzip, used for compressing FASTQ files. [pigz GitHub Repository](https://github.com/madler/pigz)

Make sure all paths are correctly set, and the required reference files are in place.

## Installing CCS

CCS can be installed directly using conda:

1. Install Miniconda or Anaconda if it is not already installed. Anaconda Installation Guide

2. Install CCS via conda without creating a new environment:

```bash
conda install bioconda::pbccs
```
3. Verify the installation by checking the CCS version:

```bash
ccs --version
```
    
## Steps for HiFi Read Generation

### 1. Simulate High-Fidelity (HiFi) Reads using PBSIM3

First, generate HiFi reads with `pbsim3` using the specified genome file and parameters.

```bash
pbsim \
  --strategy wgs \
  --method qshmm \
  --qshmm pbsim3/data/QSHMM-RSII.model \
  --depth 15 \
  --genome hg38.fa \
  --pass-num 7
```

Parameters:
- `--strategy wgs`: Specifies the Whole Genome Sequencing (WGS) strategy for simulating the reads.
- `--method qshmm`: Uses the QSHMM method for read simulation.
- `--qshmm`: Path to the QSHMM model file (should be inside the pbsim3 directory).
- `--depth 15`: The sequencing depth for simulating the reads.
- `--genome`: Path to the reference genome file (hg38.fa in this case).
- `--pass-num 7`: Specifies the number of passes for the HiFi reads.

### 2. Convert SAM Files to BAM Files

After generating the reads, convert the generated SAM files to BAM format using samtools. The following loop processes all SAM files (`sd_0001.sam`, `sd_0002.sam`, ..., `sd_0046.sam`):

```bash
for i in $(seq -f "%04g" 1 46); do 
  input_file="sd_$i.sam"
  output_file="sd_$i.bam"
  if [ -f "$input_file" ]; then 
    echo "Processing $input_file to $output_file"
    samtools view -b -o "$output_file" "$input_file"
  fi
done
```

Explanation:

- Loops through all `sd_XXXX.sam` files.
- Uses samtools view to convert each SAM file to its corresponding BAM file.

### 3. Merge Sorted BAM Files

Merge all the BAM files into a single file:

```bash
samtools merge sd_merged.bam sd_*.bam
```
### 4. Sort the merged BAM File

Next, sort the BAM file to prepare them for downstream processing:

```bash
samtools sort -o sd_merged_sorted.bam sd_merged.bam
```

### 5. Generate Circular Consensus Sequences (CCS)

Ensure CCS is installed as described above. Then generate consensus sequences from the merged BAM file:

```bash
ccs sd_merged_sorted.bam sd_merged_sorted.ccs.bam
```

### 6. Convert the CCS BAM to FASTQ Format

Convert the CCS BAM file into a FASTQ file using samtools:

```bash
samtools bam2fq sd_merged_sorted.ccs.bam > hg38.fastq
```

### 7. Compress the FASTQ File

Finally, compress the generated FASTQ file using `pigz` to speed up compression with multiple threads:

```bash
pigz -p 8 hg38.fastq
```

## Single Command Execution

```bash
prefix=hg38; ref=${prefix}.fa; fq=${prefix}.fq; count=46; pbsim --strategy wgs --method qshmm --qshmm pbsim3/data/QSHMM-RSII.model --depth 15 --genome ${ref} --pass-num 7 --prefix ${prefix}; rm -f *.ref; for i in $(seq -f "%04g" 1 ${count}); do if [ -f "${prefix}_${i}.sam" ]; then echo "Processing ${prefix}_${i}.sam"; samtools view -b -o ${prefix}_${i}.bam ${prefix}_${i}.sam; rm -f ${prefix}_${i}.sam; fi done; samtools merge ${prefix}_merged.bam ${prefix}_*.bam; for i in $(seq -f "%04g" 1 ${count}); do rm -f "${prefix}_${i}.bam"; done; samtools sort -o ${prefix}_merged_sorted.bam ${prefix}_merged.bam; rm -f ${prefix}_merged.bam; ccs ${prefix}_merged_sorted.bam ${prefix}_merged_sorted.ccs.bam; rm -f ${prefix}_merged_sorted.bam; samtools bam2fq ${prefix}_merged_sorted.ccs.bam > ${fq}; pigz -p 8 ${fq}; rm -f ${prefix}_merged_sorted.ccs.bam ${prefix}_merged_sorted.ccs.bam.pbi;
```

# License

This project follows the licensing policies of the data sources and tools used.
