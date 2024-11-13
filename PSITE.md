# Tumor Chain Simulation with PSiTE

This guide provides step-by-step instructions for simulating a tumor evolution chain using [PSiTE](https://github.com/hchyang/PSiTE). 
The simulation includes generating a coalescent tree, filtering variants, creating normal genomes, simulating tumor mutations, and assembling final genomes of tumor cells.

## Prerequisites

1. **Download Reference Files**:
   - **HG001_GRCh38_1_22_v4.2.1_all.vcf.gz** (Germline variants for the normal genome)
     - Available at: `[Download from NCBI](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/SupplementaryFiles/HG001_GRCh38_1_22_v4.2.1_all.vcf.gz)`
   - **COSMIC Non-Coding Variants** (Somatic mutation reference for cancer cells)
     - [Download from COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v100/noncodingvariantsvcf) (login required).

2. **Install PSiTE**:
   Ensure that `PSiTE` is correctly installed and accessible from the command line.

3. **Reference Genome**:
   - Path: `human_v38_no_chr.fasta`

Note: Your fasta should contain numerical values for chromosomes. If they don't, you can run the following script to remove the `chr` from the fasta:

```bash
python3 -c "import sys; [sys.stdout.write(line.replace('chr', '') if line.startswith('>') else line) for line in open(sys.argv[1])]" input.fasta > human_v38_no_chr.fasta
```

## Steps for Simulation

### 1. Generate a Coalescent Tree

Generate a coalescent tree of 10 tumor cells to model the evolutionary history of tumor cell mutations.

```bash
ms 10 1 -T -G 1 | tail -n1 > ms_tree.txt
```

- This uses the `ms` tool to simulate a simple tree structure of 10 tumor cells, saving it to `ms_tree.txt`.

### 2. Filter Germline Variants

Filter the germline VCF file to obtain only phased SNPs (germline variants) required to simulate the genome of a normal individual.

```bash
zcat HG001_GRCh38_1_22_v4.2.1_all.vcf.gz | awk '/^#/ || ($NF~/^[01]\|[01]/ && length($4)==1 && length($5)==1)' | gzip -c > HG001_GRCh38_1_22_v4.2.1_all.phased.vcf.gz
```

- This command extracts phased SNPs, filtering single nucleotide variants (SNVs) from the VCF file.

### 3. Simulate Normal Genome

Integrate germline variants into the reference genome to create a simulated normal genome of a female individual.

```bash
python3 psite.py vcf2fa -r human_v38_no_chr.fasta -v HG001_GRCh38_1_22_v4.2.1_all.phased.vcf.gz --autosomes 1..22 --sex_chr X,X -o normal_fa
```

- **Output**: `normal_fa` directory containing the simulated normal genome sequences.

### 4. Simulate Somatic Mutations in Tumor Chain

Using the coalescent tree and normal genome, simulate somatic (tumor) mutations.

```bash
python3 psite.py phylovar -t ms_tree.txt --config cfg_template_female.yaml --purity 0.8 --sex_chr X,X --trunk_length 2.0 --prune 0.05 --chain tumor_chain --map map
```

- **Output**: `tumor_chain`, representing the mutation history and tumor evolution.

### 5. Generate Tumor Cell Genomes

Build the full genomes of each tumor cell based on the simulated tumor chain.

```bash
python3 psite.py chain2fa -c tumor_chain -n normal_fa/normal.parental_0.fa,normal_fa/normal.parental_1.fa -o tumor_fa --cores 16
```

- **Output**: `tumor_fa` directory containing the genomes for each tumor cell.

### 6. Merge Diploid Files

Combine diploid parental files for each tumor node to create a single merged file for each cell.

```bash
nodes=$(ls ./tumor_fa | grep -oP 'node\d+' | sort | uniq); for node in $nodes; do parental_0="./tumor_fa/${node}.parental_0.fa"; parental_1="./tumor_fa/${node}.parental_1.fa"; if [[ -f $parental_0 && -f $parental_1 ]]; then cat $parental_0 $parental_1 > ./tumor_fa/${node}.fa; echo "Merged $parental_0 and $parental_1 into ./tumor_fa/${node}.fa"; fi done
```

- **Output**: Merged `.fa` files in `tumor_fa`, each representing a complete genome for a tumor cell.

### Single Command Execution

Alternatively, run the entire pipeline as a single command:

```bash
python3 psite.py vcf2fa -r human_v38_no_chr.fasta -v HG001_GRCh38_1_22_v4.2.1_all.phased.vcf.gz --autosomes 1..22 --sex_chr X,X -o normal_fa; python3 psite.py phylovar -t ms_tree.txt --config cfg_template_female.yaml --purity 0.8 --sex_chr X,X --trunk_length 2.0 --prune 0.05 --chain tumor_chain --map map; python3 psite.py chain2fa -c tumor_chain -n normal_fa/normal.parental_0.fa,normal_fa/normal.parental_1.fa -o tumor_fa --cores 16; nodes=$(ls ./tumor_fa | grep -oP 'node\d+' | sort | uniq); for node in $nodes; do parental_0="./tumor_fa/${node}.parental_0.fa"; parental_1="./tumor_fa/${node}.parental_1.fa"; if [[ -f $parental_0 && -f $parental_1 ]]; then cat $parental_0 $parental_1 > ./tumor_fa/${node}.fa; echo "Merged $parental_0 and $parental_1 into ./tumor_fa/${node}.fa"; fi done; mv tumor_fa fa; rm fa/*parental*
```

This command will execute the full pipeline, ending with merged tumor genome files in the `fa` directory.

## License

This project follows the licensing policies of the data sources and tools used.
