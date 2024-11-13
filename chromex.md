# Chromosome Extractor Script

This script extracts specified chromosomes from a FASTA file directly from the command line. 
It is designed to be compact and is intended to be run as a single-line command.

## Prerequisites
- Python 3
- [Biopython](https://biopython.org/) library (can be installed with `pip install biopython`)

## Script Description

The script uses Biopython's `SeqIO` module to parse a FASTA file and filter out sequences that match specific chromosome IDs. 
It outputs the selected chromosomes in FASTA format directly to the terminal or to a specified output file.

## Usage

To run this script, use the following command in the terminal:

```bash
python3 -c "import sys; from Bio import SeqIO; [print(f'>{rec.id}\n{rec.seq}') for rec in SeqIO.parse(sys.argv[1], 'fasta') if rec.id in sys.argv[2:]]" input.fasta chr1 chr2 chrX > output.fasta
```

### Parameters
- `input.fasta`: Path to the input FASTA file containing all chromosomes.
- `chr1 chr2 chrX`: Space-separated list of chromosome names (headers) you want to extract.

### Notes
- The output is in FASTA format.
- You can replace `output.fasta` with any file name you choose to save the result.

## License
This script is provided under the MIT License.
