# Constructing a Variation Graph Using vg by Chunks

When constructing a variation graph using the `vg` tool, processing the entire genome [vg-issue](https://github.com/vgteam/vg/issues/4152) at once leads to crashes. To address this, the genome is split into smaller chunks, processed individually, and then merged back together. Below is the step-by-step guide to constructing the variation graph by chunks.

## Step 1: Splitting the Genome into Chunks

The genome is divided into chunks of 10Mb each to ensure manageable processing. The following script reads the FASTA index file (`hg38.fa.fai`) and generates chunked genomic regions:

```bash
CHUNK_SIZE=10000000; \
rm -f hg38.chunks.txt; \
while read -r chr chr_length _; do \
    for ((i = 0; i * CHUNK_SIZE < chr_length; i++)); do \
        start=$((i * CHUNK_SIZE + 1)); \
        end=$(((i + 1) * CHUNK_SIZE)); \
        [[ $end -gt $chr_length ]] && end=$chr_length; \
        echo "${chr}:${start}-${end}" >> hg38.chunks.txt; \
    done; \
done < ../hg38.fa.fai
```

This script iterates through each chromosome, dividing it into non-overlapping 10Mb chunks and saving them to `hg38.chunks.txt`.

## Step 2: Constructing Graph Chunks in Parallel

Using multiple threads, each chunk is processed separately to construct its respective variation graph. The following loop processes chunks in parallel using different thread counts (1, 2, 4, 8, and 16):

```bash
for t in 1 2 4 8 16; do \
    cores=$(seq -s, 0 $((t-1))); \
    i=0; \
    /bin/time -v taskset -c ${cores} bash -c "while read -r region; do output_file='chunk.t${t}.'\"\${i}\"'.vg'; vg construct -r ../hg38.fa -v ../pggb.vcf.gz -f -R \"\$region\" -t ${t} > \"\$output_file\"; ((i++)); done < hg38.chunks.txt" > hg38.pggb.vg.t${t}.out 2>&1; \
done
```

This script:
1. Iterates through different thread counts (`t` values: 1, 2, 4, 8, 16) to optimize performance.
2. Uses `taskset` to assign CPU cores for parallel processing.
3. Constructs each chunk separately using `vg construct`, saving outputs as `chunk.t${t}.*.vg`.
4. Logs runtime details for performance evaluation.

## Step 3: Merging the Processed Chunks

Once all chunks have been processed, they are combined into a single variation graph:

```bash
for t in 1 2 4 8 16; do \
    /bin/time -v vg combine -p chunk.t${t}.*.vg > hg38.pggb.t${t}.vg; rm chunk.t${t}.*.vg; \
done
```

This step:
1. Merges all chunked `.vg` files into a single `hg38.pggb.t${t}.vg` file using `vg combine`.
2. Removes temporary chunk files to free up disk space.

## Summary

This method efficiently constructs a variation graph by:
- Splitting the genome into 10Mb chunks.
- Processing each chunk in parallel using `vg construct`.
- Merging the processed chunks back into a single graph.

This approach helps overcome crashing failures. Yet, be aware that it fails to leverage parallelism.

