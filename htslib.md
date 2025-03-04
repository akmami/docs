# Using `HTSLIB` Library

## Table of Contents
- [Using TBI Index](#using-tbi-index)
- [Working with BGZF](#working-with-bgzf)
- [Handling BAM Files](#handling-bam-files)

---

## Using TBI Index

### Overview
This document provides an example of how to use HTSlib to query a VCF file (.vcf.gz) using its TBI index. The provided C program demonstrates how to open a VCF file, load its index, and retrieve records within a specified genomic region.

### Prerequisites
Ensure you have the following installed:
- **HTSlib** (https://github.com/samtools/htslib)
- A **VCF file** compressed using `bgzip`
- A **TBI index** created using `tabix`

### Including HTSlib and Initializing
To use HTSlib, include the necessary headers:
```c
#include <stdio.h>
#include <stdlib.h>
#include "htslib/tbx.h"
```

### Opening the VCF File
To read the compressed VCF file:
```c
htsFile *fp = hts_open("sample.vcf.gz", "r");
if (!fp) {
    fprintf(stderr, "Error: Unable to open VCF file\n");
    return 1;
}
```

### Loading the TBI Index
To load the index file:
```c
tbx_t *tbx = tbx_index_load("sample.vcf.gz");
if (!tbx) {
    fprintf(stderr, "Error: Could not load index\n");
    hts_close(fp);
    return 1;
}
```

### Querying a Genomic Region
To search a region:
```c
hts_itr_t *itr = tbx_itr_queryi(tbx, chrom-1, start, end);
kstring_t str = {0, 0, NULL};
if (itr) {
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
        printf("%s\n", str.s);
    }
}
```

### Iterating Over the Results
Use `tbx_itr_next` to iterate through the results:
```c
while (tbx_itr_next(fp, tbx, itr, &str) >= 0) {
    printf("%s\n", str.s);
}
```

### Cleanup and Freeing Resources
After processing, clean up the allocated resources:
```c
tbx_itr_destroy(itr);
tbx_destroy(tbx);
hts_close(fp);
```

## Notes
- Use `bgzip` for compression before indexing:
  ```sh
  bgzip -c sample.vcf > sample.vcf.gz
  ```
- Ensure the `.vcf.gz` file is indexed with `tabix`:
  ```sh
  tabix -p vcf sample.vcf.gz
  ```
- Chromosomes are 1-based in input but 0-based internally in `tbx_itr_queryi`.

## References
- HTSlib Documentation: https://www.htslib.org/
- Tabix User Guide: http://www.htslib.org/doc/tabix.html

---
