# Profiling with perf

## Introduction
`perf` is a powerful Linux profiling tool used to analyze CPU performance, I/O statistics, cache usage, and memory access patterns. It is particularly useful for understanding the behavior of parallelized programs.

## Installing perf
On Debian-based systems (Ubuntu):
```sh
sudo apt-get install linux-tools-common linux-tools-generic linux-tools-`uname -r`
```
On RHEL-based systems:
```sh
sudo yum install perf
```

## Basic Usage
To record performance data while running a program:
```sh
perf record -o perf.data ./my_program
```
To generate a summary report:
```sh
perf report -i perf.data
```

## Profiling I/O Performance
To measure disk I/O statistics:
```sh
perf stat -e block:block_rq_issue,block:block_rq_complete ./my_program
```
This provides insights into read/write requests and their completion times.

## Cache Usage Statistics
To measure different levels of cache performance:
```sh
perf stat -e L1-dcache-loads,L1-dcache-load-misses ./my_program
perf stat -e L2_rqsts.all,L2_rqsts.miss ./my_program
perf stat -e LLC-loads,LLC-load-misses ./my_program
```
- **L1-dcache-loads**: Number of L1 data cache loads.
- **L1-dcache-load-misses**: L1 cache load misses.
- **L2_rqsts.all**: All L2 cache requests.
- **L2_rqsts.miss**: L2 cache misses.
- **LLC-loads**: Last Level Cache (LLC) loads.
- **LLC-load-misses**: LLC load misses.

## RAM Access Statistics
To analyze RAM accesses:
```sh
perf stat -e dTLB-loads,dTLB-load-misses ./my_program
```
- **dTLB-loads**: Data Translation Lookaside Buffer (dTLB) loads.
- **dTLB-load-misses**: TLB misses resulting in memory access.

## Parallelization Analysis
To analyze CPU utilization and thread performance:
```sh
perf stat -e task-clock,context-switches,cpu-migrations,page-faults ./my_program
```
- **task-clock**: Total time spent on CPU.
- **context-switches**: Number of context switches.
- **cpu-migrations**: Number of times a process migrated between CPUs.
- **page-faults**: Number of memory page faults.

To track thread execution:
```sh
perf record -g ./my_program
perf report --stdio
```
The `-g` flag records call graphs, helping to analyze thread behavior.

## Cleaning Up
To remove profiling data files:
```sh
rm perf.data
```

## Conclusion
Using `perf`, you can gain deep insights into I/O performance, cache efficiency, RAM access, and parallelization bottlenecks, enabling efficient performance tuning for your programs.

