# Profiling with gprof

## Introduction
`gprof` is a performance analysis tool for profiling applications in Linux. It provides function call statistics and execution time details, helping to identify performance bottlenecks in a program.

## Compilation with Profiling Enabled
To use `gprof`, compile your program with the `-pg` flag:

```sh
 g++ -pg -o my_program my_program.cpp
```

## Running the Program
Execute the compiled binary as usual:

```sh
 ./my_program
```

This will generate a file named `gmon.out`, which contains profiling data.

## Generating the Profiling Report
Use `gprof` to analyze the `gmon.out` file and generate a report:

```sh
 gprof my_program gmon.out > profile.txt
```

## Understanding the Output
The `profile.txt` file contains two main sections:

1. **Flat Profile**: Lists the total execution time of each function.
2. **Call Graph**: Displays function call hierarchy and time spent in each function.

## Cleaning Up
After profiling, you can remove the generated `gmon.out` file to free up space:

```sh
 rm gmon.out
```

## Additional Options
- To get more detailed output:
  ```sh
  gprof -b my_program gmon.out
  ```
- To visualize the call graph using `gprof2dot` (if installed):
  ```sh
  gprof my_program gmon.out | gprof2dot | dot -Tpng -o call_graph.png
  ```

## Conclusion
Using `gprof` can help optimize performance by identifying bottlenecks in the program. By analyzing the generated reports, you can make informed decisions to improve efficiency.

