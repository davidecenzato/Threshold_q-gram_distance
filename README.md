# Threshold q-gram distance

This tool implements the *threshold* *q*-gram distance (```tqd```), an alignment free distance measure on strings which is similar to the *q*-gram distance by Ukkonen but uses reduced information on the multiplicities of the *q*-grams.

### Requirements

* A modern C++17 compiler such as `g++` version 8.4 or higher.
* A Linux or macOS 64-bit operating system.
* CMake installed on your system.

### Download and Compile

To clone the repository, run the following commands:

```console
git clone https://github.com/davidecenzato/Threshold_q-gram_distance.git
cd Threshold_q-gram_distance
cmake ..
make
```

### Usage
Compute the pairwise distance matrices of the input sequences for several combinations of input parameters:

```console
tqd-dna [options]
Options:
-h          Print usage info.
-i <arg>    Directory path containing the input FASTA files. (REQUIRED)
-q <arg>    q-gram lengths list. (Def. 8)
-t <arg>    Threshold values list. (Def. 1)
-o <arg>    Output directory path for pairwise distance matrices. (REQUIRED)
 ```

  ```
./sources/tqd-dna -i input_data -q 7,8,9 -t 1,2,6 -o output_matrices
 ```
This command computes the pairwise distance matrices for all combinations of  
`q = {7, 8, 9}` and `t = {1, 2, 6}` using the input files in the `input_data` directory.  
The resulting TSV-formatted matrices are saved in the `output_matrices` folder.

> **Note:** Each FASTA file in the input directory is supposed to contain exactly one sequence,  
> or the set of reads corresponding to a single sequence.

### External resources

* [sdsl-lite](https://github.com/simongog/sdsl-lite.git)
* [robin-map](https://github.com/Tessil/robin-map)

### Authors

- [**Davide Cenzato**](https://github.com/davidecenzato)
- **Giuditta Franco**
- **Zsuzsanna Lipták**
- [**Alessio Milanese**](https://github.com/AlessioMilanese)

### Publication

**Davide Cenzato**, **Giuditta Franco**, **Zsuzsanna Lipták**, **Alessio Milanese**. The threshold *q*-gram distance: a simple, efficient, and effective distance measure for genomic sequence comparison. *Submitted*, 2025.

### Contacts

If you notice any bugs, please feel free to report them by opening a Git issue or by contacting us at davide[dot]cenzato[at]unive[dot]it.