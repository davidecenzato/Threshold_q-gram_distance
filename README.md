Threshold q-gram distance
========
The *threshold* *q*-gram distance (```tqd```) measures the similarity between two sequences using the concept of *q*-grams, and is able to capture the hapax (uniquely occurring substring) and repeat content in the sequences.

Pre-requisites
--------------
The *threshold* *q*-gram distance requires:
* Python 3
* g++

Installation
--------------
```bash
git clone https://github.com/AlessioMilanese/Threshold_q-gram_distance.git
cd Threshold_q-gram_distance
./setup
```

Note: in the following examples we assume that the python script ```tqd``` is in the system path.


Simple examples
--------------
Here is a simple example on how to obtain the ```tqd``` between two fasta files:

```bash
tqd file1.fasta file2.fasta -q 6 -threshold 1
```

If you want to directly input the string sequences:
```bash
tqd ATGGATCAGTC CTGGATCAGAC -q 3 -threshold 1 -strings
```
