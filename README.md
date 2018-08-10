# bwsd: 

This software is an implementation of the algorithms described in \[[1](https://github.com/felipelouza/bwsd#references)\] to compute all-pairs of **Burrows-Wheeler similarity distributions (BWSD)** for a string collection.

Given a collection of _d_ strings, _bwsd_ computes a matrix M<sub>dxd</sub> with all pairs of BWSD-based distances.

The Burrows-Wheeler transform (BWT) and the Document array (DA) are computed using [gsaca-k](https://github.com/felipelouza/gsa-is/) \[2\].

## install

```sh
git clone https://github.com/felipelouza/bwsd
cd bwsd
make compile
```

_Note:_ our implementaion needs [sdsl-lite](https://github.com/simongog/sdsl-lite).

## run

Given a collection of _d_ strings in a single file FILE.

```sh
./bwsd [options] FILE d
```

Available options:

```sh
-h    this help message
-A a  preferred algorithm to use (default is alg. 1 BIT_sd)
-B b  BWSD-based distance to compute, options: 1. expectation (default), 2. Shannon entropy
-T t  use t parallel threads (default 1)
-o    write output matrix to FILE.output.bin
-p    print the output matrix (for debug)
-v    verbose output
```
_Notes:_ 
- Supported extensions are _.txt_, _.fasta_ and _.fastq_.

## quick test

To run a test with d=10 strings from [dataset/input.100.txt](https://github.com/felipelouza/bwsd/blob/master/dataset/input.100.txt) using Alg. 1 `-A 1`, writing the output to disk `option -o`, type:

```sh
./bwsd dataset/input.100.txt 10 -A 1 -o
```
_Note:_ output matrix is written to `FILE.output.bin`

## options

To see the resulting M<sub>dxd</sub>, computing distances based on the Shannon entropy of BWSD `option -B 2`, use `option -p`:

```sh
./bwsd dataset/input.100.txt 10 -A 1 -o -B 2 -p
```

Result:

```sh
## BWSD_BIT_sd ##
writing 360 bytes to: input.100.txt.output.bin
0.00	0.00	0.68	1.46	2.16	2.13	1.75	2.16	1.90	1.30	
0.00	0.00	0.68	1.46	2.16	2.13	1.75	2.16	1.90	1.30	
0.68	0.68	0.00	1.61	1.81	1.22	1.75	1.75	1.61	1.30	
1.46	1.46	1.61	0.00	2.25	2.06	2.16	2.16	0.96	1.75	
2.16	2.16	1.81	2.25	0.00	2.32	0.94	2.02	2.52	1.69	
2.13	2.13	1.22	2.06	2.32	0.00	1.79	2.25	1.96	1.46	
1.75	1.75	1.75	2.16	0.94	1.79	0.00	1.63	2.52	1.39	
2.16	2.16	1.75	2.16	2.02	2.25	1.63	0.00	1.75	1.51	
1.90	1.90	1.61	0.96	2.52	1.96	2.52	1.75	0.00	2.00	
1.30	1.30	1.30	1.75	1.69	1.46	1.39	1.51	2.00	0.00	
```

_Notes:_ 
- We compute only the (d<sup>2</sup>-d)/2 entries of M<sub>dxd</sub> (upper triangular matrix), which can be accessed as in [here](https://github.com/felipelouza/bwsd/blob/master/main.cpp#L312)


## alternatives (Alg. 1)

Using uncompressed bitvectors (bit\_vector):

```sh
make SD_VECTOR=0
./bwsd [options] FILE d -A 1
```

Using wavelet tree:

```sh
make WT=1
./bwsd [options] FILE d -A 1
```

_Note:_ Alg. 1 uses compressed bitvectors (BIT\_sd) as default.

## debug

To see a more detailed execution use:

```sh
make clean
make DEBUG=1
```
## citation

Please, if you use this tool in an academic setting cite the following paper \[1\]:

    @inproceedings{LouzaTGZ18,
      author    = {Louza, Felipe A. and Telles, Guilherme P. and Gog, Simon and Zhao, Liang},
      title     = {Computing Burrows-Wheeler Similarity Distributions for String Collections},
      booktitle = {String Processing and Information Retrieval - 25th International Symposium,
                  {SPIRE} 2018, Proceedings},
      pages     = {1--12},
      year      = {2018},
      series    = {Lecture Notes in Computer Science},
      volume    = {},
      publisher = {Springer}
    }

## references

\[1\] 
Louza, Felipe A., Telles, Guilherme P., Gog, Simon, Liang Zhao (2018): 
Computing Burrows-Wheeler Similarity Distributions for String Collections. 
To appear in Proc. SPIRE 2018: 1-12.

\[2\] 
Louza, Felipe A., Gog, Simon, Telles, Guilherme P. (2017). 
Inducing enhanced suffix arrays for string collections. 
Theor. Comput. Sci. 678: 22-39, [github](https://github.com/felipelouza/gsa-is).

