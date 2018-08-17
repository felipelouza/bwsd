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
- If the input is changed, please run `make remove DIR=dataset/`, to rebuild the BWTs.

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
0.000	0.000	0.684	1.459	2.156	2.128	1.750	2.156	1.896	1.299
0.000	0.000	0.684	1.459	2.156	2.128	1.750	2.156	1.896	1.299
0.684	0.684	0.000	1.614	1.811	1.224	1.750	1.750	1.614	1.299
1.459	1.459	1.614	0.000	2.250	2.059	2.156	2.156	0.959	1.750
2.156	2.156	1.811	2.250	0.000	2.322	0.944	2.020	2.522	1.686
2.128	2.128	1.224	2.059	2.322	0.000	1.792	2.252	1.961	1.459
1.750	1.750	1.750	2.156	0.944	1.792	0.000	1.627	2.522	1.392
2.156	2.156	1.750	2.156	2.020	2.252	1.627	0.000	1.753	1.506
1.896	1.896	1.614	0.959	2.522	1.961	2.522	1.753	0.000	2.000
1.299	1.299	1.299	1.750	1.686	1.459	1.392	1.506	2.000	0.000
```

_Notes:_ 
- We compute only the (d<sup>2</sup>-d)/2 entries of M<sub>dxd</sub> (upper triangular matrix), which can be accessed as [here](https://github.com/felipelouza/bwsd/blob/master/main.cpp#L312)


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

