# bwsd: Burrows-Wheeler Similarity Distributions for String Collections

This code is an implementation of the algorithms described in \[[1](https://github.com/felipelouza/bwsd#references)\] to compute all-pairs of **Burrows-Wheeler similarity distributions (BWSD)** for a string collection.

Given a collection of _d_ strings, _bwsd_ outputs `option -o`:

* A matrix M<sub>dxd</sub> with all pairs of BWSD-based distances.

_Note:_ we compute only the (d<sup>2</sup>-d)/2 entries of M<sub>dxd</sub> (upper triangular matrix), which can be accessed as in [here](https://github.com/felipelouza/bwsd/blob/master/main.cpp#L312)

## install

```sh
git clone https://github.com/felipelouza/bwsd
cd bwsd
make compile
```

_Note:_ all algorithms need [sdsl-lite](https://github.com/simongog/sdsl-lite).

## run

Given a collection of _d_ strings in a single file FILE.

```sh
./bwsd [options] FILE d
```

Available options:

```sh
-h    this help message
-A a  preferred algorithm to use (default is alg. 1 BIT_sd)
-B b  BWSD-based distance to compute, options: 1. expectation (default), 2. shannon entropy
-o    write output matrix to FILE.output.bin
-p    print the output matrix (for debug)
-v    verbose output
```
_Notes:_ 
- supported extensions are _.txt_, _.fasta_ and _.fastq_.

## quick test

To run a test with d=10 strings from _dataset/input.100.txt_ using Alg. 1 `-A 1`, writing the output to disk `option -o`, type:

```sh
./bwsd dataset/input.100.txt 10 -A 1 -o
```
_Note:_ output matrix is written to `FILE.output.bin`

## options

To see the resulting matrix M<sub>dxd</sub>, use `option -p`:

```sh
./bwsd dataset/input.100.txt 10 -A 1 -o -p
```

```sh
## BWSD_BIT ##
0.00	0.00	-0.50	-0.67	-0.80	-0.80	-0.75	-0.80	-0.80	-0.67	
0.00	0.00	-0.50	-0.67	-0.80	-0.80	-0.75	-0.80	-0.80	-0.67	
-0.50	-0.50	0.00	-0.75	-0.75	-0.67	-0.75	-0.75	-0.75	-0.67	
-0.67	-0.67	-0.75	0.00	-0.80	-0.80	-0.80	-0.80	-0.67	-0.75	
-0.80	-0.80	-0.75	-0.80	0.00	-0.80	-0.67	-0.80	-0.83	-0.75	
-0.80	-0.80	-0.67	-0.80	-0.80	0.00	-0.75	-0.80	-0.80	-0.67	
-0.75	-0.75	-0.75	-0.80	-0.67	-0.75	0.00	-0.80	-0.83	-0.67	
-0.80	-0.80	-0.75	-0.80	-0.80	-0.80	-0.80	0.00	-0.75	-0.80	
-0.80	-0.80	-0.75	-0.67	-0.83	-0.80	-0.83	-0.75	0.00	-0.80	
-0.67	-0.67	-0.67	-0.75	-0.75	-0.67	-0.67	-0.80	-0.80	0.00
```

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
## external resources

We have included the source codes of the following algorithms:

* gSACA-K+DA: The Burrows-Wheeler transform (BWT) and the Document array (DA) were computed using [gsaca-k](https://github.com/felipelouza/gsa-is/).


## References

\[1\] 
Louza, F. A., & Telles, G. P. & Gog, S. & Zhao, L.: Computing Burrows-Wheeler Similarity Distributions for String Collections, 2018, To appear in Proc. SPIRE, 1-12. 

