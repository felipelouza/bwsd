# bwsd 

This code is an implementation of three algorithms [1] to compute all-pairs of **Burrows-Wheeler similarity distributions (BWSD)** for a string collection.

The Burrows-Wheeler transform (BWT) and the Document array (DA) for string collections are computed using [gsaca-k](https://github.com/felipelouza/gsa-is/).

## install

```sh
git clone https://github.com/felipelouza/bwsd
cd bwsd
make compile
```

Note: all algorithms need [sdsl-lite](https://github.com/simongog/sdsl-lite).

## run

To run a test with K=5 strings from _dataset/input.100.txt_, type:

### Alg. 1

```sh
./bwsd dataset/input.100.txt 5 -M 1
```

### Alg. 2

```sh
./bwsd dataset/input.100.txt 5 -M 2
```

### Straightforward

```sh
./bwsd dataset/input.100.txt 5 -M 3
```

## options

Computing distance D\_M:

```sh
make OUTPUT=1
```

Computing distance D\_E:

```sh
make OUTPUT=2
```

To see the output matrix:

```sh
make DEBUG=1
./bwsd dataset/input.100.txt 5 -M 1
```


## alternatives (Alg. 1)

Using uncompressed bitvectors (bit\_vector):

```sh
make SD_VECTOR=0
./bwsd dataset/input.100.txt 5 -M 1
```

Using wavelet tree:

```sh
make WT=1
./bwsd dataset/input.100.txt 5 -M 1
```

Note: Alg. 1 uses compressed bitvectors (Bit\_sd) as default.

## References

\[1\] 
Louza, F. A., & Telles, G. P. & Gog, S. & Zhao, L.: Computing Burrows-Wheeler Similarity Distributions for String Collections, 2018, To appear in Proc. SPIRE, 1-12. 
