# bwsd 

This code is an implementation of three algorithms [1] that compute all-pairs of Burrows-Wheeler similarity distributions (BWSD) for a string collection.

# run

To run a test with K=5 strings from DIR=dataset INPUT=input.100.txt type:

```sh
make compile
```

Note: all algorithms need [sdsl-lite](https://github.com/simongog/sdsl-lite).

### Alg. 1

```sh
make run DIR=dataset/ INPUT=input.100.txt K=5 MODE=1
```

### Alg. 2

```sh
make run DIR=dataset/ INPUT=input.100.txt K=5 MODE=2
```

### Straightforward

```sh
make run DIR=dataset/ INPUT=input.100.txt K=5 MODE=3
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
make run DIR=dataset/ INPUT=input.100.txt K=5 MODE=3
```


## alternatives (Alg. 1)

Using sparse bitvectors (sd\_vector):

```sh
make SD_VECTOR=1
make run DIR=dataset/ INPUT=input.100.txt K=5 MODE=1
```

Using wavelet tree:

```sh
make WT=1
make run DIR=dataset/ INPUT=input.100.txt K=5 MODE=1
```

## References

\[1\] 
Louza, F. A., & Telles, G. P. & Gog, S. & Zhao, L.: Computing Burrows-Wheeler Similarity Distributions for String Collections, 2018, To appear in Proc. SPIRE, 1-12. 
