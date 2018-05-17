# all-bwsd 


# run

To run a test with K=5 strings from DIR=dataset INPUT=input.100.txt type:

```sh
make
make run DIR=dataset/ INPUT=input.100.txt K=5 
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

Using sparse bitvectors (sd\_vector):

```sh
make SD_VECTOR=1
```

Wavelet tree:

```sh
make WT=1
```

