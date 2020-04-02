# E2PE Singularity

Singularity instance for E2P2 (prediction of enzymatic functions).

E2P2 V4 forked by Ludo Cottret.

## Build

```
sudo bash
singularity build e2p2-singularity e2p2.singularity
```

## Usage

Display the options:

```
e2p2-singularity -help
```

Note: You don't need the following parameters since everything is packaged in the container:
- -b 
- -bb 
- -j 
- -ps 
- -r 
- -pp 

Use the -be option to set the blastp evalue threshold (default: 1e-2)


Run on a multifasta file:
```
e2p2-singularity -i proteome.fas -be "1e-6"
```

