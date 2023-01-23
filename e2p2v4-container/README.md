# E2P2 V4 Singularity Container

Singularity container for E2P2 version 4 (prediction of enzymatic functions).
E2P2 Source: 20221206 Master version <https://github.com/carnegie/E2P2>  (including PRIAM jar file)

## Version
v4 (20221206)

## Build
```bash
sudo bash
singularity build e2p2v4.sif Singularity
```

## Usage

### Display the options:
```bash
./e2p2v4.sif -help
# or
singularity exec e2p2v4.sif E2P2 -help
```
Note: You don't need the following parameters since everything is packaged in the container:
- -b
- -bb
- -j
- -ps
- -r
- -pp

Use the -be option to set the blastp evalue threshold (default: 1e-2)

### Run on a multifasta file:
```bash
./e2p2v4.sif  -i proteome.fas -be "1e-6"
# or
singularity exec e2p2v4.sif E2P2   -i proteome.fas -be "1e-6"
```

## Maintainers
 - Sebastien.Carrere@inrae.fr
 - Ludovic.Cottret@inrae.fr


