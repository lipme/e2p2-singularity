# E2P2 V3.1 Singularity Container

Singularity container for E2P2 version 3.1 (prediction of enzymatic functions).

## Version
v3.1

## Build
```bash
sudo bash
singularity build e2p2v3.sif Singularity
```

## Usage

### Display the options:
```bash
./e2p2v3.sif -h
# or
singularity exec e2p2v3.sif E2P2 -h
```
Use the -e option to set the blastp evalue threshold (default: 1e-2)

### Run on a multifasta file:
```bash
./e2p2v3.sif  -i proteome.fas -e "1e-6"
# or
singularity exec e2p2v3.sif E2P2   -i proteome.fas -e "1e-6"
```

## Maintainers
 - Sebastien.Carrere@inrae.fr
 - Ludovic.Cottret@inrae.fr


