# toRNAdo

This repository provides functionalities for measuring the distance of a given RNA folding to the set of Nussinov-optimal or Zuker-optimal foldings of the RNA string. Specifically, it provides a command line tool that allows the user to specify an RNA string, a folding for that RNA string (in .dbn format), and an algorithm for predicting the optimal RNA foldings (either Nussinov or Zuker's algorithm), and then it outputs the maximum distance between the given folding and any optimal folding as well as a distance vector that tells the number of optimal foldings there are that are some distance away from the given folding. The following section explains how to use the command line tool in more detail.

## Set Up
It is advisable for the user to set up a virtual environment to install the correct packages. We provide instructions for setting up a virtual environment using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html):
- Create a new Conda environment with python version: `conda create -n newenv python=3.9`
- Activate Conda environment: `conda activate newenv`
- Install matplotlib for drawing the distance vector as a histogram: `conda install matplotlib`

## Command Line Tool

To display help messages, run:
```bash
python tornado.py --help
```
The required arguments are:
- `rna_file`: filepath for the RNA string
- `{N,Z}`: select an algorithm to get optimal foldings: N for Nussinov, Z for Zuker

The optional arguments are:
- `-dbn <DBN_FILE>`, `--dbn_file <DBN_FILE>`: Specify the filepath where the dbn string is stored. The dbn string is the dot-bracket notation that represents the folding we want to compare against.
- `-s <INTERNAL_LOOP_SIZE>`, `--internal_loop_size <INTERNAL_LOOP_SIZE>`: Specify a positive integer that bounds the maximum size of an internal loop in Zuker's algorithm for speeding up the runtime. If unspecified, there is no constraint on the internal loop size.
- `-d`, `--display_max_distance`: Display a folding of maximum distance from the given folding
- `-p`, `--plot`: Plot the histogram for the distance vector
- `-f <PLOT_FILE>`, `--plot_file <PLOT_FILE>`: Specify the filepath for storing the histogram plot of the distance vector

For example, if we want to find the maximum distance and distance vector for a given folding stored in `examples/dbn_string.txt` and compare it with all Zuker-optimal solutions for the RNA string contained in `examples/rna_string.txt`, then we run 
```bash
python tornado.py Z examples/rna_string.txt -dbn examples/dbn_string.txt
```
If we want to display a folding of maximum distance in .dbn format and also store a plot for the distance vector in `examples/plot_Zuker.png`, then we run
```bash
python tornado.py Z examples/rna_string.txt -dbn examples/dbn_string.txt -d -p -f examples/plot_Zuker
```

## File Structure

The core algorithms files are located in `algs` folder: 
- `zuker_backtrack.py` contains code for Zuker DP algorithm.
- `zuker_distance.py` contains code for computing the maximum distance and distance vector between a given folding and Zuker-optimal foldings.
- `nussinov.py` contains code for Nussinov DP algorithm and also code for computing the maximum distance and distance vector between a given folding and Nussinov-optimal foldings.

## Visualize RNA foldings

VARNA is an applet that allows you to visualize RNA secondary structure from dot-bracket-notation strings: https://varna.lri.fr/index.php?lang=en&page=downloads&css=varna
