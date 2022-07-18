# rna-diameter

Example uses are shown in `example.py`

Core files are: 
- `zuker_backtrack.py` contains code for Zuker DP algorithm that records breadcrumbs
- `zuker_distance.py` contains code for computing max distance and distance vector for solutions yielded by Zuker's algorithm
- `nussinov.py` contains code for Nussinov DP algorithm and also code for computing max distance 
and distance vector for solutions yielded by Nussinov's algorithm.

Testing files are:
- `test_zuker_backtrack.py`
- `test_nussinov.py`

To run examples:
- Create a new Conda environment: `conda create -n newenv python=3.9`
- Activate Conda environment: `conda activate newenv`
- Run file: `python example.py`
  
