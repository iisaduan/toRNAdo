# rna-diameter



Core Algorithms files are: 
- `zuker_backtrack.py` contains code for Zuker DP algorithm that records breadcrumbs
- `zuker_distance.py` contains code for computing max distance and distance vector for solutions yielded by Zuker's algorithm
- `nussinov.py` contains code for Nussinov DP algorithm and also code for computing max distance 
and distance vector for solutions yielded by Nussinov's algorithm.

Testing files are:
- `test_zuker_backtrack.py`
- `test_nussinov.py`

To run testing files:
- Use `python -m test.test_zuker_backtrack`

To run examples:
- Create a new Conda environment: `conda create -n newenv python=3.9`
- Activate Conda environment: `conda activate newenv`
- Run an example using the Command Line Tool: in project root directory, run `python -m tornado rna_string.txt --dbn_file dbn_string.txt Z --d --h --f newplot`

To visualize RNA foldings:
- VARNA is an applet that allows you to visualize RNA secondary structure from dot-bracket-notation strings: https://varna.lri.fr/index.php?lang=en&page=downloads&css=varna
