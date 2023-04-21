# Read T-Walk HDF5 results 

* Extract parameters
* Split into walkers/chains
* Convert to arviz inference object
* Compute diagnostics R-hat parameter and ESS
* Plot with arviz and corner

# Basic usage

    pip install -r requirements.txt
    python3 reader.py /path/to/twalk-result.hdf5
