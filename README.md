# h5logiqa
Accessing Hi-C data stored in LOGIQA HDF5 files

## Requirements

The API requires to following dependencies:
- Python (&ge; 3) with the *h5py*, *numpy*, and *scipy* modules
- Gnuplot (&ge; 4.6.6, not tested with more recent versions)

## Usage

Basic usage below:

### Loading the module

    >>> from h5logiqa import H5Dset
    
### Opening an HDF5 file

    >>> mat = H5Dset("/home/data/hESC.h5")
    
If you want to include genes to the figure, you need to pass the genome assembly and an SQLite database as well:

    >>> mat = H5Dset("/home/data/hESC.h5", "hg19", "/home/data/genes.db")
    
### Retrieving data for a given genomic region

    >>> mat.get("chr3", 181399711, 181459711)
    
Please note that the `get()` function does not return anything.
By default, frequency counts are retrieved. To retrieve loops' dispersion (90%, 70% or 50%) pass the desired dispersion with the `sampling` argument.

    >>> mat.get("chr3", 181399711, 181459711, sampling=70)
    
You can filter loops with a minimum contact count threshold with the `count` argument, or with a maximum dispersion threshold with the `disp` argument:

    >>> mat.get("chr3", 181399711, 181459711, sampling=50, count=5, disp=10)
    
### Plotting data

Generate a figure of the retrieved data:

    >>> mat.plot()
    /tmp/tmp.4gHaOMFGjz
    
The `plot()` method returns the path of the generated image. If Gnuplot is not in your `PATH`, pass the binary's path as follows:

    >>> mat.plot(gnuplot="/programs/gnuplot-4.6.6/src/gnuplot")
    
You can tweak the scale limit with the `vmax` argument, and generate a triangular display instead of a square one with the `astriu` argument:

    >>> mat.plot(vmax=100, astriu=True)
    