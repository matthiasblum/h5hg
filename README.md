# h5hg
Accessing NGS data (ChIP-Seq, Hi-C) stored in HDF5 files.

## Requirements

h5hg requires to following dependencies:
- Python (&ge; 3) with the *h5py*, *numpy*, and *scipy* modules.
- Gnuplot (&ge; 4.6.6, not tested with more recent versions) with the *pngcairo* terminal (only for plotting Hi-C data).

## Hi-C loops

Basic usage below:

### Loading the module

    >>> from h5hg import Matrix
    
### Opening an HDF5 file

    >>> mat = Matrix("/home/data/hESC-40kb.h5")
    
If you want to include genes to the figure, you need to pass the genome assembly and an SQLite database as well:

    >>> mat = Matrix("/home/data/hESC-40kb.h5", "hg19", "/home/data/genes.db")
    
The SQLite database can be generated with [Gégène](https://github.com/matthiasblum/gegene).
    
### Retrieving data for a given genomic region

    >>> mat.get("chr3", 181399711, 181459711)
    
The `get()` function does not return the data. They are stored internally in the `Matrix` object.
By default, frequency counts are retrieved. To retrieve loops' dispersion (90%, 70% or 50%) pass the desired dispersion with the `sampling` argument.

    >>> mat.get("chr3", 181399711, 181459711, sampling=70)
    
You can filter loops with a minimum contact count threshold with the `count` argument, or with a maximum dispersion threshold with the `disp` argument:

    >>> mat.get("chr3", 181399711, 181459711, sampling=50, count=5, disp=10)
    
Passing 0 for the `start` and `end` arguments will retrieve the entire chromosome.
    
### Plotting data

Generate a figure of the retrieved data:

    >>> mat.plot()
    /tmp/tmp.4gHaOMFGjz
    
The `plot()` method returns the path of the generated image. If Gnuplot is not in your `PATH`, pass the binary's path as follows:

    >>> mat.plot(gnuplot="/programs/gnuplot-4.6.6/src/gnuplot")
    
You can tweak the scale limit with the `vmax` argument, and generate a triangular display instead of a square one with the `astriu` argument:

    >>> mat.plot(vmax=100, astriu=True)
    
### Chaining

If you need to retrieve only one region at a time, you can chain the methods:

    >>> Matrix("/home/data/hESC-40kb.h5").get("chr1", 181399711, 181459711).plot()
    
## ChIP-Seq

### Fixed-size windows (bins)

    >>> from h5hg import BinArray
    >>> arr = BinArray("/home/data/hESC-H3K4me3.h5")
    >>> data = arr.get("chr8", 128723314, 128773314)

The `get()` function accepts the following optional arguments:
- *localqc*: if True, retrieves the localQCs in the genomic region (default: True).
- *wig*: if True, retrieves the wiggles in the genomic region (default: True).
- *wigu*: if True, Wiggles for unique reads are retrieved as well (default: True).
- *nreps*: minimal number of times a bin is passing *dRCI < 10%* in the 5-replicates mode (default: 0, disabled). Must be between 
 and 5.
 - *shrink*: reduce the resolution so `resolution × shrink < region_length` by merging genomic features (default: 0, disabled).

`data` is a dictionary with the following keys/values:
- *localqc*: list of localQC dispersions. If the intensity was lower than the background threshold or if the dispersion was greater than or equal to 10, the intensity is 0 and the dispersion -1.
- *localqcres*: int, resolution of bins.
- *wig*: list of wiggles.
- *wigu*: list of wiggles for unique reads.
- *wigres*: int, resolution of Wiggles.

If the genomic region is invalid, `data` is None.

### Peaks

    >>> from h5hg import PeakArray
    >>> data = PeakArray("/home/data/hESC-H3K4me3.bin").get("chr8", 128723314, 128773314)

`data` is a list of peaks. Each peak is stored as a tuple: (start, end, summit, -log10(pvalue)).

#### Format

Unlike Hi-C contacts, localQCs, and Wiggles, peaks are not stored in the HDF5 format, but in a custom binary format. Thus, retrieving a small subset is faster but it becomes particularly slow for large chuncks (> 10Mb).

