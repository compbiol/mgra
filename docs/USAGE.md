Usage instructions for MGRA
=============================
In this manual it is assumed that MGRA is properly installed.

    Usage: mgra [-h] [-c CONFIGURATION_FILE] [-g INPUT_GENOMES] [-o OUTPUT_DIR] [--debug] [--version]

Supported arguments:

	-h [ --help ]        		show help message and exit
	-v [ --version ]     		show program's version number and exit
	-c [ --config ] <arg>		input configure file
	-g [ --genomes ] <arg>		input file which contains genomes
	-o [ --output_dir ] <arg> 	path to the output directory
	-d [ --debug ]			enable debug output (default: False)

Examples
--------

You can try MGRA on the provided ready-to-use examples:

    mgra -c examples/Xchr/x_chr.cfg -f grimm -g examples/Xchr/xchr.txt -o examples/Xchr/out
    mgra -c examples/mam6/sim.cfg -f grimm -g examples/Xchr/blocks.txt -o examples/mam6/out

Input 
----- 
The input data are described by a configuration file and file, which contain genomes. 

###Configuration file

See *docs/CONFIG_FILE.md* file. 

###Input genomes 

MGRA currently supports two formats of synteny blocks: 

* InferCARs
* GRIMM

These formats are self-expository - please see the sample input files for details. Generally, it is required that each synteny block is present in the given genomes in a single copy.

Output
------
After running MGRA, an output directory will contain:

* folder __genomes__: detailed description of ancestor adjacencies
* folder __transformations__: detailed description of ancestor transformation
* __history_stat.txt__: Statistical information about DCJ-operations. 
* __sorted_full_history.txt__: List of sorted DCJ-operations.
* folder __debug__ :  Debug information about the progress after each stage. (if --debug was specified)

