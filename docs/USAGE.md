Usage instructions for MGRA
=============================
In this manual it is assumed that MGRA is properly installed.

    Usage: mgra  {-g FILENAME|-i FILENAME} [-s] [-d] -o DIRNAME -c FILENAME [--] [--version] [-h] 

Supported arguments:

	-h [ --help ]        			show help message and exit
	-v [ --version ]     			show program's version number and exit
	-c [ --config ] <arg>			input configure file
	-g [ --grimm_genomes ] <arg>		input file which contains genomes in grimm format
	-i [ --infercars_genomes ] <arg>	input file which contains genomes in infercars format
	-o [ --output_dir ] <arg> 		path to the output directory
	-d [ --debug ]				enable debug output (default: False)
	-s [ --saves ]				enable output of saves in json format (default: False)

Examples
--------

You can try MGRA on the provided ready-to-use examples:

    mgra -c examples/Xchr/x_chr.cfg -g examples/Xchr/xchr.txt -o examples/Xchr/out
    mgra -c examples/mam6/sim.cfg -g examples/Xchr/blocks.txt -o examples/mam6/out

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
* folder __trees__: detailed description of recovered phylogenetic trees (reservation for future)
* folder __input__: config and blocks files
* folder __debug__ :  Debug information about the progress after each stage. (if --debug was specified)
* folder __saves__ :  Saves information about the progress after each stage. (if --saves was specified) (reservation for future)
* folder __mgra.log__ :  File with log informaton.

