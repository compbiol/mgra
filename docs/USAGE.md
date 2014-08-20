Usage instructions for MGRA
=============================
In this manual it is assumed that MGRA is properly installed.

    Usage: mgra [-h] [-c CONFIGURE_FILE] [-f {grimm,infercars}] [-g FILE_WITH_GENOMES]
                     [-o OUTPUT_DIR] [--colorscheme COLORSCHEME_ARG] [--debug] [--version]

Supported arguments:
	-h [ --help ]           show help message and exit
	-v [ --version ]        show program's version number and exit
	-c [ --config ] arg     input configure file
	-f [ --format ] {grimm, infercars} input format file for genomes file
	-g [ --genomes ] arg    input file which contains genomes
	-o [ --output_dir ] arg path to the output directory
  	--colorscheme arg       colorscheme, which used in output breakpoint graph (default : use mgra colorscheme)
	-d [ --debug ]          enable debug output (default: False)
	--assembly		Reservation for future. Does not support now.

Examples
---------

You can try MGRA on the provided ready-to-use examples:

    mgra -c examples/Xchr/x_chr.cfg -f grimm -g examples/Xchr/xchr.txt -o examples/Xchr/out
    mgra -c examples/mam6/sim.cfg -f grimm -g examples/Xchr/blocks.txt -o examples/mam6/out

Input 
----- 
MGRA takes as input: 

	* Names of genomes
	* Algorithmic parameters
	* Phylogenetic tree containing all genomes in NEWICK format

All these parameters should be described in a single configuration file.
See the example of such file below.

Output
------
After running MGRA, an output directory will contain:

* folder __genomes__: detailed description of ancestor adjacencies
* folder __transformations__: detailed description of ancestor transformation
* __history_stat.txt__: Statistical information about DCJ-operations. 
* __last_graph.dot__: The breakpoint graph (without complete multi-edges) after process in the graphviz .dot format (http://graphviz.org)
* __legend.dot__: Input colored phylogenetic tree
* __sorted_full_history.txt__: List of sorted DCJ-operations.
* folder __debug__ :  Debug information about the progress after each stage. (if --debug was specified)

Configure file
--------------

See *docs/CONFIG_FILE.md* file. 

Synteny blocks file 
-------------------
MGRA currently supports two formats of synteny blocks: 

* infercars
* grimm

These formats are self-expository -- please see the sample input files for details. Generally, it is required that each synteny block is present in the given genomes in a single copy.

Output transformation file
-------------------

Output genome file
-------------------


