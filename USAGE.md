USAGE
=======

In this manual it is assumed that "MGRA" is builded and copied in current work directory.

Basic usage
------
The easiest way to run "MGRA" is to type:

	Usage: ./mgra configure_file

Directory "examples" contains two sets of mammalian genomes. For example, consider dataset Xchr

	./mgra x_chr.cfg

Configure file and synteny blocks file described below in details.

Configure file
------
MGRA requires a single parameter -- a plain text file describing the problem instance and related options ("configuration").

The configuration file should describe follow sections:

* [Genomes]
Each genome may have a multiple aliases by which it is recognized in the synteny blocks data; but the main/first alias is priority name which MGRA will use to denote multicolors and in output files.

* [Blocks]
    * format <name_format>
    Information about format of synteny blocks in synteny blocks file. About which formats MGRA supports see below.

    * file <name_file>
    Filename which contains input synteny blocks. Important note -- "MGRA" requires that file with this name located in the same directory as the configuration file.

* [Trees]
A phylogenetic tree may be given as a single tree or as a number of subtrees in 
the Newick format: http://en.wikipedia.org/wiki/Newick_format

* [Algorithm]
    * stages <number>
    This parameter determines the number of steps that the algorithm will make. About detailed description of each step refer to the article. Minimal value = 0, maximum value = 4. (e.g. stage 3, MGRA runs on first, second and third stage).
    * rounds <number>
    This parameter determines the number of rounds that the algorithm will make. About detailed description of each round refer to the article. Minimal value = 0, maximum value = 3. (e.g. round 2, MGRA runs on first and second round).
    * bruteforce <number>
    This parameter is switched on experimental stage which can resolve complex component after latest stage. The number describes maximum count vertices in component. We recommend set this number to 12.
If you set a bigger number, it will lead to long-term data processing. If you set too small number, some components may remain unprocessed. 

* [Graphs]
    * filename <name>
    This parameter sets the name of file which is used to output breakpoint graphs after each stage in .dot files (e.g if name = stage then "stage0.dot", "stage1.dot" and etc).

    * colorscheme <parameter>
    This parameret sets the colorscheme, which will be used when writing breakpoint graph after stage in dot file. Valid values ​​are the parameters that are specified in the graphviz .dot format. If not specified MGRA uses a standard scheme.
    
Synteny blocks file 
------
MGRA currently supports two formats of synteny blocks: 

* "infercars" 
* "grimm"

These formats are self-expository -- please see the sample input files for details. Generally, it is required that each synteny block is present in the given genomes in a single copy.

Output files
------
Most important information MGRA outputs to the standard output. After running MGRA, an output directory will contain:

* Ancestor genomes represented as permutations of the synteny blocks - ".gen"

The ancestor genome(s) are saved in GRIMM format.

* Transformation to each ancestor genome - "*.trs"

Each .gen file is accompanied with .trs file describing transformation of this genome along the incident branch towards the root genome.

* Input colored phylogenetic tree - "legend.dot"

* Breakpoint graphs after each stage - "stage*.dot"

The breakpoint graph (without complete multi-edges) after each stage in the graphviz .dot format (http://graphviz.org) into "stage0.dot" (initial graph), "stage1.dot" (after MGRA Stage 1), "stage2.dot" etc. To obtain image in PDF (PNG, JPG, etc. are also possible) format from .dot, one can run:

	neato -Tpdf -O stage1.dot

that will produce file "stage1.dot.pdf" in PDF format. 

* Statistical report - "stats.txt"

Statistical information about the progress after each stage.
