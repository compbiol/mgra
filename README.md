Compiling and Running MGRA
==========================

I. Introduction
---------------
For algorithmic description of MGRA refer to the following paper:

Max Alekseyev, Pavel Pevzner "Breakpoint Graphs and Ancestral Genome Reconstructions."
Genome Research, special issue "Genomics and Darwinism", 2009.

Note that MGRA got further development since the paper had been written.
In particular, the current release of MGRA has larger number of stages:
Stage 2 described in the paper corresponds to Stages 2 and 3 in the current release 
of MGRA; while a new Stage 4 is able to automatically resolve some complex breakpoints 
that previously required manual processing. 


II. Installation and compilation
--------------------------------
Compilation of MGRA sources is tested with GNU C++ (g++) ver. 4.3.
To compile MGRA, run:
	g++ -O -o mgra.bin mgra.cpp


II. Input data
--------------
MGRA binary requires a single parameter -- a plain text file describing the 
problem instance and related options ("configuration").
Sample problem configurations are provided in the `input' directory.

Most importantly, the configuration file should describe:

* given genomes: 
Each genome may have a multiple aliases by which it is recognized in the synteny
blocks data; but the main/first alias must be a single unique letter which MGRA 
will use to denote multicolors. For example, if M stands for `mouse' and R stands 
for `rat', then MR will refer to the multicolor of mouse and rat together.

* phylogenetic tree: 
A phylogenetic tree may be given as a single tree or as a number of subtrees in 
the Newick format: http://en.wikipedia.org/wiki/Newick_format
Each (sub)tree should be a binary tree with no branch distances or names of 
internal nodes specified (names are derived automatically from the leaf nodes).
However, it is possible to specify group of genomes referring to a subtree with 
unknown local topology as a whole. 
For example, for seven mammalian genomes: mouse(M), rat(R), dog(D), macaque(Q), 
human(H), and chimpanzee(C), if the opossum branch is known to be only somewhere 
between the rodents, carnivore, and primates, one can specify two subtrees: 
(MRDO,(Q,(H,C))) and ((M,R),ODQHC)describing the topology of primates and 
rodents respectively. There is no need to describe single leaf branches, 
such as the dog branch (D,OMRQHC) in this case; the leaf branches are always 
taken into consideration by default.

* synteny blocks: MGRA currently supports two formats of synteny blocks: 
`infercars' and `grimm'. These formats are self-expository -- please see the 
sample input files for details. Generally, it is required that each synteny 
blocks is present in every of the given genomes in a single copy.


III. Output data
----------------
Most important information MGRA outputs to the standard output. Statistical 
information about the progress after each stage goes to file `stats.txt'.

MGRA also outputs the breakpoint graph (without complete multi-edges) after each
stage in the graphviz .dot format (http://graphviz.org) into `stage0.dot' 
(initial graph), `stage1.dot' (after MGRA Stage 1), `stage2.dot' etc.
To obtain image in PDF (PNG, JPG, etc. are also possible) format from .dot, 
one can run:
	neato -Tpdf -O stage1.dot
that will produce file `stage1.dot.pdf' in PDF format. 

The resulting genome(s) are saved in GRIMM format in .gen files (e.g., MRD.gen).
In the default mode, each .gen file is accompanied with .trs file describing 
transformation of this genome along the incident branch towards the root genome.


IV. Contact information
-----------------------
Please send your questions/comments/bugreports to:
Max Alekseyev <maxal@cse.sc.edu>
