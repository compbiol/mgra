Configure file description
==========================

MGRA requires a plain text file describing the problem instance and related options ("configuration").
The configuration file should describe follow sections:

###[Genomes]

Each genome may have a multiple aliases by which it is recognized in the synteny blocks data; but the main/first alias is priority name which MGRA will use to denote multicolors and in output files.

###[Trees]

A phylogenetic tree may be given as a single tree or as a number of subtrees in 
the Newick format: http://en.wikipedia.org/wiki/Newick_format

###[Target]

Reservation for future. Does not support now.

###[Algorithm]

* stages [__number__]

This parameter determines the number of steps that the algorithm will make. About detailed description of each step refer to the article. Minimal value = 0, maximum value = 4. (e.g. stage 3, MGRA runs on first, second and third stage).

* rounds [__number__]

This parameter determines the number of rounds that the algorithm will make. About detailed description of each round refer to the article. Minimal value = 0, maximum value = 3. (e.g. round 2, MGRA runs on first and second round).

* bruteforce 

This parameter is switched on experimental stage which can resolve complex component after latest stage with BlossumV algorithm. It will lead to long-term data processing. 

* recostructed_tree
	
Reservation for future. Does not support now.

###[Completion]
