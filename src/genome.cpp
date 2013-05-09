#include "genome.h"

void Genome::insert(std::string gene, std::string chromosome, size_t offset, int sign_, size_t start, size_t end) { 
	coord_t p = std::make_pair(chromosome, offset);  
	main_genome.insert(std::make_pair(p, gene));
	sign[gene] = sign_;
	orf2cpan[gene] = std::make_pair(chromosome, std::make_pair(start, end));
}

bool Genome::isCircular(std::string name) const { 
	if (circular_chromosome.find(name) != circular_chromosome.end()) { 
		return true; 
	} else { 
		return false;
	} 
} 

