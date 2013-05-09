#include "genome_match.h"

std::vector<orf_t> genome_match::number_to_genome;
genome_match::gen2num genome_match::genome_to_number;   

void genome_match::init_name_genomes(const std::vector<Genome>& genomes) {
  number_to_genome.resize(genomes.size());
	
  for(size_t i = 0; i < genomes.size(); ++i) { 
    number_to_genome[i] = ProblemInstance::get_name(i);
    genome_to_number.insert(std::make_pair(number_to_genome[i], i));
    std::clog << "Genome " << ProblemInstance::get_name(i) << " blocks: " << genomes[i].size() << std::endl;
  }
} 

Mcolor genome_match::name_to_mcolor(const std::string& name) {
	Mcolor current;

	for (size_t j = 0; j < name.size(); ++j) {
		std::string t = name.substr(j, 1);
		if (genome_to_number.find(t) == genome_to_number.end()) {
			std::cerr << "ERROR: Malformed multicolor " << name << std::endl;
			exit(1);
		}
		current.insert(genome_to_number.find(t)->second);
	}

	return current;
}

std::string genome_match::mcolor_to_name(const Mcolor& S) {
	if (S.empty()) { 
		return "\\ensuremath{\\emptyset}";
	}

	std::ostringstream os;

	for(auto is = S.cbegin(); is != S.cend(); ++is) {
		const std::string& sym = number_to_genome[is->first]; 
		for(size_t i = 0; i < is->second; ++i) { 
			os << sym;
		} 
	}

	return os.str();
}     
