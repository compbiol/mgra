#ifndef GENOME_MATCH_ 
#define GENOME_MATCH_ 

#include "defined.h"
#include "pconf.h"

struct genome_match { 
  typedef structure::Genome genome_t;
  typedef structure::Mcolor mcolor_t;
  
  typedef std::unordered_map<std::string, size_t> gen2num; 
 
  static void init_name_genomes(const ProblemInstance<mcolor_t>& cfg, const std::vector<genome_t>& genomes) {
    number_to_genome.resize(genomes.size());
	
    for(size_t i = 0; i < genomes.size(); ++i) { 
      number_to_genome[i] = cfg.get_priority_name(i);
      genome_to_number.insert(std::make_pair(number_to_genome[i], i));
    }
  } 
  
  inline static bool member_name (const std::string& i) { 
    return (genome_to_number.find(i) != genome_to_number.end());
  } 
	
  inline static size_t get_number(const std::string& s) {
    assert(member_name(s));
    return genome_to_number.find(s)->second;
  }

  static mcolor_t name_to_mcolor(const std::string& name) {
    mcolor_t current;

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

  static std::string mcolor_to_name(mcolor_t const & S) {
    if (S.empty()) { 
      return "\\ensuremath{\\emptyset}";
    }

    std::ostringstream os;
    for(auto is = S.cbegin(); is != S.cend(); ++is) {
      std::string sym = number_to_genome[is->first]; 
      for(size_t i = 0; i < is->second; ++i) { 
	os << sym;
      }
    }
    return os.str();
  }
     
private: 
  static std::vector<std::string> number_to_genome;
  static gen2num genome_to_number;    		
};

  

#endif
