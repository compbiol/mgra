#ifndef READER_H_
#define READER_H_

#include "defined.h"
#include "pconf.h"

namespace reader { 
  typedef structure::Genome genome_t;
  
  std::vector<genome_t> read_genomes(const ProblemInstance<Mcolor>& cfg); 
	
  void read_infercars(const ProblemInstance<Mcolor>& cfg, std::vector<genome_t>& genome); 

  void read_grimm(const ProblemInstance<Mcolor>& cfg, std::vector<genome_t>& genome);

  std::unordered_map<std::string, std::vector<std::string> > read_cfg_file(const std::string& name_cfg_file);
} 

#endif
