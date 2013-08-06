#ifndef READER_H_
#define READER_H_

#include <cstdlib>

#include "pconf.h"
#include "genome.h"
#include "mcolor.h"

namespace reader { 
  std::vector<Genome> read_genomes(const ProblemInstance<Mcolor>& cfg); 
	
  void read_infercars(const ProblemInstance<Mcolor>& cfg, std::vector<Genome>& genome); 

  void read_grimm(const ProblemInstance<Mcolor>& cfg, std::vector<Genome>& genome);

  std::unordered_map<std::string, std::vector<std::string> > read_cfg_file(const std::string& name_cfg_file);
} 

#endif
