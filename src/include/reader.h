#ifndef READER_H_
#define READER_H_

#include "defined.h"
#include "pconf.h"

namespace reader { 
  typedef structure::Genome genome_t;
  typedef structure::Mcolor mcolor_t;
  	
  std::vector<structure::Genome> read_infercars(ProblemInstance<mcolor_t> const & cfg, fs::path const & path_to_file); 
  std::vector<structure::Genome> read_grimm(ProblemInstance<mcolor_t> const & cfg, fs::path const & path_to_file);
  std::unordered_map<std::string, std::vector<std::string> > read_cfg_file(fs::path const & path_to_file);
} 

#endif
