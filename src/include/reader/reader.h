#ifndef READER_H_
#define READER_H_

#include "defined.hpp"
#include "io/path_helper.hpp"
#include "config/config_struct.hpp"

namespace reader { 
  using genome_t = structure::Genome;
  	
  std::vector<genome_t> read_infercars(std::string const & path_to_file); 
  std::vector<genome_t> read_grimm(std::string const & path_to_file);
} 

#endif
