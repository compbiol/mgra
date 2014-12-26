#ifndef READER_H_
#define READER_H_

#include "defined.h"
#include "io/path_helper.hpp"
#include "structures/config_struct.hpp"

namespace reader { 
  typedef structure::Genome genome_t;
  	
  std::vector<genome_t> read_infercars(std::string const & path_to_file); 
  std::vector<genome_t> read_grimm(std::string const & path_to_file);
} 

#endif
