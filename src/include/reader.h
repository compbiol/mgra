#ifndef READER_H_
#define READER_H_

#include "defined.h"
#include "structures/pconf.h"

namespace reader { 
  typedef structure::Genome genome_t;
  typedef structure::Mcolor mcolor_t;
  	
  inline std::string trim(std::string s, std::string const & drop = " \t\r\n") {
    s = s.erase(s.find_last_not_of(drop) + 1);
    return s.erase(0, s.find_first_not_of(drop));
  }

  std::vector<structure::Genome> read_infercars(ProblemInstance<mcolor_t> const & cfg, fs::path const & path_to_file); 
  std::vector<structure::Genome> read_grimm(ProblemInstance<mcolor_t> const & cfg, fs::path const & path_to_file);
  std::unordered_map<std::string, std::vector<std::string> > read_cfg_file(fs::path const & path_to_file);
} 

#endif
