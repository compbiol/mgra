#ifndef GENOME_MATCH_ 
#define GENOME_MATCH_ 

#include "defined.h"
#include "structures/pconf.h"

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
  
  inline static bool member_name (std::string const & i) { 
    return (genome_to_number.find(i) != genome_to_number.end());
  } 
	
  inline static size_t get_number(std::string const & s) {
    assert(member_name(s));
    return genome_to_number.find(s)->second;
  }

  static mcolor_t name_to_mcolor(std::string const & temp) {
    mcolor_t answer; 
    if (temp[0] == '{' && temp[temp.length() - 1] == '}') { 
      std::string current = "";
      for (size_t i = 1; i < temp.length(); ++i) { 
        if (temp[i] == ',' || temp[i] == '}') { 
          if (genome_to_number.find(current) != genome_to_number.end()) {  
            answer.insert(genome_to_number.find(current)->second);
          }   
          current = "";
        } else if (std::isalpha(temp[i]) || std::isdigit(temp[i])) { 
          current += temp[i];
        } else { 
          std::cerr << "Bad format target " << temp << std::endl;
          exit(1);
        }  
      } 
    } else {
      std::cerr << "Bad format target " << temp << std::endl;
      exit(1);
    }
    return answer;
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
