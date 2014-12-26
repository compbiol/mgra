#ifndef WGENOME_H_
#define WGENOME_H_

namespace writer {

template <class genome_t>
struct Wgenome {

  explicit Wgenome(std::string const & path) 
  : m_path(path)
  { 
  }

  void save_genomes(std::vector<genome_t> const & genomes, bool isEmptyTarget) const { 
    for (auto const & genome: genomes) {
      save_genome_in_text_format(genome, isEmptyTarget);
    } 
  } 

private: 
  void save_genome_in_text_format(genome_t const & genome, bool isEmptyTarget) const;

private: 
  std::string m_path;

private:
  DECL_LOGGER("GenomeWriter");
};

} 

template <class genome_t>
void writer::Wgenome<genome_t>::save_genome_in_text_format(genome_t const & genome, bool isEmptyTarget) const { 
  std::string const & outname = genome.get_name();  
  std::ofstream out(path::append_path(m_path, (outname + ".gen")));

  out << "# Genome " << genome.get_name() << std::endl;

  std::string chr_title;

  if (isEmptyTarget) { 
    chr_title = "chromosome"; 
  } else { 
    chr_title = "CAR";
  } 

  size_t number_circular = 0; 
  size_t length_circular = 0;

  for(auto const & chromosome : genome) {
    out << std::endl;

    if (chromosome.second.is_circular()) {
      ++number_circular;
      length_circular += chromosome.second.size();
      out << "# circular ";
    } else {
      out << "# linear ";
    }

    out << chr_title << " of length " << chromosome.second.size() << " follows:" << std::endl;

    if (chromosome.second.begin()->second.second == 1 || (--chromosome.second.end())->second.second == 1) {
      for(auto const & block : chromosome.second) {
        out << ((block.second.second == -1)?"-":"+") << block.second.first << " ";
      }
    } else {
      for(auto ip = chromosome.second.crbegin(); ip != chromosome.second.crend(); ++ip) {
      	out << ((ip->second.second == -1)?"+":"-") << ip->second.first << " ";
      }
    }

    if (chromosome.second.is_circular()) {			
      out << "@" << std::endl;
    } else { 
      out << "$" << std::endl;
    } 
  }

  {
    std::ostringstream ost; 
    ost << "Reconstructed genome " << outname << " has " << genome.count_chromosome() << " " 
        << chr_title << "s" << " out of which " << number_circular << " are circular of total length " << length_circular;
    INFO(ost.str())        
  }  
  
  out << std::endl << "# Reconstructed genome " << outname << " has " << genome.count_chromosome() << " " << chr_title << "s" << std::endl;
  if (number_circular) {
    out << "#\tout of which " << number_circular << " are circular of total length " << length_circular << std::endl;
  }
  out.close();
}

#endif
