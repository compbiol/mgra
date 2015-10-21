#include "reader/reader.h"

std::vector<structure::Genome> reader::read_infercars(std::string const & path_to_file) {
  std::ifstream input(path_to_file);

  if (!input) {
    std::cerr << "Unable to open " << path_to_file << std::endl;
    exit(1);
  }

  std::vector<genome_t> genomes(cfg::get().get_count_genomes());    
  std::string gene;
  std::vector<std::string> chromosome(cfg::get().get_count_genomes());
  std::vector<std::string> sign(cfg::get().get_count_genomes()); 
  std::vector<size_t> start_block(cfg::get().get_count_genomes());  
  std::vector<size_t> end_block(cfg::get().get_count_genomes()); 
  std::vector<size_t> count_block(cfg::get().get_count_genomes());

  while(!input.eof()) {
    std::string line;
    std::getline(input, line);
    boost::trim(line);
    if (line.empty()) {
      if (gene.empty()) { 
        continue;
      } 

      bool ambig = false;
      for(size_t i = 0; (i < count_block.size()) && !ambig; ++i) { 
        ambig = (count_block[i] != 1);
      } 	
      
      if (ambig) {
        std::cerr << "Ambiguous block: " << gene << std::endl;
      } else {
        for(size_t i = 0; i < count_block.size(); ++i) { 
          if (count_block[i] == 1) {
            int sign_ = (sign[i] == "+") ? +1: -1;
            genomes[i].insert(gene, chromosome[i], (start_block[i] + end_block[i]) / 2, sign_); 
            //, std::min(start_block[i], end_block[i]), std::max(start_block[i], end_block[i]));         	   
          }
        } 
        gene.clear();
      } 
    } else if(line[0] == '>') {
      gene = line.substr(1);
      std::fill(count_block.begin(), count_block.end(), 0);
    } else if (line[0] != '#') { 
      line[line.find(".")] = ' ';
      line[line.find(":")] = ' ';
      line[line.find("-")] = ' ';
      std::istringstream istr(line);
      std::string genome_name;
      istr >> genome_name;
      if (cfg::get().is_genome_name(genome_name)) {
        size_t k = cfg::get().get_genome_number(genome_name);
        istr >> chromosome[k] >> start_block[k] >> end_block[k] >> sign[k];
        ++count_block[k];
      } else { 
        std::cerr << "Unknown genome name: " << genome_name << std::endl;
      } 
    }
  } 
  input.close();

  return genomes;
}


std::vector<structure::Genome> reader::read_grimm(std::string const & path_to_file) {
  std::ifstream input(path_to_file);
    
  if (!input) {
    std::cerr << "Unable to open " << path_to_file << std::endl;
    exit(1);
  }

  size_t nchr = 0;
  size_t number_genome = 0;
  std::vector<genome_t> genomes(cfg::get().get_count_genomes());   
  auto inserter_lambda = [&] (std::string gene, std::string const & chr, size_t offset) {
    int sign = (gene[0] == '-')? -1: +1; 
    if (gene[0] == '-' || gene[0] == '+') {
      gene = gene.substr(1);
    }
    genomes[number_genome].insert(gene, chr, offset, sign); //, genord, genord); 
  };

  while(!input.eof()) {
    std::string line;
    std::getline(input, line);
    boost::trim(line);
	
    if (line[0] == '>') {		
      line = boost::trim_left_copy(line.substr(1));
      if (cfg::get().is_genome_name(line)) {
        number_genome = cfg::get().get_genome_number(line);
        nchr = 0; 
      } else { 
        std::cerr << "Unknown genome: " << line << std::endl;
      } 
    } else if (!line.empty() && (line[0] != '#')) {
      ++nchr;			
      std::string chr = "chr" + std::to_string(nchr);
      std::istringstream istr(line);
      size_t genord = 0;
      while (!istr.eof()) {
      	std::string gene;
      	istr >> gene;	
      	if (gene.empty()) {  
      	  break;
      	} else if (gene[0] == '$') {
      	  break;
        } else if (gene.find("$") != std::string::npos) { 
          inserter_lambda(gene.substr(0, gene.find("$")), chr, ++genord);
          break;
      	} else if (gene[0] == '@') {
      	  genomes[number_genome].registrate_circular_chr(chr);
      	  break;
        } else if (gene.find("@") != std::string::npos) { 
          inserter_lambda(gene.substr(0, gene.find("@")), chr, ++genord);
          genomes[number_genome].registrate_circular_chr(chr);
          break;
      	} else {
          inserter_lambda(gene, chr, ++genord);
      	} 
      }
    } 
  }

  input.close();
  return genomes;
}
