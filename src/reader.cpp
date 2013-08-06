#include "reader.h"

std::vector<Genome> reader::read_genomes(const ProblemInstance<Mcolor>& cfg) { 

  if (cfg.get_count_genomes() < 2) {
    std::cerr << "ERROR: at least two input genomes required" << std::endl;
    exit(1);
  }
	
  std::vector<Genome> genomes(cfg.get_count_genomes());

  /*Reading file in current format*/
  if (cfg.get_blk_format() == "infercars") {
    reader::read_infercars(cfg, genomes);
  } else if (cfg.get_blk_format() == "grimm") {
    reader::read_grimm(cfg, genomes);
  } else {
    std::cerr << "ERROR: unknown synteny blocks format" << cfg.get_blk_format() << std::endl;
    exit(1);
  }

  return genomes;
} 

void reader::read_infercars(const ProblemInstance<Mcolor>& cfg, std::vector<Genome>& genome) {
	
  std::ifstream input(cfg.get_blk_file().c_str());
  if(!input) {
    std::cerr << "Unable to open " << cfg.get_blk_file() << std::endl;
    exit(1);
  }
    
  std::string gene;
  std::vector<std::string> chromosome(cfg.get_count_genomes());
  std::vector<std::string> sign(cfg.get_count_genomes()); 
  std::vector<size_t> start_block(cfg.get_count_genomes());  
  std::vector<size_t> end_block(cfg.get_count_genomes()); 
  std::vector<size_t> count_block(cfg.get_count_genomes());

  while(!input.eof()) {
    std::string line;
    std::getline(input, line);
    line = trim(line);
    if (line.empty()) {
      if (gene.empty()) { 
	continue;
      } 

      bool ambig = false;
      for(int i = 0; (i < count_block.size()) && !ambig; ++i) { 
	if (count_block[i] != 1) { 
	  ambig = true;
	}  
      } 	
      
      if (ambig) {
	std::cerr << "Ambiguous block: " << gene << std::endl;
      } else {
	for(int i = 0; i < count_block.size(); ++i) { 
	  if (count_block[i] == 1) {
	    int sign_ = (sign[i] == "+") ? +1: -1;
	    genome[i].insert(gene, chromosome[i], (start_block[i] + end_block[i]) / 2, sign_, std::min(start_block[i], end_block[i]), std::max(start_block[i], end_block[i]));         	   
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
      if (cfg.is_genome_name(genome_name)) {
	size_t k = cfg.get_genome_number(genome_name);
	istr >> chromosome[k] >> start_block[k] >> end_block[k] >> sign[k];
	++count_block[k];
      } else { 
	std::cerr << "Unknown genome name: " << genome_name << std::endl;
      } 
    }
  } 
  input.close();
}


void reader::read_grimm(const ProblemInstance<Mcolor>& cfg, std::vector<Genome>& genome) {
  std::ifstream input(cfg.get_blk_file().c_str());;
    
  if (!input) {
    std::cerr << "Unable to open " << cfg.get_blk_file() << std::endl;
    exit(1);
  }

  size_t nchr = 0;
  size_t number_genome = 0;
  while(!input.eof()) {
    std::string line;
    std::getline(input, line);
    line = trim(line);
	
    if (line[0] == '>') {		
      line = trim(line.substr(1));
      if (cfg.is_genome_name(line)) {
	number_genome = cfg.get_genome_number(line);
	nchr = 0; 
      } else { 
	std::clog << "Unknown genome: " << line << std::endl;
      } 
    } else if (!line.empty() && (line[0] != '#')) {
      ++nchr;			
      std::string chr = "chr" + toString<size_t>(nchr);
      std::istringstream istr(line);
      size_t genord = 0;
      while (!istr.eof()) {
	std::string gene;
	istr >> gene;	
	if (gene.empty()) {  
	  break;
	} else if (gene[0] == '$') { //FIXME
	  break;
	} else if (gene[0] == '@') { //FIXME
	  genome[number_genome].registrate_circular_chr(chr);
	  break;
	} else {
	  int sign = (gene[0] == '-')? -1: +1; 
	  ++genord;
	  if (gene[0] == '-' || gene[0] == '+') {
	    gene = gene.substr(1);
	  }
	  genome[number_genome].insert(gene, chr, genord, sign, genord, genord); 
	} 
      }
    } 
  }
  input.close();
}

std::unordered_map<std::string, std::vector<std::string> > reader::read_cfg_file(const std::string& name_cfg_file) { 
  std::unordered_map<std::string, std::vector<std::string> > problem_config;
	
  std::ifstream input(name_cfg_file.c_str());
	
  if (!input) {
    std::cerr << "ERROR: Cannot open " << name_cfg_file << std::endl;
    exit(1); 
  }
	
  std::string section;
  while(!input.eof()) {
    std::string line;
    std::getline(input, line);
    line = trim(line);

    if (line[0] == '[' && line[line.size() - 1] == ']') {
      section = line;
    } else if (!line.empty() && (line[0] != '#')) {
      problem_config[section].push_back(line);
    }
  } 
  input.close();
	
  return problem_config;
} 

