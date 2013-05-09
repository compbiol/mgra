#ifndef WDOTS_H_
#define WDOTS_H_

#include <string>
#include <vector>
#include <fstream>

namespace writer {
struct Wdots { 
	Wdots(std::string name_file);

	void save_dot(bool flag);
	void write_legend_dot(const std::vector<Genome>& genomes, const std::vector<std::string>& output); 
private: 
	std::ofstream output;
}; 
} 

#endif

