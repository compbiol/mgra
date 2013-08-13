#ifndef WGENOME_H_
#define WGENOME_H_

namespace writer {

template <class genome_t>
struct Wgenome { 
  void save_genom_in_text_format(const std::string& outname, const Genome& AllChr, bool isEmptyTarget);
};
 
} 

template <class genome_t>
void writer::Wgenome<genome_t>::save_genom_in_text_format(const std::string& outname, const Genome& AllChr, bool isEmptyTarget) { 
	std::ofstream out((outname + ".gen").c_str());

	out << "# Genome " << outname << std::endl;

	std::string ChrTitle;

	if (isEmptyTarget) { 
		ChrTitle = "chromosome"; 
	} else { 
		ChrTitle = "CAR";
	} 

	size_t ncirc = 0; 
	size_t lcirc = 0;

	for(const auto &chromosome : AllChr) {
		out << std::endl;

		if (chromosome.second.is_circular()) {
			++ncirc;
			lcirc += chromosome.second.size();
			out << "# circular ";
		} else {
			out << "# linear ";
		}

		out << ChrTitle << " of length " << chromosome.second.size() << " follows:" << std::endl;

		if (chromosome.second.begin()->second.second == 1 || (--chromosome.second.end())->second.second == 1) {
			for(const auto &block : chromosome.second) {
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

	out << std::endl << "# Reconstructed genome " << outname << " has " << AllChr.count_chromosome() << " " << ChrTitle << "s" << std::endl;
	//std::cout << std::endl << "Reconstructed genome " << outname << " has " << AllChr.count_chromosome() << " " << ChrTitle << "s" << std::endl; //FIXME

	if (ncirc) {
		out << "#\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl;
		//std::cout << "\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl; //FIXME
	}

	out.close();
}

#endif
