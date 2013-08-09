#ifndef WGENOME_H_
#define WGENOME_H_

namespace writer {

template <class genome_t>
struct Wgenome { 
  void save_genom_in_text_format(const std::string& outname, const genome_t& AllChr, bool isEmptyTarget);
};
 
} 

//rename and move to namespace writer
template <class genome_t>
void writer::Wgenome<genome_t>::save_genom_in_text_format(const std::string& outname, const genome_t& AllChr, bool isEmptyTarget) { 
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

		const path_t& path = chromosome.first;

		if (chromosome.second) {
			++ncirc;
			lcirc += path.size();
			out << "# circular ";
		} else {
			out << "# linear ";
		}

		out << ChrTitle << " of length " << path.size() << " follows:" << std::endl;

		if ((*path.begin())[0] == '+' || (*path.rbegin())[0] == '+') {
			for(const auto &block : path) {
				out << block<< " ";
			}
		} else {
			for(auto ip = path.rbegin(); ip != path.rend(); ++ip) {
				vertex_t e = *ip;
				if (e[0] == '-') { 
					e[0] = '+'; 
				} else if (e[0] == '+') { 
					e[0] = '-';
				} 
				out << e << " ";
			}
		}

		out << "$" << std::endl;
	}

	out << std::endl << "# Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << std::endl;
	std::cout << std::endl << "Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << std::endl; //FIXME

	if (ncirc) {
		out << "#\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl;
		std::cout << "\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl; //FIXME
	}

	out.close();
}


#endif
