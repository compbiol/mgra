#include "writer/Wstats.h"

writer::Wstats::Wstats(std::string name_file): write_parametres(5) {
	ofstat.open(name_file); 
} 

void writer::Wstats::print_all_statistics(int stage, Statistics<mbgraph_with_history<Mcolor> >& info, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& graph) { 
	if (stage == 0) { 
		println("Initial graph:");
#ifdef VERSION2
		println("... Unique blocks: " + toString(graph.count_vertex() / 2));
		println("... Vertex: " + toString(graph.count_vertex()));
#endif
	}  else { 
		println("After Stage " + toString(stage) + " graph:");
	} 

#ifdef VERSION2
	print_duplication_statistics(info.count_all()); 
#endif
	print_complete_edges(graph);
	print_connected_components(graph);
	print_rear_characters(info.get_compl_stat()); 


//#ifndef VERSION2
	//print_estimated_dist(stage, cfg, graph);
//#endif
	print_fair_edges(graph, info);
	print_not_compl_characters(info.get_no_compl_stat()); 
} 

void writer::Wstats::print_duplication_statistics(const std::vector<size_t>& answer) {
  ofstat << "... Duplication vertex: " << answer[0] << std::endl;
  ofstat << "... Insertion/deletion vertex: " << answer[1] << std::endl;
  ofstat << "... Count self loop: " << answer[2] << std::endl;
  ofstat << "... Colors is not one-to-one match: " << answer[3] << std::endl; 
}

////////////////////////////////////////////////////////
void writer::Wstats::print_connected_components(const mbgraph_with_history<Mcolor>& MBG) {
	// count connected components in the graph
	equivalence<vertex_t> C;
	
	for(auto is = MBG.begin_vertices(); is != MBG.end_vertices(); ++is) {
		C.addrel(*is, *is);
	} 

	for(auto lc = MBG.begin_local_graphs(); lc != MBG.end_local_graphs(); ++lc) { 
		for(auto il = lc->cbegin(); il != lc->cend(); ++il) {
			if (il->first != Infty && il->second != Infty) { 
				C.addrel(il->first, il->second);	
			} 
		}
	}

	C.update();

	std::map<vertex_t, std::set<vertex_t> > cls = C.get_eclasses<std::set<vertex_t> >();

	std::map<size_t, size_t> stx;
	for(auto ic = cls.begin(); ic != cls.end(); ++ic) {
		++stx[ic->second.size()];
	}

	ofstat << "... connected components:";
	for(auto is = stx.begin(); is != stx.end(); ++is) {
		ofstat << " " << is->first << "^" << is->second;
	}
	ofstat << std::endl;
}

void writer::Wstats::histStat(const mbgraph_with_history<Mcolor>& graph) { //FIXME
	ofstat << std::endl << "Total number of 2-breaks: " << graph.get_count_2breaks() << std::endl;

 	std::map<Mcolor, size_t> n2br;

	for(auto il = graph.crbegin_2break_history(); il != graph.crend_2break_history(); ++il) {
		++n2br[il->get_mcolor()];
	}

	for(auto im = n2br.begin(); im != n2br.end(); ++im) {
		ofstat << genome_match::mcolor_to_name(im->first) << "\t" << im->second << std::endl;
	}
	
	ofstat << std::endl;
#ifdef VERSION2
	ofstat << std::endl << "Total number of insertion/deletions events: " << graph.get_count_insdel() << std::endl;

 	std::map<Mcolor, size_t> ninsdel;

	for(auto il = graph.begin_insdel_history(); il != graph.end_insdel_history(); ++il) {
		++ninsdel[il->get_mcolor()];
	}

	for(auto im = ninsdel.begin(); im != ninsdel.end(); ++im) {
		ofstat << genome_match::mcolor_to_name(im->first) << "\t" << im->second << std::endl;
	}
	
	ofstat << std::endl;

	ofstat << std::endl << "Total number of (reverse) tandem duplication: " << graph.get_count_tandem_duplication() << std::endl;

 	std::map<Mcolor, size_t> ntd;

	for(auto il = graph.begin_tandem_duplication_history(); il != graph.end_tandem_duplication_history(); ++il) {
		++ntd[il->get_mcolor()];
	}

	for(auto im = ntd.begin(); im != ntd.end(); ++im) {
		ofstat << genome_match::mcolor_to_name(im->first) << "\t" << im->second << std::endl;
	}
	
	ofstat << std::endl;
#endif
}
////////////////////////////////////////////////////////
void writer::Wstats::print_fair_edges(const mbgraph_with_history<Mcolor>& MBG, Statistics<mbgraph_with_history<Mcolor>>& info) {
	// output H-subgraphs count
	ofstat << std::endl << "% Fair multi-edges count: " << std::endl << std::endl;

	std::map<std::pair<Mcolor, Mcolor>, size_t> Hcount = info.get_Hsubgraph();
	std::list<Mcolor> HCrow; 	//the set of multicolors that appear in H-subraphs

#ifndef HG_TONLY // adding other colors
	std::set<Mcolor> proc;
	for(auto ih = Hcount.cbegin(); ih != Hcount.cend(); ++ih) {
		if (proc.find(ih->first.first) == proc.end()) {
			HCrow.push_back(ih->first.first);
			proc.insert(ih->first.first);
		}

		if (proc.find(ih->first.second) == proc.end()) {
			HCrow.push_back(ih->first.second);
			proc.insert(ih->first.second);
		}
	}
#endif
	ofstat << "\\begin{table}[h]" << std::endl;
	ofstat << "\\centering \\begin{tabular}{|c||"; //FIXME::why twice
	for(int i = 0; i < HCrow.size(); ++i) { 
		ofstat << "c|";
	} 
	ofstat << "}" << std::endl;
	ofstat << "\\hline" << std::endl;

	for(auto ic = HCrow.cbegin(); ic != HCrow.cend(); ++ic) {
		ofstat << " & ${";
		if (MBG.is_T_consistent_color(*ic)) {
			ofstat << "\\bf ";
		} 
		ofstat <<  genome_match::mcolor_to_name(*ic) << "+}$";
	}

	println("\\\\");
	println("\\hline \\hline");

	for(auto Q1 = HCrow.cbegin(); Q1 != HCrow.cend(); ++Q1) {
		ofstat << "${";
		if (MBG.is_T_consistent_color(*Q1)) {
			ofstat << "\\bf ";
		} 
		ofstat << genome_match::mcolor_to_name(*Q1) << "+}$"; 

		for(auto Q2 = HCrow.cbegin(); Q2 != HCrow.cend(); ++Q2) {
			ofstat << " & ";

			if (Hcount.find(std::make_pair(*Q1, *Q2)) != Hcount.end()) {
				ofstat << " "; //"${";

				if (MBG.are_adjacent_branches(*Q1, *Q2)) { 
					ofstat << "{\\cellcolor[gray]{.9}}";
				} 

#ifdef HG_TONLY
				if (member(BranchLen, genome_match::mcolor_to_name(*Q1)) && member(BranchLen,  genome_match::mcolor_to_name(*Q2))) {
					long t = (long) (1e5 * Hcount[std::make_pair(*Q1, *Q2)] / (double) BranchLen[genome_match::mcolor_to_name(*Q1)] / (double) BranchLen[genome_match::mcolor_to_name(*Q2)]); 
					ofstat << "$" << t << "\\cdot 10^{-5}$ ";
				}
#else
				ofstat << Hcount[std::make_pair(*Q1, *Q2)] << " "; 
#endif
			}

			if (*Q1 == *Q2) {
				ofstat << "$\\star$";
			}
		}
		ofstat << " \\\\" << std::endl;
		ofstat << "\\hline" << std::endl;
	}
	print_close_table(false);
}

void writer::Wstats::print_complete_edges(const mbgraph_with_history<Mcolor>& graph) { 
	size_t nc = 0;
	ofstat << "... complete multiedges:";
	for(auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
		Mularcs<Mcolor> M = graph.get_adjacent_multiedges(*it);
		if (M.size() == 1 && M.cbegin()->second == graph.get_complete_color() && (*it < M.cbegin()->first || M.cbegin()->first == Infty)) {
			ofstat << " " << *it << "~" << M.cbegin()->first;
			++nc;
		}
    	}
	ofstat << "\t(total: " << nc << ")" << std::endl;
} 

void writer::Wstats::print_rear_characters(const std::vector<std::string>& info) { 
	ofstat << std::endl;
	ofstat << "% Rearrangement characters:" << std::endl << std::endl;
	print_start_table(5); 
	ofstat << "Multicolors & multiedges & simple vertices & simple multiedges & simple paths+cycles & irreg. multiedges\\\\" << std::endl;
	ofstat << "\\hline" << std::endl;

	for(auto im = info.cbegin(); im != info.cend(); ++im) {
		ofstat << *im << "\\\\" << std::endl;
	}
	print_close_table();	
} 

void writer::Wstats::print_not_compl_characters(const std::vector<std::string>& info) { 
	if (!info.empty()) { 
		ofstat << std::endl;
		ofstat << "% Edges for duplication block:" << std::endl << std::endl;
		print_start_table(1);
		ofstat << "Multicolors & multiedges\\\\" << std::endl;
		ofstat << "\\hline" << std::endl;
		for(auto im = info.cbegin(); im != info.cend(); ++im) {
			ofstat << *im << "\\\\" << std::endl;
		}
		print_close_table();	
	} 
} 

void writer::Wstats::print_estimated_dist(size_t stage, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& graph) { 
	ofstat << "% Estimated distances:" << std::endl << std::endl;
	/*print_start_table(graph.count_local_graphs());
	ofstat << "Stage " << stage;
	for(size_t i = 0; i < cfg.get_count_genomes(); ++i) { 
		ofstat << " & " << cfg.get_priority_name(i); 
	} 
	ofstat << " \\\\" << std::endl << "\\hline" << std::endl;

	size_t i = 0; 
	for(size_t i = 0; i < cfg.get_count_genomes(); ++i) { 
		ofstat << cfg.get_priority_name(i) << " & "; 
		
		for(size_t j = 0; j < graph.count_local_graphs(); ++j) {
			if (j > i) {
				ofstat << genome_dist(graph.get_local_graph(i), graph.get_local_graph(j), graph.get_obverce_graph())[2]; //FIXME
			}

			if (j == graph.count_local_graphs() - 1) { 
				ofstat << "\\\\";
			} else { 
				ofstat << " & ";			
			} 
		}
		ofstat << std::endl;
	}*/
	print_close_table();
} 

void writer::Wstats::print_start_table(size_t count_column) { 
	ofstat << "\\begin{table}[h]" << std::endl;
	ofstat << "\\centering \\begin{tabular}{|c|";
	for(size_t i = 0; i < count_column; ++i) { 
		ofstat << "c|";
	} 
	ofstat << "}" << std::endl;
	ofstat << "\\hline" << std::endl;	
} 

void writer::Wstats::print_close_table(bool flag) { 
	if (flag) { 
		ofstat << "\\hline" << std::endl;
	} 
	ofstat << "\\end{tabular}" << std::endl;
	ofstat << "\\end{table}" << std::endl;
	ofstat << std::endl;
} 


