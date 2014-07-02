#include "writer/Wstats.h"

void writer::Wstats::print_all_statistics(size_t stage, Statistics<mbgraph_with_history<mcolor_t> >& info, const mbgraph_with_history<mcolor_t>& graph) { 
  if (stage == 0) { 
    ofstat << "Initial graph:" << std::endl;
    ofstat << "... Unique blocks: " << std::to_string(graph.size() / 2) << std::endl;
    ofstat << "... Vertex: " << std::to_string(graph.size()) << std::endl;
  }  else { 
    ofstat << "After Stage " << std::to_string(stage) << " graph:" << std::endl;
  } 

  print_vertex_statistics(info.get_vertex_statistics()); 
 
  print_complete_edges(info.get_complete_edge());
  print_connected_components(graph);

  print_rear_characters(info.get_compl_stat()); 
  print_indel_statistics(info.get_indel_stat());
	
//print_fair_edges(graph, info);
//print_estimated_dist(stage, cfg, graph);

} 

void writer::Wstats::print_complete_edges(std::vector<arc_t> const & edges) { 
  ofstat << "... complete multiedges:";
  for(auto const & edge : edges) {
    ofstat << " " << edge.first << "~" << edge.second;
  } 
  ofstat << "\t(total: " << edges.size() << ")" << std::endl;
} 

void writer::Wstats::print_connected_components(const mbgraph_with_history<mcolor_t>& graph) {
  const std::map<vertex_t, std::set<vertex_t> >& classes = graph.split_on_components(false);
  std::map<size_t, size_t> stx;

  for(const auto &component : classes) {
    ++stx[component.second.size()];
  }

  ofstat << "... connected components:";
  for(const auto &sc : stx) {
    ofstat << " " << sc.first << "^" << sc.second;
  }
  ofstat << std::endl;
}

void writer::Wstats::print_rear_characters(const std::vector<std::string>& info) { 
  ofstat << std::endl;
  ofstat << "% Rearrangement characters:" << std::endl << std::endl;
  print_start_table(6); 
  ofstat << "Multicolors & multiedges & simple vertices & simple multiedges & simple paths+cycles & irreg. multiedges\\\\" << std::endl;
  ofstat << "\\hline" << std::endl;

  for(auto im = info.cbegin(); im != info.cend(); ++im) {
    ofstat << *im << "\\\\" << std::endl;
  }
  print_close_table();	
} 

void writer::Wstats::print_indel_statistics(const std::vector<std::pair<std::pair<mcolor_t, mcolor_t>, std::array<size_t, 3> > >& indels) {
  ofstat << std::endl << "% Insertion/Deletion characters: " << std::endl << std::endl;
  
  print_start_table(5); 
  ofstat << "insert multicolor Q + \\bar{Q} & number edges & number operations & size of split Q & size of split \\bar{Q} \\\\" << std::endl;
  ofstat << "\\hline" << std::endl;

  for (auto line = indels.crbegin(); line != indels.crend(); ++line) { 
    const std::array<size_t, 3>& temp = line->second;
    ofstat << "{" << genome_match::mcolor_to_name(line->first.first) << " + " << genome_match::mcolor_to_name(line->first.second) << "} & " 
	<< temp[0]	<< " & " << (temp[0] * std::min(temp[1], temp[2])) << " & " << temp[1] << " & " << temp[2] << "\\\\" << std::endl;
  } 

  print_close_table();
} 

void writer::Wstats::print_history_statistics(const mbgraph_with_history<mcolor_t>& graph, const edges_t& bad_edges) {
  std::map<mcolor_t, size_t> n2br;
  std::map<mcolor_t, size_t> nins;
  std::map<mcolor_t, size_t> ndel; 
  size_t tbr = 0;
  size_t ins = 0; 
  size_t del = 0;

  for(auto il = graph.crbegin_2break_history(); il != graph.crend_2break_history(); ++il) {
    if (il->get_arc(0).first == Infty || il->get_arc(0).second == Infty 
        || il->get_arc(1).first == Infty || il->get_arc(1).second == Infty) {
      ++n2br[il->get_mcolor()];
      ++tbr;
    } else { 
      const vertex_t& p = il->get_arc(0).first;
      const vertex_t& q = il->get_arc(0).second;
      const vertex_t& x = il->get_arc(1).first;
      const vertex_t& y = il->get_arc(1).second;
      if (p == graph.get_obverse_vertex(x) && bad_edges.defined(p, x)) { 
        ++ndel[il->get_mcolor()];
        ++del;
      } else if (q == graph.get_obverse_vertex(y) && bad_edges.defined(q, y)) {  
        ++ndel[il->get_mcolor()];
        ++del;
      } else if (p == graph.get_obverse_vertex(q) && bad_edges.defined(p, q)) {  
        ++nins[il->get_mcolor()];
        ++ins;
      } else if (x == graph.get_obverse_vertex(y) && bad_edges.defined(y, x)) { 
        ++nins[il->get_mcolor()];
        ++ins;
      } else { 
        ++n2br[il->get_mcolor()];
        ++tbr;
      }
    }
  }

  ofstat << std::endl << "Total number of 2-breaks: " << tbr << std::endl;

  for(const auto &event : n2br) {
    ofstat << genome_match::mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;

  ofstat << std::endl << "Total number of insertion events: " << ins << std::endl;
  for(const auto &event : nins) {
    ofstat << genome_match::mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;

  ofstat << std::endl << "Total number of deletion events: " << del << std::endl;
  for(const auto &event : ndel) {
    ofstat << genome_match::mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;

  size_t c_td = 0; 
  std::map<mcolor_t, size_t> ntd;
  for(auto il = graph.begin_tandem_duplication_history(); il != graph.end_tandem_duplication_history(); ++il) {
    ++c_td;
    ++ntd[il->get_mcolor()];
  }

  ofstat << std::endl << "Total number of (reverse) tandem duplication: " << c_td << std::endl;
  for(auto im = ntd.begin(); im != ntd.end(); ++im) {
    ofstat << genome_match::mcolor_to_name(im->first) << "\t" << im->second << std::endl;
  }
	
  ofstat << std::endl;
}

/*void writer::Wstats::print_postponed_deletion_statistics(const std::map<arc_t, Mcolor>& postponed_deletions) {
  ofstat << std::endl << "Total number of postponed deletions: " << postponed_deletions.size() << std::endl;

  std::map<Mcolor, size_t> npdel; 
  for (const auto &deletion: postponed_deletions) {
    ++npdel[deletion.second]; 
  } 

  for(const auto &event : npdel) {
    ofstat << genome_match::mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;
}

void writer::Wstats::print_bad_complete_edges(const mbgraph_with_history<Mcolor>& graph, const std::multimap<arc_t, Mcolor>& insertions) {
  size_t bad_complete_edges = 0;
  std::unordered_set<vertex_t > processed; 
  std::map<Mcolor, size_t> bad_compl; 
  
  for (const auto &a1 : graph) {  
    if (processed.count(a1) != 0) { 
      continue;
    } 
    const vertex_t& a2 = graph.get_obverse_vertex(a1);
    Mularcs<Mcolor> mularcs = graph.get_adjacent_multiedges(a1);

    if (mularcs.get_multicolor(a2) == graph.get_complete_color())  {
      ++bad_complete_edges;
      if (insertions.count(std::make_pair(a1, a2)) != 0) { 
        auto colors = insertions.equal_range(std::make_pair(a1, a2)); 
        for (auto it = colors.first; it != colors.second; ++it) {  
          ++bad_compl[it->second];
        } 
      } else { 
        auto colors = insertions.equal_range(std::make_pair(a2, a1)); 
        for (auto it = colors.first; it != colors.second; ++it) {  
          ++bad_compl[it->second];
        } 
      } 
    }
    processed.insert(a1); 
    processed.insert(a2);
  }

  ofstat << std::endl << "Bad complete edges: " << bad_complete_edges << std::endl;
  for(const auto &event : bad_compl) {
    ofstat << genome_match::mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;
} 
*/

/*void writer::Wstats::print_estimated_dist(size_t stage, const ProblemInstance<Mcolor>& cfg, const mbgraph_with_history<Mcolor>& graph) { 
	ofstat << "% Estimated distances:" << std::endl << std::endl;
	print_start_table(graph.count_local_graphs());
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
				ofstat << genome_dist(graph.get_local_graph(i), graph.get_local_graph(j), graph.get_obverce_graph())[2]; 
			}

			if (j == graph.count_local_graphs() - 1) { 
				ofstat << "\\\\";
			} else { 
				ofstat << " & ";			
			} 
		}
		ofstat << std::endl;
	}
	print_close_table();
} */

/*
void writer::Wstats::print_fair_edges(const mbgraph_with_history<Mcolor>& MBG, Statistics<mbgraph_with_history<Mcolor>>& info) {
	// output H-subgraphs count
	ofstat << std::endl << "% Fair multi-edges count: " << std::endl << std::endl;

	std::map<std::pair<Mcolor, Mcolor>, size_t> Hcount = info.get_Hsubgraph();
	std::list<Mcolor> HCrow; 	//the set of multicolors that appear in H-subraphs

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

	//print_start_table(HCrow.size());
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

	ofstat << "\\\\" << std::endl;
	ofstat << "\\hline \\hline" << std::endl;

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

				ofstat << Hcount[std::make_pair(*Q1, *Q2)] << " "; 
			}

			if (*Q1 == *Q2) {
				ofstat << "$\\star$";
			}
		}
		ofstat << " \\\\" << std::endl;
		ofstat << "\\hline" << std::endl;
	}
	print_close_table();
}
*/
