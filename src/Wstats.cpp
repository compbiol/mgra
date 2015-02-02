#include "writer/Wstats.h"

/*void writer::Wstats::print_all_statistics(size_t stage, Statistics<GraphPack<mcolor_t> >& info, const GraphPack<mcolor_t>& graph) { 
  if (stage == 0) { 
    ofstat << "Initial graph:" << std::endl;
    ofstat << "... Unique blocks: " << std::to_string(graph.graph.size() / 2) << std::endl;
    ofstat << "... Vertex: " << std::to_string(graph.graph.size()) << std::endl;
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
} */

void writer::Wstats::print_complete_edges(std::vector<edge_t> const & edges) { 
  ofstat << "... complete multiedges:";
  for(auto const & edge : edges) {
    ofstat << " " << edge.first << "~" << edge.second;
  } 
  ofstat << "\t(total: " << edges.size() << ")" << std::endl;
} 

void writer::Wstats::print_connected_components(const GraphPack<mcolor_t>& graph) {
  auto components = graph.split_on_components(false); 
  std::map<vertex_t, std::set<vertex_t> > const & classes = components.get_eclasses<std::set<vertex_t> >();
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
    ofstat << "{" << cfg::get().mcolor_to_name(line->first.first) << " + " << cfg::get().mcolor_to_name(line->first.second) << "} & " 
	<< temp[0]	<< " & " << (temp[0] * std::min(temp[1], temp[2])) << " & " << temp[1] << " & " << temp[2] << "\\\\" << std::endl;
  } 

  print_close_table();
} 

void writer::Wstats::print_history_statistics(GraphPack<mcolor_t> const & graph) {
  auto const & bad_edges = graph.get_bad_edges();
  std::map<mcolor_t, size_t> n2br;
  std::map<mcolor_t, size_t> nins;
  std::map<mcolor_t, size_t> ndel; 
  size_t tbr = 0;
  size_t ins = 0; 
  size_t del = 0;

  for(auto il = graph.history.rbegin(); il != graph.history.rend(); ++il) {
    vertex_t const & p = il->get_vertex(0);
    vertex_t const & q = il->get_vertex(1);
    vertex_t const & x = il->get_vertex(2);
    vertex_t const & y = il->get_vertex(3);

    if (p == Infty || q == Infty || x == Infty || y == Infty) {
      ++n2br[il->get_mcolor()];
      ++tbr;
    } else { 
      if (p == graph.graph.get_obverse_vertex(x) && bad_edges.defined(p, x)) { 
        ++ndel[il->get_mcolor()];
        ++del;
      } else if (q == graph.graph.get_obverse_vertex(y) && bad_edges.defined(q, y)) {  
        ++ndel[il->get_mcolor()];
        ++del;
      } else if (p == graph.graph.get_obverse_vertex(q) && bad_edges.defined(p, q)) {  
        ++nins[il->get_mcolor()];
        ++ins;
      } else if (x == graph.graph.get_obverse_vertex(y) && bad_edges.defined(y, x)) { 
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
    ofstat << cfg::get().mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;

  ofstat << std::endl << "Total number of insertion events: " << ins << std::endl;
  for(const auto &event : nins) {
    ofstat << cfg::get().mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;

  ofstat << std::endl << "Total number of deletion events: " << del << std::endl;
  for(const auto &event : ndel) {
    ofstat << cfg::get().mcolor_to_name(event.first) << "\t" << event.second << std::endl;
  }
  ofstat << std::endl;

  /*size_t c_td = 0; 
  std::map<mcolor_t, size_t> ntd;
  for(auto il = graph.begin_tandem_duplication_history(); il != graph.end_tandem_duplication_history(); ++il) {
    ++c_td;
    ++ntd[il->get_mcolor()];
  }*/

  //ofstat << std::endl << "Total number of (reverse) tandem duplication: " << c_td << std::endl;
  /*for(auto im = ntd.begin(); im != ntd.end(); ++im) {
    ofstat << genome_match::mcolor_to_name(im->first) << "\t" << im->second << std::endl;
  }
	
  ofstat << std::endl;*/
}

/*void writer::Wstats::print_postponed_deletion_statistics(const std::map<edge_t, Mcolor>& postponed_deletions) {
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

*/
