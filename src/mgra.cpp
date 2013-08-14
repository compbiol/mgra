/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008 - 2013 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#include "reader.h"
#include "algo/Algorithm.h"
#include "Wgenome.h"

std::vector<std::string> genome_match::number_to_genome;
genome_match::gen2num genome_match::genome_to_number;   

typedef std::list<TwoBreak<Mcolor> > transform_t;

#include "RecoveredGenomes.h"

void tell_root_besides(const mbgraph_with_history<Mcolor>& graph) {
  // tell where the root resides
  std::clog << "the root resides in between:";

  std::set<Mcolor> T(graph.cbegin_T_consistent_color(), graph.cend_T_consistent_color()); 

  for (auto it = T.begin(); it != T.end(); ++it) {
    for (auto jt = it; ++jt != T.end(); ) {
      Mcolor C(*it, *jt, Mcolor::Intersection);
      if (C.size() == it->size()) {
	T.erase(it++);
	jt = it;
	continue;
      }
      if (C.size() == jt->size()) {
	T.erase(jt++);
	--jt;
      }
    }
    std::clog << " " << genome_match::mcolor_to_name(*it);
  }
  std::clog << std::endl;
}


////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
  std::cout << "MGRA (Multiple Genome Rearrangements & Ancestors) ver. 1.5" << std::endl;
  std::cout << "(c) 2008,12 by Max Alekseyev <maxal@cse.sc.edu>" << std::endl;
  std::cout << "Distributed under GNU GENERAL PUBLIC LICENSE license." << std::endl;
  std::cout << std::endl;

  /*reading flags*/
  std::string name_cfg_file;
  if (argc != 2) {
    std::cerr << "Usage: mgra <ProblemConfiguration>" << std::endl;
    return 1;
  } else { 
    name_cfg_file = argv[1];
  } 

  /*Reading problem configuration*/
  ProblemInstance<Mcolor> PI(reader::read_cfg_file(name_cfg_file)); 

  std::vector<Genome> genomes = reader::read_genomes(PI);
  
  genome_match::init_name_genomes(PI, genomes); //FIXME: IT'S DEBUG

  for(size_t i = 0; i < genomes.size(); ++i) { 
	std::clog << "Genome " << PI.get_priority_name(i) << " blocks: " << genomes[i].size() << std::endl;
  } 

  std::shared_ptr<mbgraph_with_history<Mcolor> > graph(new mbgraph_with_history<Mcolor>(genomes, PI)); 

  std::clog << "vecT-consistent colors: " << graph->count_vec_T_consitent_color() << std::endl;
  for (auto id = graph->cbegin_T_consistent_color(); id != graph->cend_T_consistent_color(); ++id) {
    std::clog << "\t" << genome_match::mcolor_to_name(*id);  ///FIXME: CHANGE
  }
  std::clog << std::endl;

  tell_root_besides(*graph); 	

  Algorithm<mbgraph_with_history<Mcolor> > main_algo(graph);

  main_algo.convert_to_identity_bgraph(PI); 
 
  if (PI.get_target().empty()) {
    for(auto it = graph->cbegin_local_graphs(); it != graph->cend_local_graphs() - 1; ++it) { 
      if (*it != *(it + 1)) {
	std::cout << "T-transformation is not complete. Cannot reconstruct genomes." << std::endl; //FIXME clog
	exit(1);
      }
    }
  } 

  RecoveredGenomes<mbgraph_with_history<Mcolor> > reductant(*graph, PI.get_target());	
	
  if (!PI.get_target().empty()) {
    Genome genome = reductant.get_genomes()[0]; 
    writer::Wgenome<std::set<std::pair<path_t, bool> > > writer_genome; 
    writer_genome.save_genom_in_text_format(genome_match::mcolor_to_name(PI.get_target()), genome, PI.get_target().empty());
  } else {  /* empty target */
    size_t i = 0;
    //std::vector<Genome> genomes = reductant.get_genomes();
    std::vector<transform_t> RT = reductant.get_history();
    for (auto im = graph->cbegin_T_consistent_color(); im != graph->cend_T_consistent_color(); ++im, ++i) {
      Genome genome; 
      reductant.splitchr(reductant.get_graph(i), genome, reductant.pg_empty);

      writer::Wgenome<std::set<std::pair<path_t, bool> > > writer_genome; 
      writer_genome.save_genom_in_text_format(genome_match::mcolor_to_name(*im), genome, PI.get_target().empty());
        
      //writer_genome.save_genom_in_text_format(genome_match::mcolor_to_name(*im), genomes[i], PI.get_target().empty());
        
      std::ofstream tr((genome_match::mcolor_to_name(*im) + ".trs").c_str());
      for(const auto &event : RT[i]) {
	tr << event.get_arc(0).first << " " << event.get_arc(0).second << "\t" 
	   << event.get_arc(1).first << " " << event.get_arc(1).second << "\t" 
  	   << genome_match::mcolor_to_name(event.get_mcolor()) << std::endl;
      }
      tr.close(); 
    }
  }

  return 0;
}

