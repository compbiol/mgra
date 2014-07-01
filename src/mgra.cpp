/* 
** Module: MGRA 2.0 main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008 - 2014 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/
#include "algo/Algorithm.h"
#include "writer/Wgenome.h"

#include "reader.h"
#include "RecoveredGenomes.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/exception.hpp>


std::vector<std::string> genome_match::number_to_genome;
genome_match::gen2num genome_match::genome_to_number;   

void tell_root_besides(mbgraph_with_history<structure::Mcolor> const & graph) {
  // tell where the root resides
  std::clog << "the root resides in between:";

  std::set<structure::Mcolor> T(graph.cbegin_T_consistent_color(), graph.cend_T_consistent_color()); 

  for (auto it = T.begin(); it != T.end(); ++it) {
    for (auto jt = it; ++jt != T.end(); ) {
      structure::Mcolor C(*it, *jt, structure::Mcolor::Intersection);
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

int main(int argc, char* argv[]) {
  /*reading flags*/
  std::string const VERSION("2.2");

  std::string cfg_file;
  std::string type; 
  std::string blocks_file;
  std::string out_path; 
  std::string graph_name;
  std::string colorsheme;

  namespace po = boost::program_options;
  po::options_description desc("Options"); 

  desc.add_options()
    ("help,h", "Show help message")
    ("version,v", "Print version program")
    ("config,c", po::value<std::string>(&cfg_file), "Input configure file")
    ("format,f", po::value<std::string>(&type)->required(), "Input format file in genomes file")
    ("genomes,g", po::value<std::string>(&blocks_file)->required(), "Input file contains genomes")
    ("output_dir,o", po::value<std::string>(&out_path)->required(), "Output directory")
    ("filename", po::value<std::string>(&graph_name)->default_value("stage"), "File name for output breakpoint graph")
    ("colorsheme", po::value<std::string>(), "Colorsheme, which used in output breakpoint graph")
    ;

  po::positional_options_description p;
  p.add("config", -1);

  po::variables_map vm; 
  try { 
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm); 

    if ( vm.count("help") ) {
      std::cout << "MGRA (Multiple Genome Rearrangements & Ancestors) version " << VERSION << std::endl;
      std::cout << "(c) 2008-2014 by Shuai Jiang, Pavel Avdeyev, Max Alekseyev" << std::endl;
      std::cout << "Distributed under GNU GENERAL PUBLIC LICENSE license." << std::endl;
      std::cout << std::endl;
      std::cout << desc << std::endl;
      return 0;
    }

    if ( vm.count("version") ) {
      std::cout << "MGRA (Multiple Genome Rearrangements & Ancestors) version " << VERSION << std::endl;
      return 0;
    }

    po::notify(vm);

    if ( !vm.count("config") ) { 
      throw po::required_option("--config");
    }
    
  } catch (po::required_option& e) {
    std::cerr << "ERROR: " << e.what() << std::endl; 
    std::cerr << desc << std::endl;
    return 1;
  } catch (po::error& e) { 
    std::cerr << "ERROR: " << e.what() << std::endl; 
    return 1;
  }

  if ( vm.count("colorsheme") ) { 
    colorsheme = vm["colorsheme"].as<std::string>();
  }

/*
  typedef structure::Genome genome_t;
  
  //Reading problem configuration
  ProblemInstance<structure::Mcolor> cfg(reader::read_cfg_file(name_cfg_file)); 

  std::vector<genome_t> genomes = reader::read_genomes(cfg);
  
  genome_match::init_name_genomes(cfg, genomes); //FIXME: IT'S DEBUG

  for(size_t i = 0; i < genomes.size(); ++i) { 
    std::clog << "Genome " << cfg.get_priority_name(i) << " blocks: " << genomes[i].size() << std::endl;
  } 

  std::shared_ptr<mbgraph_with_history<structure::Mcolor> > graph(new mbgraph_with_history<structure::Mcolor>(genomes, cfg)); 

  std::clog << "vecT-consistent colors: " << graph->count_vec_T_consitent_color() << std::endl;
  for (auto id = graph->cbegin_T_consistent_color(); id != graph->cend_T_consistent_color(); ++id) {
    std::clog << "\t" << genome_match::mcolor_to_name(*id);  ///FIXME: CHANGE
  }
  std::clog << std::endl;

  tell_root_besides(*graph); 	

  if (cfg.is_reconstructed_trees()) {
    std::clog << "Start algorithm for reconstructed trees" << std::endl;
  } 

  std::clog << "Start algorithm for convert from breakpoint graph to identity breakpoint graph" << std::endl;

  Algorithm<mbgraph_with_history<structure::Mcolor> > main_algo(graph, cfg.get_size_component_in_brutforce(), cfg.get_max_number_of_split());
  main_algo.convert_to_identity_bgraph(cfg); 

  if (cfg.get_target().empty()) {
    for(auto it = graph->cbegin_local_graphs(); it != graph->cend_local_graphs() - 1; ++it) { 
      if (*it != *(it + 1)) {
	std::clog << "T-transformation is not complete. Cannot reconstruct genomes." << std::endl; 
        exit(1);
      }
    }
  } 

  auto bad_edges = main_algo.get_bad_edges();

  std::clog << "Start reconstruct genomes." << std::endl;
  RecoveredGenomes<mbgraph_with_history<structure::Mcolor> > reductant(*graph, cfg, bad_edges); 

  std::clog << "Save history in files." << std::endl;
  if (cfg.get_target().empty()) {
    size_t i = 0;
    auto recover_transformation = reductant.get_history();
    for (auto im = graph->cbegin_T_consistent_color(); im != graph->cend_T_consistent_color(); ++im, ++i) {
      std::ofstream tr((cfg.mcolor_to_name(*im) + ".trs").c_str());
      for(auto const & event : recover_transformation[i]) {
        vertex_t const & p = event.get_arc(0).first;
        vertex_t const & q = event.get_arc(0).second;
        vertex_t const & x = event.get_arc(1).first;
        vertex_t const & y = event.get_arc(1).second;
      
	tr << "(" << p << ", " << q << ") x (" << x << ", " << y << ") " << genome_match::mcolor_to_name(event.get_mcolor()); 

	if (p != Infty && q != Infty && x != Infty && y != Infty) { 
	  if (p == graph->get_obverse_vertex(x) && bad_edges.defined(p, x)) { 
	    tr << " # deletion"; 
          } else if (q == graph->get_obverse_vertex(y) && bad_edges.defined(q, y)) {
	    tr << " # deletion"; 
	  } else if (p == graph->get_obverse_vertex(q) && bad_edges.defined(p, q)) { 
	    tr << " # insertion"; 
          } else if (x == graph->get_obverse_vertex(y) && bad_edges.defined(x, y)) {
	    tr << " # insertion"; 
	  } 
	}
	tr << std::endl;  
      }
      tr.close(); 
    } 
  }

  std::clog << "Save ancestor genomes in files." << std::endl; 
  writer::Wgenome<genome_t> writer_genome;
  writer_genome.save_genomes(reductant.get_genomes(), cfg.get_target().empty()); 
*/ 
  return 0;
}

