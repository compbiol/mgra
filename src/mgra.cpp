/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
*/

#include "algo/Algorithm.h"
#include "writer/Wgenome.h"

#include "reader.h"
#include "RecoveredGenomes.h"
#include "RecoveredInfo.hpp"

#include <boost/program_options.hpp>

std::vector<std::string> genome_match::number_to_genome;
genome_match::gen2num genome_match::genome_to_number;   

void tell_root_besides(BreakpointGraph<structure::Mcolor> const & graph) {
  // tell where the root resides
  std::clog << "the root resides in between:";

  std::set<structure::Mcolor> colors(graph.cbegin_vec_T_consistent_color(), graph.cend_vec_T_consistent_color()); 

  for (auto it = colors.begin(); it != colors.end(); ++it) {
    for (auto jt = it; ++jt != colors.end(); ) {
      structure::Mcolor C(*it, *jt, structure::Mcolor::Intersection);
      if (C.size() == it->size()) {
        colors.erase(it++);
        jt = it;
        continue;
      }
      if (C.size() == jt->size()) {
        colors.erase(jt++);
        --jt;
      }
    }
    std::clog << " " << genome_match::mcolor_to_name(*it);
  }
  std::clog << std::endl;
}

bool organize_output_directory( fs::path const & path ) { 
  auto creater_lambda = [](fs::path const & directory) -> bool {
    sys::error_code error; 
    if ( fs::exists(directory, error) ) {
      if ( !fs::is_directory(directory, error) ) {
        return false;
      }
    } else {
      bool flag = fs::create_directory(directory, error); 

      if ( !flag || error != 0 ) { 
        return false;
      }
    }
    return true; 
  };

  fs::path debug_dir = path / "debug"; 
  fs::path genomes_dir = path / "genomes"; 
  fs::path transformation_dir = path / "transformations";

  return creater_lambda(debug_dir) && creater_lambda(genomes_dir) && creater_lambda(transformation_dir);
}

std::vector<BreakpointGraph<structure::Mcolor>::twobreak_t> save_total_history(fs::path const & path, BreakpointGraph<structure::Mcolor> const & graph, 
        ProblemInstance<structure::Mcolor> const & cfg) {
  typedef BreakpointGraph<structure::Mcolor>::twobreak_t twobreak_t;
  
  std::vector<twobreak_t> transform;

  std::cerr << "::::Save total history:::: " << std::endl;
  for(auto il = graph.cbegin_2break_history(); il != graph.cend_2break_history(); ++il) {
    auto twobreak = il->get_canonical_twobreak();
    transform.push_back(twobreak); 
  } 

  for (size_t i = 0; i < transform.size() - 1; ++i) { 
    for (size_t j = 0; j < transform.size() - i - 1; ++j) { 
      if ((transform[j] > transform[j + 1]) && transform[j].is_independent(transform[j + 1])) {
        std::swap(transform[j], transform[j + 1]); 
      }
    }
  }

  fs::ofstream tr(path / "sorted_full_history.txt");
  
  for (auto il = transform.cbegin(); il != transform.cend(); ++il) { 
    vertex_t const & p = il->get_vertex(0);
    vertex_t const & q = il->get_vertex(1);
    vertex_t const & x = il->get_vertex(2);
    vertex_t const & y = il->get_vertex(3);

    tr << p << " " << q << " " << x << " " << y << " " << cfg.mcolor_to_name(il->get_mcolor()) << std::endl; 
  }

  tr.close();

  return transform;
}

int main(int argc, char* argv[]) {
  std::string const VERSION(std::to_string(MGRA_VERSION_MAJOR) + "." + std::to_string(MGRA_VERSION_MINOR) + "."
    + std::to_string(MGRA_VERSION_PATCH));

  std::string cfg_file;
  std::string type; 
  std::string blocks_file;
  std::string out_path; 
  bool debug = false; 

  /*Reading flags*/
  namespace po = boost::program_options;
  po::options_description desc("Options"); 

  desc.add_options()
    ("help,h", "Show help message")
    ("version,v", "Print version program")
    ("config,c", po::value<std::string>(&cfg_file), "Input configure file")
    ("format,f", po::value<std::string>(&type)->required(), "Input format file in genomes file")
    ("genomes,g", po::value<std::string>(&blocks_file)->required(), "Input file contains genomes")
    ("output_dir,o", po::value<std::string>(&out_path)->required(), "Output directory")
    ("debug,d", "Switch on debug output")
    ;

  po::positional_options_description p;
  p.add("config", -1);

  po::variables_map vm; 
  try { 
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm); 

    if ( vm.count("help") ) {
      std::cout << "MGRA (Multiple Genome Rearrangements & Ancestors) version " << VERSION << std::endl;
      std::cout << "(c) 2008-2014 by Pavel Avdeyev, Shuai Jiang, Max Alekseyev" << std::endl;
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

  if ( vm.count("debug") ) { 
    debug = true; 
  }

  /*Check paths for cfg file, block file and output directory*/
  sys::error_code error;
  fs::path out_path_directory(out_path);
  fs::path in_path_cfg(cfg_file);
  fs::path in_path_blocks(blocks_file);

  if ( !fs::exists(in_path_cfg, error) || !fs::is_regular_file(in_path_cfg, error) ) { 
    std::cerr << "ERROR: Problem with configuration file: " << cfg_file << std::endl;
  }

  if ( !fs::exists(in_path_blocks, error) || !fs::is_regular_file(in_path_blocks, error) ) { 
    std::cerr << "ERROR: Problem with genomes file: " << blocks_file << std::endl;
  }

  if ( fs::exists(out_path_directory, error) ) {
    if ( !fs::is_directory(out_path_directory, error) ) {
      std::cerr << "ERROR: " << out_path << " is not directory" << std::endl;
      return 1;
    }
  } else {
    bool flag = fs::create_directory(out_path_directory, error); 

    if ( !flag || error != 0 ) { 
      std::cerr << "ERROR: Problem to create " << out_path << " directory" << std::endl; 
      return 1;
    }
  }

  out_path_directory = fs::canonical(out_path_directory);
  bool is_create = organize_output_directory(out_path_directory); 
  if ( !is_create ) {
    std::cerr << "ERROR: problem with organize output directory " << out_path_directory << std::endl;
    return 1; 
  }  

  /*Reading problem configuration and genomes*/
  typedef structure::Genome genome_t;
  typedef structure::Mcolor mcolor_t;
  typedef BreakpointGraph<mcolor_t> graph_t;

  ProblemInstance<mcolor_t> cfg(reader::read_cfg_file(in_path_cfg)); 

  if (cfg.get_count_genomes() < 2) {
    std::cerr << "ERROR: at least two input genomes required" << std::endl;
    return 1;
  }
  
  std::vector<genome_t> genomes; 
  if (type == "infercars") {
    genomes = reader::read_infercars(cfg, in_path_blocks);
  } else if (type == "grimm") {
    genomes = reader::read_grimm(cfg, in_path_blocks);
  } else {
    std::cerr << "ERROR: unknown synteny blocks format" << type << std::endl;
    return 1;
  }

  /*Do job work*/
  //if ( debug ) { 
  genome_match::init_name_genomes(cfg, genomes); //FIXME: IT'S DEBUG
  //}

  for(size_t i = 0; i < genomes.size(); ++i) { 
    std::clog << "Genome " << cfg.get_priority_name(i) << " blocks: " << genomes[i].size() << std::endl;
  } 


  std::shared_ptr<graph_t> graph(new graph_t(genomes, cfg)); 

  std::clog << "vecT-consistent colors: " << graph->count_vec_T_consitent_color() << std::endl;
  for (auto id = graph->cbegin_vec_T_consistent_color(); id != graph->cend_vec_T_consistent_color(); ++id) {
    std::clog << "\t" << genome_match::mcolor_to_name(*id);  ///FIXME: CHANGE
  }
  std::clog << std::endl;

  tell_root_besides(*graph); 	

  std::clog << "Start algorithm for convert from breakpoint graph to identity breakpoint graph" << std::endl;
  Algorithm<graph_t> main_algo(graph, cfg);
  main_algo.init_writers(out_path_directory, "stage", debug);
  main_algo.convert_to_identity_bgraph(); 

  if (cfg.get_target().empty() && !graph->is_identity()) {
    std::clog << "T-transformation is not complete. Cannot reconstruct genomes." << std::endl; 
    return 1;
  } 
  
  bool consist = graph->check_consistency_graph();
  if (cfg.get_target().empty() && !consist) {
    std::cerr << "We have problem with edges, corresponding postponed deletions." << std::endl;
    std::cerr << "If you have indentity breakpoint graph after stages, please contact us." << std::endl;
    return 1;
  } 

  //FIXME: go to mgra cpp
  std::cerr << "start to change history " << std::endl;

  graph->change_history();
  
  writer::Wstats write_stats;
  write_stats.open(out_path_directory, "history_stats.txt");
  write_stats.print_history_statistics(*graph);

  //std::shared_ptr<graph_t> new_graph(new graph_t(genomes, cfg)); 
  std::shared_ptr<graph_t> new_graph(new graph_t(genomes, cfg)); 
  Algorithm<graph_t> alg(new_graph, cfg);  
  Algorithm<graph_t>::Balance balance(new_graph);
  balance.do_action();
  //writer::Wdots<graph_t, ProblemInstance<mcolor_t> > write_dots(cfg);
  //write_dots.save_dot(*new_graph, cfg, 100);

  for (auto br = graph->cbegin_2break_history(); br != graph->cend_2break_history(); ++br) {
    //std::cerr << br->get_vertex(0) << " " << br->get_vertex(1) << " " << br->get_vertex(2) 
    //<< " " << br->get_vertex(3) << " " << genome_match::mcolor_to_name(br->get_mcolor()) << std::endl;
    /*if (br->get_vertex(0) == "190h" && br->get_vertex(1) == "423t" && br->get_vertex(2) == "120t"
        && br->get_vertex(3) == "946h") { 
      auto color = new_graph->get_all_multicolor_edge("190h", "423t");
      std::cerr << genome_match::mcolor_to_name(color)  << std::endl;
      color = new_graph->get_all_multicolor_edge("120t", "946h");
      std::cerr << genome_match::mcolor_to_name(color)  << std::endl;
    } */
    new_graph->apply(*br);
  }

  save_total_history(out_path_directory, *graph, cfg);

  auto bad_edges = graph->get_bad_edges();
  std::clog << "Start reconstruct genomes." << std::endl;
  RecoveredInfo<graph_t> reductant(*graph, cfg); 
  
  std::clog << "Save history in files." << std::endl;
  if (cfg.get_target().empty()) {
    size_t i = 0;
    auto recover_transformation = reductant.get_history();
    for (auto im = graph->cbegin_vec_T_consistent_color(); im != graph->cend_vec_T_consistent_color(); ++im, ++i) {
      std::string namefile = cfg.mcolor_to_name(*im) + ".trs";
      fs::ofstream tr(out_path_directory / "transformations" / namefile);
      for(auto const & event : recover_transformation[i]) {
        vertex_t const & p = event.get_vertex(0);
        vertex_t const & q = event.get_vertex(1);
        vertex_t const & x = event.get_vertex(2);
        vertex_t const & y = event.get_vertex(3);
      
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
  writer::Wgenome<genome_t> writer_genome(out_path_directory / "genomes");
  writer_genome.save_genomes(reductant.get_genomes(), cfg.get_target().empty()); 

  return 0;
}
