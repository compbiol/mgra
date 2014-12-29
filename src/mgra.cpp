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

#include "io/path_helper.hpp"
#include "logger/logger.hpp"
#include "logger/log_writers.hpp"

#include "reader.h"
#include "RecoveredInfo.hpp"

#include "tclap/CmdLine.h"

bool organize_output_directory(std::string const & path, bool is_debug) { 
  auto creater_lambda = [](std::string const & directory) -> bool {
    if (path::check_existence(directory)) {
      if (path::FileExists(directory)) return false; 
    } else {
      return path::make_dir(directory);
    }
    return true; 
  };

  std::string genomes_dir = path::append_path(path, "genomes"); 
  std::string transformation_dir = path::append_path(path, "transformations");
  bool res = creater_lambda(genomes_dir) && creater_lambda(transformation_dir);

  if (is_debug) {
    std::string debug_dir =  path::append_path(path, "debug"); 
    res = res && creater_lambda(debug_dir);
  } 

  return res;
}

void create_console_logger(std::string const & file_path, std::string const & log_filename) {
  using namespace logging;
  /*
   * Need to create proporties file for logger and use it in build. 
   */
  std::string log_file = path::append_path(file_path, log_filename);
  
  logger *lg = create_logger();
  lg->add_writer(std::make_shared<file_writer>(log_file));
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int argc, char* argv[]) {
  std::string const VERSION(std::to_string(MGRA_VERSION_MAJOR) + "." + std::to_string(MGRA_VERSION_MINOR) + "."
    + std::to_string(MGRA_VERSION_PATCH));
  std::string const LOGGER_FILENAME = "mgra.log";

  /*Reading flags*/
  try {  

  TCLAP::CmdLine cmd("MGRA (Multiple Genome Rearrangements & Ancestors) (c) 2008-2014 by Pavel Avdeyev, Shuai Jiang, Max Alekseyev. Distributed under GNU GENERAL PUBLIC LICENSE license.", ' ', VERSION);

  TCLAP::ValueArg<std::string> path_to_cfg_file_arg("c",
    "config",
    "Input configure file",
    true,
    "",
    "filename", 
    cmd);

  TCLAP::ValueArg<std::string> path_to_blocks_grimm_file_arg("g",
    "grimm_genomes",
    "Input file contains genomes in GRIMM format",
    true,
    "", 
    "filename");

  TCLAP::ValueArg<std::string> path_to_blocks_infercars_file_arg("i",
    "infercars_genomes",
    "Input file contains genomes in InferCARs format",
    true, 
    "", 
    "filename");

  cmd.xorAdd( path_to_blocks_grimm_file_arg, path_to_blocks_infercars_file_arg );
  
  TCLAP::ValueArg<std::string> output_arg("o",
    "output",
    "Output directory",
    true,
    "",
    "dirname", 
    cmd);
  
  TCLAP::SwitchArg debug_arg("d",
    "debug",
    "Switch on debug output",
    cmd,
    false);
  
  cmd.parse(argc, argv); 

  /*Check paths for cfg file, block file and output directory*/
  path::CheckFileExistenceFATAL(path_to_cfg_file_arg.getValue());
  if (path_to_blocks_grimm_file_arg.isSet()) { 
    path::CheckFileExistenceFATAL(path_to_blocks_grimm_file_arg.getValue());
  } else if (path_to_blocks_infercars_file_arg.isSet()) { 
    path::CheckFileExistenceFATAL(path_to_blocks_infercars_file_arg.getValue());
  }
  
  std::string out_path_directory = path::make_full_path(output_arg.getValue());
  if (path::check_existence(out_path_directory)) {
    if (path::FileExists(out_path_directory)) {
      std::cerr << "ERROR: " << out_path_directory << " is not directory" << std::endl;
      return 1;
    }
  } else {
    if (!path::make_dir(out_path_directory)) { 
      std::cerr << "ERROR: Problem to create " << out_path_directory << " directory" << std::endl; 
      return 1;
    }
  }

  if (!organize_output_directory(out_path_directory, debug_arg.getValue())) {
    std::cerr << "ERROR: problem with organize output directory " << out_path_directory << std::endl;
    return 1; 
  }  

  create_console_logger(out_path_directory, LOGGER_FILENAME);

  /*Reading problem configuration and genomes*/
  typedef structure::Genome genome_t;
  typedef structure::Mcolor mcolor_t;
  typedef BreakpointGraph<mcolor_t> graph_t;

  cfg::create_instance(path_to_cfg_file_arg.getValue());
  
  if (cfg::get().get_count_genomes() < 2) {
    ERROR("At least two input genomes required")
    return 1;
  }
  
  std::vector<genome_t> genomes; 
  if (path_to_blocks_infercars_file_arg.isSet()) {
    genomes = reader::read_infercars(path_to_blocks_infercars_file_arg.getValue());
  } else if (path_to_blocks_grimm_file_arg.isSet()) {
    genomes = reader::read_grimm(path_to_blocks_grimm_file_arg.getValue());
  } 

  /*Do job work*/
  for(size_t i = 0; i < genomes.size(); ++i) { 
    std::ostringstream out; 
    out << "Download genome " << cfg::get().get_priority_name(i) << " with " << genomes[i].size() << " blocks.";
    INFO(out.str())
  } 

  INFO("Start build graph")
  std::shared_ptr<graph_t> graph(new graph_t(genomes)); 
  INFO("End build graph")

  { 
    std::ostringstream out; 
    out << "Determine " << graph->count_vec_T_consitent_color() << " \\vec{T}-consistent colors in tree:\n"; 
    for (auto id = graph->cbegin_vec_T_consistent_color(); id != graph->cend_vec_T_consistent_color(); ++id) {
      out << cfg::get().mcolor_to_name(*id) << " ";  
    }
    INFO(out.str());
  } 

  INFO("Start algorithm for convert from breakpoint graph to identity breakpoint graph");
  Algorithm<graph_t> main_algo(graph);
  main_algo.init_writers(out_path_directory, "stage", debug_arg.getValue());
  main_algo.convert_to_identity_bgraph(); 

  if (!cfg::get().is_target_build && !graph->is_identity()) {
    INFO("T-transformation is not complete. Cannot reconstruct genomes.")
    return 1;
  } 
  
  bool consist = graph->check_consistency_graph();
  if (!cfg::get().is_target_build && !consist) {
    INFO("We have problem with edges, corresponding postponed deletions.")
    INFO("If you have indentity breakpoint graph after stages, please contact us.")
    return 1;
  } 

  //FIXME: go from mgra cpp
  INFO("Start to replace cloning to 2-breaks")
  graph->change_history();
  INFO("Finish to replace cloning to 2-breaks")

  writer::Wstats write_stats;
  write_stats.open(out_path_directory, "history_stats.txt");
  write_stats.print_history_statistics(*graph);

  std::shared_ptr<graph_t> new_graph(new graph_t(genomes)); 
  Algorithm<graph_t> alg(new_graph);  
  Algorithm<graph_t>::Balance balance(new_graph);
  balance.do_action();
  
  INFO("Check that history is correct");
  for (auto br = graph->cbegin_2break_history(); br != graph->cend_2break_history(); ++br) {
    new_graph->apply(*br);
  }

  auto bad_edges = graph->get_bad_edges();
  INFO("Start linearization genomes.")
  RecoveredInfo<graph_t> reductant(*graph); 
  INFO("Finsih linearization genomes.")

  INFO("Save transformations in files.")
  if (!cfg::get().is_target_build) {
    size_t i = 0;
    auto recover_transformation = reductant.get_history();
    for (auto im = graph->cbegin_vec_T_consistent_color(); im != graph->cend_vec_T_consistent_color(); ++im, ++i) {
      std::string namefile = cfg::get().mcolor_to_name(*im) + ".trs";
      std::string new_path = path::append_path(out_path_directory, path::append_path("transformations", namefile));
      std::ofstream tr(new_path);
      for(auto const & event : recover_transformation[i]) {
        vertex_t const & p = event.get_vertex(0);
        vertex_t const & q = event.get_vertex(1);
        vertex_t const & x = event.get_vertex(2);
        vertex_t const & y = event.get_vertex(3);
      
      	tr << "(" << p << ", " << q << ") x (" << x << ", " << y << ") " << cfg::get().mcolor_to_name(event.get_mcolor()); 

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

  INFO("Save ancestor genomes in files.")
  writer::Wgenome<genome_t> writer_genome(path::append_path(out_path_directory, "genomes"));
  writer_genome.save_genomes(reductant.get_genomes(), !cfg::get().is_target_build); 

  std::string path_to_logfile = path::append_path(out_path_directory, LOGGER_FILENAME);
  INFO("MGRA log can be found here " << path_to_logfile)
  INFO("Thank you for using MGRA!")

  } catch (TCLAP::ArgException &e) { 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return 1;
  }

  return 0;
}
