/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
*/

#include "command_line_parsing.hpp"

#include "algo/Algorithms.hpp"
#include "algo/recover_tree/recover_tree_task.hpp"

#include "logger/log_writers.hpp"

#include "writer/txt_genome.hpp"
#include "writer/txt_transform.hpp"


<<<<<<< HEAD
  if (is_debug) {
    std::string debug_dir =  path::append_path(path, "debug"); 
    res = res && creater_lambda(debug_dir);
  } 

  return res;
}

void create_logger(std::string const & file_path, std::string const & log_filename) {
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

  TCLAP::CmdLine cmd("MGRA (Multiple Genome Rearrangements & Ancestors) (c) 2008-2015 by Pavel Avdeyev, Shuai Jiang, Max Alekseyev. Distributed under GNU GENERAL PUBLIC LICENSE license.", ' ', VERSION);
  
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
=======
int main(int argc, char** argv) {
  if (parse_config_from_command_line(argc, argv)) {
    std::cerr << "ERROR: error while parsing command line arguments" << std::endl;
    return 1;
>>>>>>> recover-tree
  }

  if (validate_application_config()) {
    std::cerr << "ERROR: error while validating config" << std::endl;
    return 1;
  }

  if (!organize_output_directory()) {
    std::cerr << "ERROR: problem with organize output directory "
        << cfg::get().out_path_directory << std::endl;
    return 1;
  }

  create_logger_from_config();

//    Reading problem configuration and genomes
  using genome_t = structure::Genome;
  using mcolor_t = structure::Mcolor;
  using graph_pack_t = GraphPack<mcolor_t>;

<<<<<<< HEAD
  INFO("Parse configure file")
  cfg::get_writable().is_debug = debug_arg.getValue();
  cfg::get_writable().out_path_to_debug_dir = path::append_path(out_path_directory, "debug");
  cfg::get_writable().how_build = default_algo;
  cfg::create_instance(path_to_cfg_file_arg.getValue());
  if (cfg::get().get_count_genomes() < 2) {
    ERROR("At least two input genomes required")
    return 1;
  }
  
  INFO("Parse genomes file")
  std::vector<genome_t> genomes; 
  if (path_to_blocks_infercars_file_arg.isSet()) {
    genomes = reader::read_infercars(path_to_blocks_infercars_file_arg.getValue());
  } else if (path_to_blocks_grimm_file_arg.isSet()) {
    genomes = reader::read_grimm(path_to_blocks_grimm_file_arg.getValue());
  } 

  for (size_t i = 0; i < genomes.size(); ++i) { 
    std::ostringstream out; 
=======
  INFO("Parse genomes file")
  std::vector<genome_t> genomes;
  if (cfg::get().block_file_type == infercars) {
    genomes = reader::read_infercars(cfg::get().blocks_file_path);
  } else if (cfg::get().block_file_type == grimm) {
    genomes = reader::read_grimm(cfg::get().blocks_file_path);
  }

  for (size_t i = 0; i < genomes.size(); ++i) {
    std::ostringstream out;
>>>>>>> recover-tree
    out << "Download genome " << cfg::get().get_priority_name(i) << " with " << genomes[i].size() << " blocks.";
    INFO(out.str())
  }

//    Do job
  INFO("Start build graph")
  graph_pack_t graph_pack(genomes);
  INFO("End build graph")


<<<<<<< HEAD
  boost::optional<algo::RecoveredInformation<graph_pack_t>::AncestorInformation> result; 
  if (cfg::get().how_build == default_algo) { 
    result = main_algorithm(graph_pack);
  } else if (cfg::get().how_build == wgd_algo) { 
    result = wgd_algorithm(graph_pack);
  } 

  /*Save different output information in files*/
  if (result) { 
    using ancestor_information_t = algo::RecoveredInformation<graph_pack_t>::AncestorInformation;
    ancestor_information_t info = *result; 

    if (!info.transformations.empty()) {
      INFO("Save transformations in files.")

      writer::TXT_transformation<graph_pack_t> writer_transform(out_path_directory, graph_pack); 

      for (auto const & transformation : info.transformations) { 
        writer_transform.save_transformation(transformation.first, transformation.second);
        writer_transform.save_reverse_transformation(transformation.first, transformation.second);
      }    
    } 

    if (!info.genomes.empty()) {
      INFO("Save ancestor genomes in files.")
      writer::TXT_genome<genome_t> writer_genome(path::append_path(out_path_directory, "genomes"));
      writer_genome.save_genomes(info.genomes); 
    } 
  } 
  
  std::string path_to_logfile = path::append_path(out_path_directory, LOGGER_FILENAME);
  INFO("MGRA log can be found here " << path_to_logfile)
  INFO("Thank you for using MGRA!")
  if (!result) { 
    return 1;
  }
=======
  if (cfg::get().is_recover_tree) {
    algo::recover_tree_task(graph_pack);
  } else {
    {
      std::ostringstream out;
      out << "Determine " << graph_pack.multicolors.count_vec_T_consitent_color() << " \\vec{T}-consistent colors in tree:\n";
      for (auto id = graph_pack.multicolors.cbegin_vec_T_consistent_color(); id != graph_pack.multicolors.cend_vec_T_consistent_color(); ++id) {
        out << cfg::get().mcolor_to_name(*id) << " ";
      }
      INFO(out.str())
    }

    bool result = false; //wgd_algorithm(graph_pack);
    if (cfg::get().how_build == default_algo) {
      result = main_algorithm(graph_pack);
    } else if (cfg::get().how_build == target_algo) {
    }

    if (!result) {
      return 1;
    }

    INFO("Start linearization genomes.")
    RecoveredInfo<graph_pack_t> reductant(graph_pack);
    INFO("Finish linearization genomes.")

    INFO("Save transformations in files.")
    if (cfg::get().how_build == default_algo) {
      writer::Wtransformation<graph_pack_t> writer_transform(cfg::get().out_path_directory, graph_pack);
      auto recover_transformations = reductant.get_history();
      for (auto const& transformation : recover_transformations) {
        writer_transform.save_transformation(transformation.first, transformation.second);
        writer_transform.save_reverse_transformation(transformation.first, transformation.second);
      }
    }
>>>>>>> recover-tree

    INFO("Save ancestor genomes in files.")
    writer::Wgenome<genome_t> writer_genome(cfg::get().genomes_path);
    writer_genome.save_genomes(reductant.get_genomes(), (cfg::get().how_build == default_algo));
  }
  INFO("MGRA log can be found here " << cfg::get().logger_path)
  INFO("Thank you for using MGRA!")

  return 0;
}
