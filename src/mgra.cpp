/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
*/

#include "tclap/CmdLine.h"

#include "reader.h"

#include "algo/Algorithms.hpp"
#include "algo/recover_tree/bruteforce_recover_tree_algorithm.hpp"

#include "io/path_helper.hpp"
#include "logger/logger.hpp"
#include "logger/log_writers.hpp"

#include "RecoveredInfo.hpp"

#include "writer/Wgenome.h"
#include "writer/Wtransform.hpp"

bool organize_output_directory(std::string const &path, bool is_debug) {
  auto creater_lambda = [](std::string const &directory) -> bool {
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
    std::string debug_dir = path::append_path(path, "debug");
    res = res && creater_lambda(debug_dir);
  }

  return res;
}

void create_logger(std::string const &file_path, std::string const &log_filename) {
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

int parse_configure_file(
    TCLAP::SwitchArg &debug_arg,
    std::string const &out_path_directory,
    TCLAP::SwitchArg &target_arg,
    TCLAP::ValueArg<std::string> &path_to_cfg_file_arg) {
  INFO("Parse configure file");
  cfg::get_writable().is_debug = debug_arg.getValue();
  cfg::get_writable().out_path_to_debug_dir = path::append_path(out_path_directory, "debug");
  if (target_arg.getValue()) {
    cfg::get_writable().how_build = target_algo;
  } else {
    cfg::get_writable().how_build = default_algo;
  }

  cfg::create_instance(path_to_cfg_file_arg.getValue());

  if (cfg::get().get_count_genomes() < 2) {
    ERROR("At least two input genomes required")
    return 1;
  }
  return 0;
}

int main(int argc, char **argv) {
  std::string const VERSION(std::to_string(MGRA_VERSION_MAJOR) + "." + std::to_string(MGRA_VERSION_MINOR) + "."
      + std::to_string(MGRA_VERSION_PATCH));
  std::string const LOGGER_FILENAME = "mgra.log";

//  Reading flags
  try {
    //TODO: Move all TCLAP parsing into separate fucntion + write all data straight to config && validate it

    TCLAP::CmdLine cmd("MGRA (Multiple Genome Rearrangements & Ancestors) (c) 2008-2015 by Pavel Avdeyev, Shuai Jiang, Max Alekseyev. Distributed under GNU GENERAL PUBLIC LICENSE license.", ' ', VERSION);

    TCLAP::SwitchArg target_arg("t",
        "target",
        "Switch on target reconstruction algorithm",
        cmd,
        false);

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

    cmd.xorAdd(path_to_blocks_grimm_file_arg, path_to_blocks_infercars_file_arg);

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

    TCLAP::SwitchArg recover_tree_arg("r",
        "recover_tree",
        "Recover phylogenetic tree instead of ancestral genomes",
        cmd,
        false);

    cmd.parse(argc, const_cast<const char *const *>(argv));

//    Check paths for cfg file, block file and output directory
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

    create_logger(out_path_directory, LOGGER_FILENAME);

//    Reading problem configuration and genomes
    using genome_t = structure::Genome;
    using mcolor_t = structure::Mcolor;
    using graph_pack_t = GraphPack<mcolor_t>;

    if (parse_configure_file(debug_arg, out_path_directory, target_arg, path_to_cfg_file_arg)) {
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
      out << "Download genome " << cfg::get().get_priority_name(i) << " with " << genomes[i].size() << " blocks.";
      INFO(out.str())
    }

//    Do job
    INFO("Start build graph")
    graph_pack_t graph_pack(genomes);
    INFO("End build graph")

    cfg::get_writable().is_recover_tree = recover_tree_arg.getValue();

    if (cfg::get().is_recover_tree) {
      INFO("Starting tree reovery")
      //Recover tree here
      graph_pack.update_graph_statistics();

      typename algo::RecoverTreeAlgorithm<graph_pack_t>::algo_ptr recover_tree_algoritm(
          new algo::BruteforceRecoverTreeAlgorithm<graph_pack_t>(graph_pack));

      auto result_tree = recover_tree_algoritm->recover_tree();

      std::cout << result_tree;
      //TODO: perform dumping
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
        writer::Wtransformation<graph_pack_t> writer_transform(out_path_directory, graph_pack);
        auto recover_transformations = reductant.get_history();
        for (auto const &transformation : recover_transformations) {
          writer_transform.save_transformation(transformation.first, transformation.second);
          writer_transform.save_reverse_transformation(transformation.first, transformation.second);
        }
      }

      INFO("Save ancestor genomes in files.")
      writer::Wgenome<genome_t> writer_genome(path::append_path(out_path_directory, "genomes"));
      writer_genome.save_genomes(reductant.get_genomes(), (cfg::get().how_build == default_algo));
    }
    std::string path_to_logfile = path::append_path(out_path_directory, LOGGER_FILENAME);
    INFO("MGRA log can be found here " << path_to_logfile)
    INFO("Thank you for using MGRA!")

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
