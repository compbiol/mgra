//
// Created by Nikita Kartashov on 06/04/2015.
//

#include "io/path_helper.hpp"
#include "tclap/CmdLine.h"

#include "command_line_parsing.hpp"

int parse_config_from_command_line(int argc, char** argv) {
  try {
    std::string const VERSION(std::to_string(MGRA_VERSION_MAJOR) + "." + std::to_string(MGRA_VERSION_MINOR) + "." + std::to_string(MGRA_VERSION_PATCH));

    TCLAP::CmdLine cmd("MGRA (Multiple Genome Rearrangements & Ancestors) (c) 2008-2015 by Pavel Avdeyev, Nikita Kartashov, Shuai Jiang, Max Alekseyev. Distributed under GNU GENERAL PUBLIC LICENSE license.",
        ' ', VERSION);

    /*TCLAP::SwitchArg recover_tree_arg("r",
        "recover_tree",
        "Recover phylogenetic tree instead of ancestral genomes",
        cmd,
        false);
    */

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

    TCLAP::SwitchArg saves_arg("s",
        "saves",
        "Switch on saves output",
        cmd,
        false);

    cmd.parse(argc, const_cast<const char* const*>(argv));

    cfg::get_writable().is_debug = debug_arg.getValue();
    cfg::get_writable().is_saves = saves_arg.getValue();
    cfg::get_writable().how_build = default_algo;
    
    // Disable
    cfg::get_writable().is_recover_tree = false; //recover_tree_arg.getValue();
    cfg::get_writable().recover_tree_statistic = distribution;

    cfg::create_instance(path_to_cfg_file_arg.getValue());

    cfg::get_writable().config_file_path = path_to_cfg_file_arg.getValue();
    if (path_to_blocks_grimm_file_arg.isSet()) {
      cfg::get_writable().blocks_file_path = path_to_blocks_grimm_file_arg.getValue();
      cfg::get_writable().block_file_type = grimm;
    } else if (path_to_blocks_infercars_file_arg.isSet()) {
      cfg::get_writable().blocks_file_path = path_to_blocks_infercars_file_arg.getValue();
      cfg::get_writable().block_file_type = infercars;
    }

    cfg::get_writable().out_path_directory = path::make_full_path(output_arg.getValue());
  
    if (cfg::get().get_count_genomes() < 2) {
      std::cerr << "At least two input genomes required" << std::endl;
      return 1;
    }

    return 0;
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }
}

int validate_application_config() {
  path::CheckFileExistenceFATAL(cfg::get().config_file_path);
  path::CheckFileExistenceFATAL(cfg::get().blocks_file_path);

  if (path::check_existence(cfg::get().out_path_directory)) {
    if (path::FileExists(cfg::get().out_path_directory)) {
      std::cerr << "ERROR: " << cfg::get().out_path_directory << " is not directory" << std::endl;
      return 1;
    }
  } else {
    if (!path::make_dir(cfg::get().out_path_directory)) {
      std::cerr << "ERROR: Problem to create " << cfg::get().out_path_directory << " directory" << std::endl;
      return 1;
    }
  }
  return 0;
}

bool create_dir_if_not_exists(std::string const& directory) {
  if (path::check_existence(directory)) {
    if (path::FileExists(directory)) return false;
  } else {
    return path::make_dir(directory);
  }
  return true;
}

void create_logger_from_config() {
  using namespace logging;
  /*
   * Need to create properties file for logger and use it in build.
   */
  logger* lg = create_logger();
  lg->add_writer(std::make_shared<file_writer>(cfg::get().out_path_to_logger_file));
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

bool organize_output_directory() {
  std::string const LOGGER_FILENAME = "mgra.log";
  std::string const INPUT_DIRNAME = "input";
  std::string const DEBUG_DIRNAME = "debug";
  std::string const SAVES_DIRNAME = "saves";
  std::string const GENOMES_DIRNAME = "genomes";
  std::string const TRANS_DIRNAME = "transformations";

  // Init path to different directories.
  cfg::get_writable().out_path_to_logger_file = path::append_path(cfg::get_writable().out_path_directory, LOGGER_FILENAME);
  cfg::get_writable().out_path_to_input_dir = path::append_path(cfg::get().out_path_directory, INPUT_DIRNAME);
  cfg::get_writable().out_path_to_debug_dir = path::append_path(cfg::get().out_path_directory, DEBUG_DIRNAME);
  cfg::get_writable().out_path_to_saves_dir = path::append_path(cfg::get().out_path_directory, SAVES_DIRNAME);
  cfg::get_writable().out_path_to_genomes_dir = path::append_path(cfg::get().out_path_directory, GENOMES_DIRNAME);
  cfg::get_writable().out_path_to_transfomations_dir = path::append_path(cfg::get().out_path_directory, TRANS_DIRNAME);
  cfg::get_writable().trees_path = path::append_path(cfg::get().out_path_directory, "trees");
  cfg::get_writable().tree_summary_path = path::append_path(cfg::get().trees_path, "summary.newick");

  // Create different directories.
  bool result = create_dir_if_not_exists(cfg::get().out_path_to_genomes_dir) &&
                create_dir_if_not_exists(cfg::get().out_path_to_transfomations_dir) &&
                create_dir_if_not_exists(cfg::get().out_path_to_input_dir) &&
                create_dir_if_not_exists(cfg::get().trees_path);

  if (cfg::get().is_debug) {
    result = result && create_dir_if_not_exists(cfg::get().out_path_to_debug_dir);
  }

  if (cfg::get().is_saves) {
    result = result && create_dir_if_not_exists(cfg::get().out_path_to_saves_dir);
  }

  // Copy 
  path::copy_file(cfg::get().blocks_file_path, path::append_path(cfg::get().out_path_to_input_dir, "blocks.txt"));
  cfg::get_writable().blocks_file_path = path::append_path(cfg::get().out_path_to_input_dir, "blocks.txt");
  path::copy_file(cfg::get_writable().config_file_path, path::append_path(cfg::get().out_path_to_input_dir, "config.txt"));
  cfg::get_writable().config_file_path = path::append_path(cfg::get().out_path_to_input_dir, "config.txt");
  return result;
}