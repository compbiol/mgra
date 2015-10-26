//
// Created by pavel on 10/21/15.
//

#include "io/path_helper.hpp"
#include "tclap/CmdLine.h"

#include "command_line_parsing.hpp"

int parse_config_from_command_line(int argc, char** argv) {
  try {
    std::string const VERSION(std::to_string(MGRA_VERSION_MAJOR) + "." + std::to_string(MGRA_VERSION_MINOR) + "." + std::to_string(MGRA_VERSION_PATCH));

    TCLAP::CmdLine cmd("Recovered tree method for MGRA (Multiple Genome Rearrangements & Ancestors) (c) 2008-2015 by Pavel Avdeyev, Nikita Kartashov, Max Alekseyev. Distributed under GNU GENERAL PUBLIC LICENSE license.",
        ' ', VERSION);

    TCLAP::ValueArg<std::string> path_to_cfg_file_arg("c",
        "config",
        "Input configure file",
        true,
        "",
        "filename",
        cmd);

    cmd.parse(argc, const_cast<const char* const*>(argv));
    std::ifstream cfg_file(path_to_cfg_file_arg.getValue());
    cfg::create_instance(cfg_file);

    return 0;
  } catch (TCLAP::ArgException& e) {
    std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }
}

int validate_application_config() {
  path::CheckFileExistenceFATAL(cfg::get().config_file_path);

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

void create_logger_from_config() {
  using namespace logging;
  //Need to create properties file for logger and use it in build.
  logger* lg = create_logger();
  lg->add_writer(std::make_shared<file_writer>(cfg::get().out_path_to_logger_file));
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

bool organize_output_directory() {
  return true;
}