//
// Created by Nikita Kartashov on 06/04/2015.
//

#include "command_line_parsing.hpp"
#include "tclap/CmdLine.h"


int parse_config_from_command_line(int argc, char **argv) {
    try {
        std::string const VERSION(std::to_string(MGRA_VERSION_MAJOR) + "." + std::to_string(MGRA_VERSION_MINOR) + "." +
                                  std::to_string(MGRA_VERSION_PATCH));

        TCLAP::CmdLine cmd(
                "MGRA (Multiple Genome Rearrangements & Ancestors) (c) 2008-2015 by Pavel Avdeyev, Nikita Kartashov, Shuai Jiang, Max Alekseyev. Distributed under GNU GENERAL PUBLIC LICENSE license.",
                ' ', VERSION);

        TCLAP::ValueArg<std::string> path_to_cfg_file_arg("c",
                                                          "config",
                                                          "Input configure file",
                                                          true,
                                                          "",
                                                          "filename",
                                                          cmd);

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

        cmd.parse(argc, const_cast<const char *const *>(argv));

        cfg::get_writable().config_file_path = path_to_cfg_file_arg.getValue();
        cfg::get_writable().out_path_directory = path::make_full_path(output_arg.getValue());
        cfg::get_writable().how_build = config::default_algo;
        cfg::get_writable().is_debug = debug_arg.getValue();
        cfg::get_writable().is_saves = saves_arg.getValue();

        std::ifstream cfg_file(path_to_cfg_file_arg.getValue());
        cfg::create_instance(cfg_file);

        return 0;
    } catch (TCLAP::ArgException &e) {
        std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }
}

int validate_application_config() {
    for (auto file = cfg::get().path_to_blocks_file.cbegin(); file != cfg::get().path_to_blocks_file.end(); ++file) {
        path::CheckFileExistenceFATAL(*file);
    }

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
    // Need to create properties file for logger and use it in build.
    logger *lg = create_logger();
    lg->add_writer(std::make_shared<file_writer>(cfg::get().out_path_to_logger_file));
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

bool organize_output_directory() {
    // Create different directories.
    bool result = create_dir_if_not_exists(cfg::get().out_path_to_genomes_dir) &&
                  create_dir_if_not_exists(cfg::get().out_path_to_transfomations_dir) &&
                  create_dir_if_not_exists(cfg::get().out_path_to_input_dir);

    if (cfg::get().is_debug) {
        result = result && create_dir_if_not_exists(cfg::get().out_path_to_debug_dir);
    }


    if (cfg::get().is_saves) {
        result = result && create_dir_if_not_exists(cfg::get().out_path_to_saves_dir);
    }


    //Copy input files
    for (auto file = cfg::get_writable().path_to_blocks_file.begin(); file != cfg::get_writable().path_to_blocks_file.end(); ++file) {
        path::copy_file(*file, path::append_path(cfg::get().out_path_to_input_dir, path::filename(*file)));
        *file = path::append_path(cfg::get().out_path_to_input_dir, path::filename(*file));
    }

    //Save config
    cfg::get_writable().config_file_path = path::append_path(cfg::get().out_path_to_input_dir, "config.txt");
    std::ofstream save_cfg(cfg::get_writable().config_file_path);
    Json::StyledStreamWriter cfg_writer; cfg_writer.write(save_cfg, cfg::get().save());

    return result;
}