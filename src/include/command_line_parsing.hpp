//
// Created by Nikita Kartashov on 06/04/2015.
//

#ifndef MGRA_ARG_PARSING_HPP_
#define MGRA_ARG_PARSING_HPP_

#include "version.hpp"
#include "defined.hpp"
#include "reader/reader.hpp"

#include "io/copy_file.hpp"
#include "logger/logger.hpp"
#include "logger/log_writers.hpp"
#include "config/config_struct.hpp"

/**
 * 
 */
int validate_application_config();

/**
 *
 */
int parse_config_from_command_line(int argc, char** argv);

/**
 * Organize output directory (creat folders, copy files) and change cfg.
 */
bool organize_output_directory();

/**
 * Create simply loger to output file. Need to parse our configuration file and trace different levels.
 */
void create_logger_from_config();

/**
 * Create directory with cheks
 */
inline bool create_dir_if_not_exists(std::string const &directory) {
    if (path::check_existence(directory)) {
        if (path::FileExists(directory)) return false;
    } else {
        return path::make_dir(directory);
    }
    return true;
}

/**
 * Parse input genomes from files
 */
inline std::vector<structure::Genome> parse_genomes() {
    std::vector<structure::Genome> genomes;

    INFO("Parse genomes file")
    if (cfg::get().block_file_type == infercars) {
        for (auto file = cfg::get().path_to_blocks_file.cbegin();
             file != cfg::get().path_to_blocks_file.cend(); ++file) {
            genomes = reader::read_infercars(*file);
        }
    } else if (cfg::get().block_file_type == grimm) {
        for (auto file = cfg::get().path_to_blocks_file.cbegin();
             file != cfg::get().path_to_blocks_file.cend(); ++file) {
            genomes = reader::read_grimm(*file);
        }
    }
    INFO("End parsing genomes file")

    for (size_t i = 0; i < genomes.size(); ++i) {
        std::ostringstream out;
        out << "Download genome " << genomes[i].get_name() << " with " << genomes[i].size() << " blocks.";
        INFO(out.str())
    }

    return genomes;
}

#endif //_MGRA_ARG_PARSING_HPP_
