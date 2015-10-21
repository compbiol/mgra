//
// Created by Nikita Kartashov on 06/04/2015.
//

#ifndef _MGRA_ARG_PARSING_HPP_
#define _MGRA_ARG_PARSING_HPP_

#include "defined.h" 

#include "logger/logger.hpp"
#include "logger/log_writers.hpp"
#include "io/copy_file.hpp"
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

#endif //_MGRA_ARG_PARSING_HPP_
