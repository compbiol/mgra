//
// Created by Nikita Kartashov on 06/04/2015.
//

#ifndef _MGRA_ARG_PARSING_HPP_
#define _MGRA_ARG_PARSING_HPP_

int parse_config_from_command_line(int argc, char** argv);
int validate_application_config();
bool organize_output_directory();
void create_logger_from_config();

#endif //_MGRA_ARG_PARSING_HPP_
