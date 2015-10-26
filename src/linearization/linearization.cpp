//
// Created by pavel on 10/26/15.
//

#include "defined.hpp"
#include "config/config_struct.hpp"
#include "command_line_parsing.hpp"

int main(int argc, char** argv) {
    if (parse_config_from_command_line(argc, argv)) {
        std::cerr << "ERROR: error while parsing command line arguments" << std::endl;
        return 1;
    }

    return 0;
}