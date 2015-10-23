//
// Created by pavel on 10/21/15.
//

#include "defined.hpp"
#include "config/config_struct.hpp"
#include "command_line_parsing.hpp"

int main(int argc, char** argv) {
    if (parse_config_from_command_line(argc, argv)) {
        std::cerr << "ERROR: error while parsing command line arguments" << std::endl;
        return 1;
    }

    /*if (validate_application_config()) {
        std::cerr << "ERROR: error while validating config" << std::endl;
        return 1;
    }

    if (!organize_output_directory()) {
        std::cerr << "ERROR: problem with organize output directory " << cfg::get().out_path_directory << std::endl;
        return 1;
    }

    create_logger_from_config();

    //Reading problem configuration and genomes
    using genome_t = structure::Genome;
    using mcolor_t = structure::Mcolor;
    using graph_pack_t = GraphPack<mcolor_t>;

    INFO("Parse genomes file")
    std::vector<genome_t> genomes;
    if (cfg::get().block_file_type == infercars) {
        genomes = reader::read_infercars(cfg::get().blocks_file_path);
    } else if (cfg::get().block_file_type == grimm) {
        genomes = reader::read_grimm(cfg::get().blocks_file_path);
    }

    for (size_t i = 0; i < genomes.size(); ++i) {
        std::ostringstream out;
        out << "Download genome " << cfg::get().get_priority_name(i) << " with " << genomes[i].size() << " blocks.";
        INFO(out.str())
    }

    //Do job
    INFO("Start build graph")
    graph_pack_t graph_pack(genomes);
    INFO("End build graph")

    INFO("Start recover trees from breakpoint graph")
    algo::recover_tree_task(graph_pack);
    INFO("Start recover trees from breakpoint graph")
    */

    return 0;
}