//
// Created by pavel on 10/21/15.
//

#include "command_line_parsing.hpp"

#include "graph/graph_pack.hpp"

#include "statistics/statistics_producer.hpp"
#include "statistics/distribution_producer.hpp"
#include "statistics/pattern_producer.hpp"

#include "algo/balance.hpp"
#include "algo/bruteforce_algorithm.hpp"
//#include "algo/dynamic_algorithm.hpp"

const size_t MAX_TREES_TO_DUMP = 3;

namespace algo {

template<class graph_pack_t>
std::vector<std::shared_ptr<structure::phyl_tree::BinaryTree<typename graph_pack_t::mcolor_type>>> recover_tree_task(
        graph_pack_t &graph_pack) {
    using namespace algo;
    using namespace structure;

    using algo_t = RecoverTreeAlgorithm<graph_pack_t>;
    using algo_ptr = typename algo_t::algo_ptr;

    //using mcolor_t = typename graph_pack_t::mcolor_type;
    //using tree_t = BinaryTree<mcolor_t>;

    INFO("Calculate statistics")
    graph_pack.update_graph_statistics();
    std::shared_ptr<algo::StatisticsProducer<graph_pack_t>> producer;
    if (cfg::get().recover_tree_statistic == patterns) {
        producer = std::make_shared<PatternProducer<graph_pack_t>>(graph_pack);
    } else {
        producer = std::make_shared<DistributionProducer<graph_pack_t> >(graph_pack);
    }

    INFO("Starting tree recovery")
    size_t const MAX_GENOMES_FOR_BRUTEFORCE = 10;
    algo_ptr recover_tree_algorithm;
    if (cfg::get().get_count_genomes() < MAX_GENOMES_FOR_BRUTEFORCE) {
        recover_tree_algorithm = std::make_shared<BruteforceRecoverTreeAlgorithm<graph_pack_t>>(graph_pack,
                                                                                                producer);
    } /*else {
    recover_tree_algorithm =
            std::make_shared<DynamicRecoverTreeAlgorithm<graph_pack_t>>(graph_pack, producer);
    }*/

    auto result_trees = recover_tree_algorithm->get_result();
    result_trees.resize(std::min(MAX_TREES_TO_DUMP, result_trees.size()));
    return result_trees;
}

}

int main(int argc, char** argv) {
    using genome_t = structure::Genome;
    //using mcolor_t = structure::Mcolor;
    //using graph_pack_t = GraphPack<mcolor_t>;

    if (parse_config_from_command_line(argc, argv)) {
        std::cerr << "ERROR: error while parsing command line arguments" << std::endl;
        return 1;
    }

    if (validate_application_config()) {
        std::cerr << "ERROR: error while validating config" << std::endl;
        return 1;
    }

    if (!organize_output_directory()) {
        std::cerr << "ERROR: problem with organize output directory " << cfg::get().out_path_directory << std::endl;
        return 1;
    }

    create_logger_from_config();
    std::vector<genome_t> genomes = parse_genomes();

    /*Build graph*/
    /*
    INFO("Start build graph")
    graph_pack_t graph_pack(genomes, cfg::get().phylotree, cfg::get().genome_number, cfg::get().target);
    algo::balance(graph_pack);
    INFO("End build graph")
    */
    /*Do main job*/
    //INFO("Start recover trees from breakpoint graph")
    //algo::recover_tree_task(graph_pack);
    //INFO("Start recover trees from breakpoint graph")

    /*Save different output information in files*/
    INFO("Dumped trees")
    //const size_t MAX_TREES_TO_DUMP = 3;
    /*result_trees.resize(std::min(MAX_TREES_TO_DUMP, result_trees.size()));
    if (cfg::get().is_debug) {
        TXT_NewickTree<tree_t> newick_tree_printer(std::clog);
        newick_tree_printer.print_trees(result_trees);
    }
    GraphDot<graph_pack_t> dot_writer;
    std::ofstream summary_file(cfg::get().tree_summary_path);
    TXT_NewickTree<tree_t> newick_tree_printer(summary_file);

    newick_tree_printer.print_trees(result_trees);

    for (size_t i = 0; i != result_trees.size(); ++i) {
        std::ofstream numbered_outfile(path::append_path(cfg::get().trees_path,
                                                         std::to_string(i) + ".dot"));
        dot_writer.save_subtrees(numbered_outfile, {*result_trees[i]});
    }*/

    return 0;
}