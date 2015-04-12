//
// Created by Nikita Kartashov on 07/04/2015.
//

#ifndef MGRA_RECOVER_TREE_TASK_HPP_
#define MGRA_RECOVER_TREE_TASK_HPP_

#include "../../../utils/logger/logger.hpp"

#include "bruteforce_recover_tree_algorithm.hpp"
#include "dynamic_recover_tree_algorithm.hpp"
#include "multiedges_count_statistics_producer.hpp"
#include "simple_path_statistics_producer.hpp"
#include "algo/main/Balance.hpp"
#include "structures/tree.hpp"
#include "writer/newick_tree_printer.hpp"

namespace algo {

  template <class graph_pack_t>
  void recover_tree_task(graph_pack_t& graph_pack) {
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using algo_t = typename algo::RecoverTreeAlgorithm<graph_pack_t>;
    using algo_ptr = typename algo_t::algo_ptr;
    using tree_t = structure::BinaryTree<mcolor_t>;

    algo::Balance<graph_pack_t> balance_stage;
    balance_stage.run(graph_pack);

    INFO("Starting tree recovery")
    graph_pack.update_graph_statistics();

    algo_ptr recover_tree_algorithm;

    const size_t MAX_GENOMES_FOR_BRUTEFORCE = 50;

    std::shared_ptr<algo::StatisticsProducer<graph_pack_t> > producer;

    if (cfg::get().recover_tree_statistic == simple_paths) {
      producer = std::make_shared<algo::SimplePathStatisticsProducer<graph_pack_t> >(graph_pack);
    } else {
      producer = std::make_shared<algo::MultiEdgesCountStatisticsProducer<graph_pack_t> >(graph_pack);
    }

    if (cfg::get().get_count_genomes() < MAX_GENOMES_FOR_BRUTEFORCE) {
      recover_tree_algorithm =
          std::make_shared<algo::BruteforceRecoverTreeAlgorithm<graph_pack_t>>(graph_pack, producer);
    } else {
      recover_tree_algorithm =
          std::make_shared<algo::DynamicRecoverTreeAlgorithm<graph_pack_t>>(graph_pack, producer);
    }

    auto result_trees = recover_tree_algorithm->recover_trees();

    INFO("Recovered trees")

    const size_t MAX_TREES_TO_DUMP = 3;
    result_trees.resize(std::min(MAX_TREES_TO_DUMP, result_trees.size()));
    writer::GraphDot<graph_pack_t> dot_writer;
    std::ofstream summary_file(cfg::get().tree_summary_path);
    writer::NewickTreePrinter<tree_t> newick_tree_printer(summary_file);

    newick_tree_printer.print_trees(result_trees);

    for (size_t i = 0; i != result_trees.size(); ++i) {
      std::ofstream numbered_outfile(path::append_path(cfg::get().trees_path,
          std::to_string(i) + ".dot"));
      dot_writer.save_subtrees(numbered_outfile, {*result_trees[i]});
    }
    INFO("Dumped trees")
  }
}

#endif //MGRA_RECOVER_TREE_TASK_HPP_
