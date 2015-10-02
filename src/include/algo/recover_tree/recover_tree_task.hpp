//
// Created by Nikita Kartashov on 07/04/2015.
//

#ifndef MGRA_RECOVER_TREE_TASK_HPP_
#define MGRA_RECOVER_TREE_TASK_HPP_

#include <sstream>

#include "../../../utils/logger/logger.hpp"

#include "bruteforce_recover_tree_algorithm.hpp"
#include "dynamic_recover_tree_algorithm.hpp"
#include "statistics/multiedges_count_statistics_producer.hpp"
#include "statistics/simple_path_statistics_producer.hpp"
#include "algo/main/Balance.hpp"
#include "structures/tree.hpp"
#include "writer/txt_newick_tree.hpp"

namespace algo {

  template <class graph_pack_t>
  void recover_tree_task(graph_pack_t& graph_pack) {
    using namespace algo;
    using namespace structure;
    using namespace writer;
    using mcolor_t = typename graph_pack_t::mcolor_type;
    using algo_t = RecoverTreeAlgorithm<graph_pack_t>;
    using algo_ptr = typename algo_t::algo_ptr;
    using tree_t = BinaryTree<mcolor_t>;

    const size_t ROUNDS_FOR_TREE_RECOVERY = 1;
    //FIXME: change recovery into stage and add
    StageManager<graph_pack_t> algorithm(ROUNDS_FOR_TREE_RECOVERY,
                                         {cfg::get().is_debug, cfg::get().out_path_to_debug_dir});

    Balance<graph_pack_t>* balance_stage = new Balance<graph_pack_t>();
    algorithm.add_stage(balance_stage);
    algorithm.run(graph_pack);

    INFO("Starting tree recovery")
    graph_pack.update_graph_statistics();

    algo_ptr recover_tree_algorithm;

    const size_t MAX_GENOMES_FOR_BRUTEFORCE = 50;

    std::shared_ptr<algo::StatisticsProducer<graph_pack_t> > producer;

    if (cfg::get().recover_tree_statistic == simple_paths) {
      producer = std::make_shared<SimplePathStatisticsProducer<graph_pack_t> >(graph_pack);
    } else {
      producer = std::make_shared<MultiEdgesCountStatisticsProducer<graph_pack_t> >(graph_pack);
    }

    if (cfg::get().get_count_genomes() < MAX_GENOMES_FOR_BRUTEFORCE) {
      recover_tree_algorithm =
          std::make_shared<BruteforceRecoverTreeAlgorithm<graph_pack_t>>(graph_pack, producer);
    } else {
      recover_tree_algorithm =
          std::make_shared<DynamicRecoverTreeAlgorithm<graph_pack_t>>(graph_pack, producer);
    }

    auto result_trees = recover_tree_algorithm->get_result();

    INFO("Recovered trees")

    const size_t MAX_TREES_TO_DUMP = 3;
    result_trees.resize(std::min(MAX_TREES_TO_DUMP, result_trees.size()));
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
    }
    INFO("Dumped trees")
  }
}

#endif //MGRA_RECOVER_TREE_TASK_HPP_
