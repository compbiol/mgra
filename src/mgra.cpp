/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
*/

#include "command_line_parsing.hpp"

#include "algo/Algorithms.hpp"
#include "algo/recover_tree/bruteforce_recover_tree_algorithm.hpp"
#include "algo/recover_tree/dynamic_recover_tree_algorithm.hpp"
#include "algo/recover_tree/simple_path_statistics_producer.hpp"
#include "algo/recover_tree/multiedges_count_statistics_producer.hpp"

#include "logger/logger.hpp"
#include "logger/log_writers.hpp"

#include "RecoveredInfo.hpp"

#include "writer/Wgenome.h"
#include "writer/Wtransform.hpp"
#include "writer/newick_tree_printer.hpp"


int main(int argc, char** argv) {
  if (parse_config_from_command_line(argc, argv)) {
    std::cerr << "ERROR: error while parsing command line arguments" << std::endl;
  }

  if (validate_application_config()) {
    std::cerr << "ERROR: error while validating config" << std::endl;
  }

  if (!organize_output_directory()) {
    std::cerr << "ERROR: problem with organize output directory "
        << cfg::get().out_path_directory << std::endl;
    return 1;
  }

  create_logger_from_config();

//  Reading flags

//    Reading problem configuration and genomes
  using genome_t = structure::Genome;
  using mcolor_t = structure::Mcolor;
  using graph_pack_t = GraphPack<mcolor_t>;
  using algo_t = typename algo::RecoverTreeAlgorithm<graph_pack_t>;
  using algo_ptr = typename algo_t::algo_ptr;
  using tree_t = structure::BinaryTree<mcolor_t>;

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

//    Do job
  INFO("Start build graph")
  graph_pack_t graph_pack(genomes);
  INFO("End build graph")


  if (cfg::get().is_recover_tree) {
    algo::Balance<graph_pack_t> balance_stage;
    balance_stage.run(graph_pack);

    INFO("Starting tree recovery")
    //Recover tree here
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
    size_t dumped_so_far = 0;
    writer::GraphDot<graph_pack_t> dot_writer;
    writer::NewickTreePrinter<tree_t> newick_tree_printer(std::clog);

    for (auto const& tree: result_trees) {
      if (dumped_so_far == MAX_TREES_TO_DUMP || dumped_so_far == result_trees.size()) {
        break;
      }

      std::ofstream numbered_outfile(path::append_path(cfg::get().trees_path,
          std::to_string(dumped_so_far) + ".dot"));
      dot_writer.save_subtrees(numbered_outfile, {*tree});
      if (cfg::get().is_debug) {
        newick_tree_printer.print_tree(tree);
      }
      ++dumped_so_far;
    }

    INFO("Dumped trees")
  } else {
    {
      std::ostringstream out;
      out << "Determine " << graph_pack.multicolors.count_vec_T_consitent_color() << " \\vec{T}-consistent colors in tree:\n";
      for (auto id = graph_pack.multicolors.cbegin_vec_T_consistent_color(); id != graph_pack.multicolors.cend_vec_T_consistent_color(); ++id) {
        out << cfg::get().mcolor_to_name(*id) << " ";
      }
      INFO(out.str())
    }

    bool result = false; //wgd_algorithm(graph_pack);
    if (cfg::get().how_build == default_algo) {
      result = main_algorithm(graph_pack);
    } else if (cfg::get().how_build == target_algo) {
    }

    if (!result) {
      return 1;
    }

    INFO("Start linearization genomes.")
    RecoveredInfo<graph_pack_t> reductant(graph_pack);
    INFO("Finish linearization genomes.")

    INFO("Save transformations in files.")
    if (cfg::get().how_build == default_algo) {
      writer::Wtransformation<graph_pack_t> writer_transform(cfg::get().out_path_directory, graph_pack);
      auto recover_transformations = reductant.get_history();
      for (auto const& transformation : recover_transformations) {
        writer_transform.save_transformation(transformation.first, transformation.second);
        writer_transform.save_reverse_transformation(transformation.first, transformation.second);
      }
    }

    INFO("Save ancestor genomes in files.")
    writer::Wgenome<genome_t> writer_genome(cfg::get().genomes_path);
    writer_genome.save_genomes(reductant.get_genomes(), (cfg::get().how_build == default_algo));
  }
  INFO("MGRA log can be found here " << cfg::get().logger_path)
  INFO("Thank you for using MGRA!")

  return 0;
}
