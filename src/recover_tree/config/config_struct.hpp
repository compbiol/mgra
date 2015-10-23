//
// Created by pavel on 10/21/15.
//

#ifndef MGRA_CONFIG_STRUCT_HPP
#define MGRA_CONFIG_STRUCT_HPP

#include "json/json.h"

template<class mcolor_t>
struct tree_config {
    using twobreak_t = event::TwoBreak<mcolor_t>;

    void load(Json::Value const &root);

private:
    void load_genomes(Json::Value const &genomes);

    void load_genome(Json::Value const &genome, size_t index);

    void load_files(Json::Value const &path_to_files);

    void load_trees(Json::Value const &trees);

    void load_tree(Json::Value const &tree);

public:
    std::string config_file_path; //Path to input config file

    std::vector<std::string> paths_to_blocks_file;

    std::string out_path_directory; // Path to output directory

    std::string out_path_to_logger_file; // Path to logfile

    std::string trees_path; //Path where recovered trees are put

    std::string tree_summary_path; //Path where summary of recovered trees in newick format is put

    recover_tree_statistic_t recover_tree_statistic;

    using phylogeny_tree_t = structure::BinaryTree<mcolor_t>;
    std::vector<phylogeny_tree_t> phylotrees;

private:
    std::vector<std::string> priority_name;

    std::unordered_map<std::string, size_t> genome_number;

    std::unordered_map<size_t, std::string> number_to_genome;
};

#include "config/config_struct_impl.hpp"

using cfg = config_common::config<tree_config<structure::Mcolor> >;

#endif //MGRA_CONFIG_STRUCT_HPP
