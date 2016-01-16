//
// Created by pavel on 10/21/15.
//

#ifndef MGRA_CONFIG_STRUCT_HPP
#define MGRA_CONFIG_STRUCT_HPP

#include "JsonCpp/json/json.h"
//#include "writer/txt_newick_tree.hpp"

enum recover_tree_statistic_t {
    patterns,
    distribution
};

template<class mcolor_t>
struct tree_config {
    using phylogeny_tree_t = structure::phyl_tree::BinaryTree<std::string>;

    /**
     * Main function for init full config from JSON format
     */
    void load(Json::Value const &root);

    /**
     * Main function for save full config in JSON format
     */
    Json::Value save() const;

    DECLARE_DELEGATE_CONST_METHOD(size_t, priority_name, get_count_genomes, size)
    /**
     * Path to config file
     */
    std::string config_file_path; //Path to input config file

    /**
     * Information about blocks file
     */
    block_file_type_t block_file_type; //Type of input dataset
    std::vector<std::string> path_to_blocks_file; //Path to input genomes file with synteny blocks

    /**
     * Paths to different out directories
     */
    std::string out_path_directory; // Path to output directory
    std::string out_path_to_logger_file; // Path to logfile
    std::string trees_path; //Path where recovered trees are put //TODO rename out_
    std::string tree_summary_path; //Path where summary of recovered trees in newick format is put

    /**
     * Information about genomes
     */
    std::vector<std::string> priority_name;
    std::unordered_map<std::string, size_t> genome_number;
    std::unordered_multimap<size_t, std::string> number_to_genome;

    /**
     * Input (sub) phylotrees
     */
    std::vector<phylogeny_tree_t> phylotrees;

    recover_tree_statistic_t recover_tree_statistic;

    std::string colorscheme; //TODO REMOVE
private:
    /**
     * Different load function
     */
    void load_genomes(Json::Value const &genomes);

    void load_genome(Json::Value const &genome, size_t index);

    void load_block_type(Json::Value const &type_files);

    void load_files(Json::Value const &path_to_files);

    void load_trees(Json::Value const &trees);

    void load_tree(Json::Value const &tree);

    void load_output_directory(Json::Value const &path_to_dir);

    void load_source_statistics(Json::Value const &source_stat);

    /**
     * Different save function
     */
    Json::Value save_genomes() const;

    Json::Value save_genome(size_t ind) const;

    Json::Value save_block_type() const;

    Json::Value save_files() const;

    Json::Value save_trees() const;

    Json::Value save_tree(size_t ind) const;

    Json::Value save_output_directory() const;

    Json::Value save_source_statistics() const;

    /**
     * Default strategy for init config
     */
    void default_rgb_colors();

    void default_directory_organization();

private:
    size_t RGBcoeff; //TODO REMOVE
    std::vector<std::string> RGBcolors; //TODO REMOVE
};

#include "config/config_struct_impl.hpp"

using cfg = config_common::config<tree_config<structure::Mcolor> >;

#endif //MGRA_CONFIG_STRUCT_HPP
