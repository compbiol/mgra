//
// Created by pavel on 1/8/16.
//

#ifndef MGRA_ABS_CONFIG_HPP
#define MGRA_ABS_CONFIG_HPP

#include "JsonCpp/json/json.h"
#include "structures/phyl_tree/tree.hpp"

enum block_file_type_t { infercars, grimm };

namespace config {

struct abs_config {
    using phylogeny_tree_t = structure::phyl_tree::BinaryTree<std::string>;

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

    /**
     * Information about genomes
     */
    std::vector<std::string> priority_name;
    std::unordered_map<std::string, size_t> genome_number;
    std::unordered_multimap<size_t, std::string> number_to_genome;

    /**
     * Input phylotree
     */
    phylogeny_tree_t phylotree;

    DECLARE_DELEGATE_CONST_METHOD(size_t, priority_name, get_count_genomes, size)

    /**
     * Main function for init full config from JSON format
     */
    virtual void load(Json::Value const &root) = 0;

    /**
     * Main function for save full config in JSON format
     */
    virtual Json::Value save() const = 0;

    virtual ~abs_config() {
    }

protected:
    /**
     * Different load function
     */
    void load_genomes(Json::Value const &genomes);

    void load_genome(Json::Value const &genome, size_t index);

    void load_block_type(Json::Value const &type_files);

    void load_files(Json::Value const &path_to_files);

    void load_tree(Json::Value const &tree);

    void load_output_directory(Json::Value const &path_to_dir);

    /**
     * Different save function
     */
    Json::Value save_genomes() const;

    Json::Value save_genome(size_t ind) const;

    Json::Value save_block_type() const;

    Json::Value save_files() const;

    Json::Value save_tree() const;

    Json::Value save_output_directory() const;
};

}

#include "abs_config_impl.hpp"

#endif //MGRA_ABS_CONFIG_HPP
