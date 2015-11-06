#ifndef MGRA_GONFIG_STRUCT_HPP
#define MGRA_GONFIG_STRUCT_HPP

#include "defined.hpp"
#include "JsonCpp/json/json.h"
#include "writer/txt_newick_tree.hpp"
#include "event/TwoBreak.hpp"
#include "event/WGD.hpp"

/**
 * Enum describe type of algorithm transforming multibreakpoint graph
 */
enum build_type {
    default_algo,
    target_algo,
    wgd_algo,
    user_algo
};

/**
 * Enum describe kind of stage which involved in algorithm
 */
enum kind_stage {
    balance_k,
    simple_path_k,
    four_cycles_k,
    fair_edge_k,
    clone_k,
    fair_clone_edge_k,
    components_k,
    change_canform_k,
    bruteforce_k,
    blossomv_k,
    completion_k
};

// This structure contains information from *.cfg file
template<class mcolor_t>
struct main_config {
    using wgd_t = event::WGD<mcolor_t>;
    using twobreak_t = event::TwoBreak<mcolor_t>;
    using phylogeny_tree_t = structure::BinaryTree<mcolor_t>;

    inline std::string const &get_RGBcolor(size_t index) const {
        return RGBcolors[RGBcoeff * index];
    }

    mcolor_t name_to_mcolor(std::string const &temp) const;

    std::string mcolor_to_name(mcolor_t const &temp) const;

    inline bool is_genome_name(std::string const &name) const {
        return (genome_number.count(name) != 0);
    }

    inline size_t get_genome_number(std::string const &name) const {
        assert(genome_number.count(name) != 0);
        return genome_number.find(name)->second;
    }

    inline std::string const &get_priority_name(size_t index) const {
        assert(index < priority_name.size());
        return priority_name[index];
    }

    DECLARE_DELEGATE_CONST_METHOD(size_t, priority_name, get_count_genomes, size)

    /**
     * Main function for init full config from JSON format
     */
    void load(Json::Value const &root);

    /**
     * Main function for save full config in JSON format
     */
    Json::Value save() const;

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
    std::string out_path_to_input_dir; // Path to input dir containing config and block
    std::string out_path_to_debug_dir; // Path to debug dir
    std::string out_path_to_saves_dir; // Path to saves dir
    std::string out_path_to_genomes_dir; //Path where resulting genomes are put
    std::string out_path_to_transfomations_dir; //Path where resulting genomes are put

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

    /**
     * Input WGD events
     */
    std::vector<wgd_t> wgds_events;

    /**
     * Input reconstruction target genome
     */
    mcolor_t target;

    /**
     * Strategy for mgra algorithm: default, target, users
     */
    build_type how_build;

    /**
     * Number of rounds (i.e. change number of split T-consistent colors)
     */
    size_t rounds;

    /**
     * Pipeline from stages which transforms multibreakpoint graph.
     */
    std::vector<kind_stage> pipeline;

    /**
     * Size of components which would be processed in bruteforce stage (time confusing).
     */
    size_t size_component_in_bruteforce;

    /**
     *
     */
    std::list<twobreak_t> completion;

    /**
     *
     */
    bool is_debug;

    /**
     *
     */
    bool is_saves;

    std::string colorscheme;

private:
    size_t RGBcoeff;
    std::vector<std::string> RGBcolors;

    std::map<mcolor_t, std::string> mcolor_name;

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

    void load_algorithm(Json::Value const &algo);

    void load_pipeline(Json::Value const &stages);

    void load_wgd_events(Json::Value const &wgds);

    void load_wgd_event(Json::Value const &wgd);

    void load_target(Json::Value const &target);

    void load_complections(Json::Value const &twobreaks);

    void load_complection(Json::Value const &twobreak);

    void load_output_directory(Json::Value const &path_to_dir);

    void load_saves(Json::Value const &enable_saves);

    void load_debug(Json::Value const &enable_debug);

    /**
     * Different save function
     */
    Json::Value save_genomes() const;

    Json::Value save_genome(size_t ind) const;

    Json::Value save_block_type() const;

    Json::Value save_files() const;

    Json::Value save_trees() const;

    Json::Value save_tree(size_t ind) const;

    Json::Value save_algorithm() const;

    Json::Value save_pipeline() const;

    Json::Value save_wgd_events() const;

    Json::Value save_wgd_event(size_t ind) const;

    Json::Value save_target() const;

    Json::Value save_complections() const;

    Json::Value save_complection(typename std::list<twobreak_t>::const_iterator const & twobreak) const;

    Json::Value save_output_directory() const;

    Json::Value save_saves() const;

    Json::Value save_debug() const;

    /**
     * Default strategy for init config
     */
    void default_rgb_colors();

    void default_directory_organization();

    void default_algorithm();

    void default_target_algorithm();
};

#include "config/config_struct_impl.hpp"

using cfg = config_common::config<main_config<structure::Mcolor> >;

#endif //MGRA_CONFIG_STRUCT_HPP
