#ifndef MGRA_GONFIG_STRUCT_HPP
#define MGRA_GONFIG_STRUCT_HPP

#include "defined.hpp"
#include "event/TwoBreak.hpp"
#include "event/wgd.hpp"

namespace config {

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
    irregular_fair_edge_k,
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
struct mgra_config : public abs_config {
    using wgd_t = event::wgd<std::string>;
    using twobreak_t = event::TwoBreak<mcolor_t>; //TODO think about twobreaks, no multicolor, just names, but vertex_ptr (no compare with names)

    inline std::string const &get_RGBcolor(size_t index) const { //FIXME delete output in dot file
        assert((RGBcoeff * index) < RGBcolors.size());
        return RGBcolors[RGBcoeff * index];
    }

    mcolor_t name_to_mcolor(std::string const &temp) const; //FIXME replace anywhere

    std::string mcolor_to_name(mcolor_t const &temp) const; //FIXME replace anywhere

    /**
     * Main function for init full config from JSON format
     */
    void load(Json::Value const &root) override;

    /**
     * Main function for save full config in JSON format
     */
    Json::Value save() const override;

    /**
     * Paths to different out directories
     */
    std::string out_path_to_input_dir; // Path to input dir containing config and block
    std::string out_path_to_debug_dir; // Path to debug dir
    std::string out_path_to_saves_dir; // Path to saves dir
    std::string out_path_to_genomes_dir; //Path where resulting genomes are put
    std::string out_path_to_transfomations_dir; //Path where resulting genomes are put

    /**
     * Input WGD events
     */
    std::map<std::pair<std::string, std::string>, wgd_t> wgds_events;

    /**
     * Input reconstruction target genome
     */
    std::string target;

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

    std::string colorscheme; //FIXME delete output in dot file

private:
    size_t RGBcoeff; //FIXME delete output in dot file
    std::vector<std::string> RGBcolors; //FIXME delete output in dot file

private:
    /**
     * Different load function
     */
    void load_algorithm(Json::Value const &algo);

    void load_pipeline(Json::Value const &stages);

    void load_wgd_events(Json::Value const &wgds);

    void load_wgd_event(Json::Value const &wgd);

    void load_target(Json::Value const &target);

    void load_complections(Json::Value const &twobreaks);

    void load_complection(Json::Value const &twobreak);

    void load_saves(Json::Value const &enable_saves);

    void load_debug(Json::Value const &enable_debug);

    /**
     * Different save function
     */
    Json::Value save_algorithm() const;

    Json::Value save_pipeline() const;

    Json::Value save_wgd_events() const;

    Json::Value save_wgd_event(size_t ind) const;

    Json::Value save_target() const;

    Json::Value save_complections() const;

    Json::Value save_complection(typename std::list<twobreak_t>::const_iterator const & twobreak) const;

    Json::Value save_saves() const;

    Json::Value save_debug() const;

    /**
     * Default strategy for init config
     */
    void default_rgb_colors();

    void default_directory_organization();

    void default_algorithm();

    void default_wgd_algorithm();

    void default_target_algorithm();

};

}

#include "config_struct_impl.hpp"

using cfg = config_common::config<config::mgra_config<structure::Mcolor>>;

#endif //MGRA_CONFIG_STRUCT_HPP
