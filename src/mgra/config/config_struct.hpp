#ifndef GONFIG_STRUCT_HPP
#define GONFIG_STRUCT_HPP

#include "defined.h"
#include "json/json.h"
#include "event/TwoBreak.hpp"

enum recover_tree_statistic_t {
  distribution,
  simple_paths
};

enum block_file_type_t {
  infercars,
  grimm
};

// This structure contains information from *.cfg file
template <class mcolor_t>
struct main_config {
  using twobreak_t = event::TwoBreak<mcolor_t>;

  mcolor_t name_to_mcolor(std::string const& temp) const;

  std::string mcolor_to_name(mcolor_t const& temp) const;

  inline bool is_genome_name(std::string const& name) const {
    return (genome_number.count(name) != 0);
  }

  inline size_t get_genome_number(std::string const& name) const {
    assert (genome_number.count(name) != 0);
    return genome_number.find(name)->second;
  }

  inline std::string const& get_priority_name(size_t index) const {
    assert(index < priority_name.size());
    return priority_name[index];
  }

  inline std::string const& get_RGBcolor(size_t index) const {
    return RGBcolors[RGBcoeff * index];
  }

  //FIXME: NEED TO REMOVE OR OPTIMIZE
  mcolor_t complete_color() const;

  DECLARE_DELEGATE_CONST_METHOD(size_t, priority_name, get_count_genomes, size)

  void load(Json::Value const & root);
  
  /** 
   * Fuction which can init all config
   */
  void parse(std::unordered_map<std::string, std::vector<std::string> > const& input);
  //void load(std::string const & path_to_file);
  //void save();

  bool is_debug;
  bool is_saves;
  block_file_type_t block_file_type;
  
  /**
  * Paths to different files
  */
  std::string config_file_path; //Path to input config file
  std::vector<std::string> path_to_blocks_file;
  std::string blocks_file_path; //Path to input genomes file with synteny blocks

  std::string out_path_directory; // Path to output directory
  std::string out_path_to_logger_file; // Path to logfile
  std::string out_path_to_input_dir; // Path to input dir containing config and block
  std::string out_path_to_debug_dir; // Path to debug dir
  std::string out_path_to_saves_dir; // Path to saves dir
  std::string out_path_to_genomes_dir; //Path where resulting genomes are put
  std::string out_path_to_transfomations_dir; //Path where resulting genomes are put
  std::string trees_path; //Path where recovered trees are put
  std::string tree_summary_path; //Path where summary of recovered trees in newick format is put

  /**
   * Different strategy for build. 
   */
  build_type how_build;

  /**
   * Switch on/off tree recovery mode
   */
  bool is_recover_tree;

  /**
  * The statistics provider used by recover tree algorithm
  */
  recover_tree_statistic_t recover_tree_statistic;
  using phylogeny_tree_t = structure::BinaryTree<mcolor_t>;
  std::vector<phylogeny_tree_t> phylotrees;

  /**
   * Number of stage and rounds 
   */
  size_t rounds;
  std::vector<algo::kind_stage> pipeline;
  bool is_linearization_ancestors; //Switch on/off stage for linearization procedure

  /**
   * Switch on/off bruteforce stage for small components
   */
  size_t size_component_in_bruteforce;

  std::list<twobreak_t> completion;

  mcolor_t target_mcolor;

  std::string colorscheme;

private:
  size_t RGBcoeff;
  std::vector<std::string> RGBcolors;

  std::vector<std::string> priority_name;
  std::unordered_map<std::string, size_t> genome_number;
  std::unordered_map<size_t, std::string> number_to_genome;
  
  std::map<mcolor_t, std::string> mcolor_name;

private:
  void default_rgb_colors();

  void default_algorithm();

  void default_target_algorithm();

  /**
   * Different load function
   */
  void load_genomes(Json::Value const& genomes);
  void load_genome(Json::Value const& genome, size_t index);
  void load_files(Json::Value const& path_to_files);
  void load_trees(Json::Value const& trees);
  void load_tree(Json::Value const& tree);
  void load_wgd_events(Json::Value const& wgds);
  void load_wgd_event(Json::Value const& wgd);
  void load_target(Json::Value const& target);
  void load_complections(Json::Value const& twobreaks); 
  void load_complection(Json::Value const& twobreak); 

  void load_output_directory(Json::Value const& path_to_dir); 
  void load_saves(Json::Value const& enable_saves); 
  void load_debug(Json::Value const& enable_debug); 
  
  /**
   * Different parse function
   */
  void parse_genomes(std::vector<std::string> const& input);

  void parse_trees(std::vector<std::string> const& input);

  void parse_algorithm(std::vector<std::string> const& input);

  void parse_target(std::vector<std::string> const& input);

  void parse_completion(std::vector<std::string> const& input);
};

#include "config/config_struct_impl.hpp"

using cfg = config_common::config<main_config<structure::Mcolor> >;

#endif
