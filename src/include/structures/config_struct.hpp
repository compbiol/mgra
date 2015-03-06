#ifndef GONFIG_STRUCT_HPP
#define GONFIG_STRUCT_HPP

#include "defined.h"

#include "event/TwoBreak.hpp"

/*This structures containes information from *.cfg file */
template<class mcolor_t>
struct main_config {
  using twobreak_t = event::TwoBreak<mcolor_t>; 

  mcolor_t name_to_mcolor(std::string const & temp) const;
  std::string mcolor_to_name(mcolor_t const & temp) const;

  inline bool is_genome_name(std::string const & name) const {
    return (genome_number.count(name) != 0);
  } 

  inline size_t get_genome_number(std::string const & name) const {
    assert (genome_number.count(name) != 0); 
    return genome_number.find(name)->second;
  } 

  inline std::string const & get_priority_name(size_t index) const {
    assert(index < priority_name.size());
    return priority_name[index];
  }

  inline std::string const & get_RGBcolor(size_t index) const { 
    return RGBcolors[RGBcoeff * index];
  } 

  DECLARE_DELEGATE_CONST_METHOD(size_t, priority_name, get_count_genomes, size)

  /*fuction which can init all config*/
  void parse(std::unordered_map<std::string, std::vector<std::string> > const & input);  
  //void save();

public:  
  bool is_debug; 
  std::string out_path_to_debug_dir; 

  /**
   * Different strategy for build. 
   */
  build_type how_build;  

  /**
   * Switch on/off tree recovery mode
   */
  bool is_recover_tree;

  using phylogeny_tree_t = structure::BinaryTree<mcolor_t>;
  std::vector<phylogeny_tree_t> phylotrees;
  
  /**
   * Number of stage and rounds 
   */
  size_t rounds; 
  std::vector<algo::kind_stage> pipeline;

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
  std::map<mcolor_t, std::string> mcolor_name;
  
private: 
  void default_rgb_colors(); 
  void default_algorithm(); 
  void default_target_algorithm();
  
  /**
   * Different parse function
   */
  void parse_genomes(std::vector<std::string> const & input);
  void parse_trees(std::vector<std::string> const & input);
  void parse_algorithm(std::vector<std::string> const & input);
  void parse_target(std::vector<std::string> const & input);
  void parse_completion(std::vector<std::string> const & input);
};

#include "structures/config_struct_impl.hpp"

using cfg = config_common::config<main_config<structure::Mcolor> >;

#endif
