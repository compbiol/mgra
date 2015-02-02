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

  /*Different parse function*/
  void parse(std::unordered_map<std::string, std::vector<std::string> > const & input); 
  void parse_genomes(std::vector<std::string> const & input);
  void parse_trees(std::vector<std::string> const & input);
  void parse_algorithm(std::vector<std::string> const & input);
  void parse_target(std::vector<std::string> const & input);
  void parse_completion(std::vector<std::string> const & input);

public:  
  /*
   * Different strategy for build. 
   */
  bool is_target_build;
  bool is_assembly_build; 

  /*
   * Switch on/off stage fo reconstruction trees
   */
  bool is_reconstructed_trees;
  
  /*
   * Number of stage and rounds 
   */
  size_t rounds; 
  size_t stages;

  /*
   * Switch on/off bruteforce stage for small components
   */
  bool is_bruteforce;
  size_t size_component_in_bruteforce;

  /*
   * Switch on/off run blossom V stage
   */
  bool is_blossom; 

  /*
   * Switch on/off linearization algorithm
   */
  bool is_linearization_algo;

  std::list<twobreak_t> completion;

  using phylogeny_tree_t = structure::BinaryTree<mcolor_t>;
  std::vector<phylogeny_tree_t> phylotrees;

  mcolor_t target_mcolor;    
  
  std::string colorscheme;

private: 
  int RGBcoeff; 
  std::vector<std::string> RGBcolors;

  std::vector<std::string> priority_name;  
  std::unordered_map<std::string, size_t> genome_number;
  std::map<mcolor_t, std::string> mcolor_name;
  
private: 
  void init_basic_rgb_colors(); 
};

#include "structures/config_struct_impl.hpp"

using cfg = config_common::config<main_config<structure::Mcolor> >;

#endif
