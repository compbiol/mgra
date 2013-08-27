#ifndef MBGRAPH_HISTORY_H_
#define MBGRAPH_HISTORY_H_

#include "mbgraph_colors.h"

#include "genome_match.h" //FIXME REMOVE LATER

#include "2break.h"
#include "Insdel.h"
#include "TandemDuplication.h"

template<class mcolor_t>
struct mbgraph_with_history : public mbgraph_with_colors<mcolor_t> { 
  typedef event::TwoBreak<mcolor_t> twobreak_t; 
  typedef event::InsDel<mcolor_t> insertion_t;
  typedef event::TandemDuplication<mcolor_t> tandem_duplication_t;
  typedef std::list<twobreak_t> transform_t;
  
  template <class conf_t>
  mbgraph_with_history(const std::vector<MBGraph::genome_t>& genomes, const conf_t& cfg)
  : mbgraph_with_colors<mcolor_t>(genomes, cfg) 
  { 
  } 

  std::vector<transform_t> get_vec_TC_events() const;

  //2-break operations
  void apply_two_break(const twobreak_t& break2, bool record = true);

  inline typename transform_t::const_reverse_iterator crbegin_2break_history() const { 
    return break2_history.crbegin();
  } 
	
  inline typename transform_t::const_reverse_iterator crend_2break_history() const { 
    return break2_history.crend(); 
  }

  //Insertion/Deletion operations
  void apply_ins_del(const insertion_t& insdel);
 
  //(Reverse) tandem duplication operations
  void apply_tandem_duplication(const tandem_duplication_t& dupl, bool record = true);

  inline typename std::list<tandem_duplication_t>::const_iterator begin_tandem_duplication_history() const { 
    return tandem_dupl_history.cbegin();
  } 
	
  inline typename std::list<tandem_duplication_t>::const_iterator end_tandem_duplication_history() const { 
    return tandem_dupl_history.cend(); 
  }
	 
private: 
  std::list<twobreak_t> break2_history;
  std::list<tandem_duplication_t> tandem_dupl_history;		
};

template<class mcolor_t>
std::vector<std::list<event::TwoBreak<mcolor_t> > > mbgraph_with_history<mcolor_t>::get_vec_TC_events() const {
    std::vector<transform_t> transformations(this->vec_T_consistent_colors.size()); 
    for(auto it = crbegin_2break_history(); it != crend_2break_history(); ++it) {
	const mcolor_t& Q = it->get_mcolor();
	size_t i = 0; 
	for(const auto &color : this->vec_T_consistent_colors) {
	  if (Q == color) {
	    transformations[i].push_front(*it);
          }
	  ++i;
	} 
    } 
    return transformations;
}
 
template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_two_break(const twobreak_t& break2, bool record) { 
  if (record) {
    break2_history.push_back(break2);
  }
 
  for (const auto &color : break2) {
    for (size_t i = 0; i < 2; ++i) {
      if (break2.get_arc(i).first != Infty || break2.get_arc(i).second != Infty) {
	erase_edge(color.first, break2.get_arc(i).first, break2.get_arc(i).second);
      } 
    }

    if (break2.get_arc(0).first != Infty || break2.get_arc(1).first != Infty) {
      add_edge(color.first, break2.get_arc(0).first, break2.get_arc(1).first);
    }

    if (break2.get_arc(0).second != Infty || break2.get_arc(1).second != Infty) {
      add_edge(color.first, break2.get_arc(0).second, break2.get_arc(1).second);
    }
  }
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_ins_del(const insertion_t& insdel) {
  for (const auto &color : insdel) {
    if (insdel.is_deletion_oper()) {
      erase_edge(color.first, insdel.get_edge().first, insdel.get_edge().second);
    } else {
      add_edge(color.first, insdel.get_edge().first, insdel.get_edge().second);
    } 
  }
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_tandem_duplication(const tandem_duplication_t& dupl, bool record) {
  if (record) {
    tandem_dupl_history.push_back(dupl);
  }
 
  if (dupl.is_deletion_oper()) {
    for (auto it = dupl.cbegin_edges(); it != (--dupl.cend_edges()); ++it) {
      for (const auto &color : dupl) {
	erase_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (const auto &color : dupl) {
	add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (const auto &color : dupl) {
	erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } else {
    for (auto it = dupl.cbegin_edges(); it != (--dupl.cend_edges()); ++it) {
      for (const auto &color : dupl) {
	add_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (const auto &color : dupl) {
	erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (const auto &color : dupl) {
	add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } 
} 

#endif
