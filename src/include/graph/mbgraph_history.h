#ifndef MBGRAPH_HISTORY_H_
#define MBGRAPH_HISTORY_H_

#include "mbgraph_colors.h"

#include "genome_match.h"

#include "2break.h"
#include "Insdel.h"
#include "TandemDuplication.h"
#include "Fake2break.h"

//#include "genome_match.h"

template<class mcolor_t>
struct mbgraph_with_history : public mbgraph_with_colors<mcolor_t> { 
  typedef event::TwoBreak<mcolor_t> twobreak_t; 
  typedef event::FakeTwoBreak<mcolor_t> fake_twobreak_t;            
  typedef event::InsDel<mcolor_t> insertion_t;
  typedef event::TandemDuplication<mcolor_t> tandem_duplication_t;
  typedef std::list<twobreak_t> transform_t;
  typedef typename mbgraph_with_colors<mcolor_t>::mularcs_t mularcs_t; 
  
  
  template <class conf_t>
  mbgraph_with_history(std::vector<MBGraph::genome_t> const & genomes, conf_t const & cfg)
  : mbgraph_with_colors<mcolor_t>(genomes, cfg) 
  , number_events(0)
  { 
  } 
 
  void change_history(); 

  inline typename transform_t::const_iterator cbegin_2break_history() const { 
    return break2_history.cbegin();
  } 
	
  inline typename transform_t::const_iterator cend_2break_history() const { 
    return break2_history.cend(); 
  }

  inline typename transform_t::const_reverse_iterator crbegin_2break_history() const { 
    return break2_history.crbegin();
  } 
	
  inline typename transform_t::const_reverse_iterator crend_2break_history() const { 
    return break2_history.crend(); 
  }

  //2-break operations
  void apply_two_break(twobreak_t const & break2, bool record = true);

  //Insertion/Deletion operations
  void apply_ins_del(insertion_t const & insdel) { 
    for (auto const &color : insdel) {
      this->add_edge(color.first, insdel.get_edge().first, insdel.get_edge().second); 
    }
  } 

  //Fork operations
  void apply_fake_twobreak(fake_twobreak_t const & fakebreak2, bool record = true);
 
  //(Reverse) tandem duplication operations
  void apply_tandem_duplication(tandem_duplication_t const & dupl, bool record = true);

  inline typename std::list<tandem_duplication_t>::const_iterator begin_tandem_duplication_history() const { 
    return tandem_dupl_history.cbegin();
  } 
	
  inline typename std::list<tandem_duplication_t>::const_iterator end_tandem_duplication_history() const { 
    return tandem_dupl_history.cend(); 
  }

private:	
  inline void registrate_real_edge(std::pair<vertex_t, mcolor_t> const & mother_edge) { 
    target_real_edges[mother_edge.first].insert(mother_edge.second);
  } 

private:
  size_t number_events;
  std::list<twobreak_t> break2_history;
  std::map<size_t, twobreak_t> br2_history;
  std::map<size_t, fake_twobreak_t> fbr2_history;
  std::list<tandem_duplication_t> tandem_dupl_history;		

  std::map<vertex_t, std::set<mcolor_t> > target_real_edges;
};

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_two_break(twobreak_t const & break2, bool record) { 
  if (record) {
    br2_history.insert(std::make_pair(number_events++, break2));
  }
 
  for (auto const & color : break2) {
    for (size_t i = 0; i < 2; ++i) {
      if (break2.get_arc(i).first != Infty || break2.get_arc(i).second != Infty) {
	this->erase_edge(color.first, break2.get_arc(i).first, break2.get_arc(i).second);
      } 
    }

    if (break2.get_arc(0).first != Infty || break2.get_arc(1).first != Infty) {
      this->add_edge(color.first, break2.get_arc(0).first, break2.get_arc(1).first);
    }

    if (break2.get_arc(0).second != Infty || break2.get_arc(1).second != Infty) {
      this->add_edge(color.first, break2.get_arc(0).second, break2.get_arc(1).second);
    }
  }

  if (!target_real_edges.empty()) { 
    auto const & checkLambda = [&] (vertex_t const & v) -> void { 
      if (target_real_edges.count(v) != 0) { 
        mularcs_t const & mularcs = this->get_adjacent_multiedges(v); 
        auto target_colors = target_real_edges.find(v)->second; 
        for (auto arc = mularcs.cbegin(); arc != mularcs.cend(); ++arc) {
          auto const & color = arc->second;
          for (auto const & target : target_colors) { 
            mcolor_t inter_color(color, target, mcolor_t::Intersection);
            if (inter_color.size() == target.size()) {
              target_real_edges[v].erase(target);
              if (target_real_edges.find(v)->second.empty()) {
	        target_real_edges.erase(v);
              }   
            }
          } 
        } 
      } 
    };

    checkLambda(break2.get_arc(0).first); 
    checkLambda(break2.get_arc(0).second); 
    checkLambda(break2.get_arc(1).first); 
    checkLambda(break2.get_arc(1).second);
  } 
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply_fake_twobreak(fake_twobreak_t const & fakebreak2, bool record) { 
  if (record) { 
    fbr2_history.insert(std::make_pair(number_events++, fakebreak2));
  } 

  auto const & central_edge = fakebreak2.get_central_arc(); 
  auto const & mother_edge = fakebreak2.get_mother_edge();

  for (auto color = mother_edge.second.cbegin(); color != mother_edge.second.cend(); ++color) {
    this->erase_edge(color->first, central_edge.second, mother_edge.first);  
    this->add_edge(color->first, central_edge.first, central_edge.second);
  } 

  mularcs_t fathers = fakebreak2.get_end_edges();
  for (auto const & arc: fathers) { 
    for (auto color = arc.second.cbegin(); color != arc.second.cend(); ++color) {
      this->erase_edge(color->first, arc.first, central_edge.first);  
      this->add_edge(color->first, arc.first, mother_edge.first);
    } 
  } 

  registrate_real_edge(mother_edge); 
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::change_history() {
  for (int index = number_events - 1; index >= 0; --index) { 
    if (br2_history.find(index) != br2_history.end()) {
      break2_history.push_front(br2_history.find(index)->second); 
    } else if (fbr2_history.find(index) != fbr2_history.end()) {
      auto const & fakebreak2 = fbr2_history.find(index)->second;
      auto const & central = fakebreak2.get_central_arc();
      auto const & mother_edge = fakebreak2.get_mother_edge();
      auto const & colors = fakebreak2.get_end_edges().get_multicolors(); 
      
      //std::cerr << "Clone " << central.first << " " << central.second << " " << mother_edge.first << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;  

      auto last_twobreak = break2_history.end(); 
      twobreak_t old_two_break;
      for (auto br = break2_history.begin(); br != break2_history.end(); ++br) {
 
        auto const change_lambda = [&] (size_t ind) -> void { 
          if (br->get_vertex(ind) == mother_edge.first) {
          //std::cerr << "Change " << br->get_arc(0).first << " " << br->get_arc(0).second << " " 
  	//<< br->get_arc(1).first << " " << br->get_arc(1).second << " " << genome_match::mcolor_to_name(br->get_mcolor()) << std::endl;
            
	    auto const & color = br->get_mcolor();
            if (!color.includes(mother_edge.second) || color == mother_edge.second) { 
//!mcolor_t(color, mother_edge.second, mcolor_t::Intersection).empty()) { // 
              bool found = false;   
              for (auto col = colors.cbegin(); col != colors.cend() && !found; ++col) {
                mcolor_t inter_color(color, *col, mcolor_t::Intersection);
                if (inter_color.size() == col->size()) {
                  old_two_break = *br;
                  //std::cerr << "Start to change " << br->get_vertex(ind) << " " << central.first << std::endl;
                  br->change_vertex(ind, central.first);
                  last_twobreak = br;
                  found = true;
                }
              }
            }   
          }
        };
  
        change_lambda(0);
        change_lambda(1);
        change_lambda(2);
        change_lambda(3);

        //std::cerr << "result " << br->get_arc(0).first << " " << br->get_arc(0).second << " " 
  	//<< br->get_arc(1).first << " " << br->get_arc(1).second << " " << genome_match::mcolor_to_name(br->get_mcolor()) << std::endl;

      } 

      if (last_twobreak != break2_history.end()) {
        ++last_twobreak;
        if (old_two_break.get_vertex(0) == mother_edge.first) {
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(2), central.second, mother_edge.first, mother_edge.second));          
        } else if (old_two_break.get_vertex(1) == mother_edge.first) { 
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(3), central.second, mother_edge.first, mother_edge.second));
        } else if (old_two_break.get_vertex(2) == mother_edge.first) { 
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(0), central.second, mother_edge.first, mother_edge.second));
        } else if (old_two_break.get_vertex(3) == mother_edge.first) { 
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(1), central.second, mother_edge.first, mother_edge.second));
        } 
      } else { 
        ;//assert(false);
        //std::cerr << "Problem with " << central.first << " " << central.second << " " << mother_edge.first << std::endl;
      } 
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
	this->erase_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (const auto &color : dupl) {
	this->add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (const auto &color : dupl) {
	this->erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } else {
    for (auto it = dupl.cbegin_edges(); it != (--dupl.cend_edges()); ++it) {
      for (const auto &color : dupl) {
	this->add_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (const auto &color : dupl) {
	this->erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (const auto &color : dupl) {
	this->add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } 
} 

#endif
