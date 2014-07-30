#ifndef GRAPHS_HISTORY_HPP
#define GRAPHS_HISTORY_HPP

#include "event/2break.h"
#include "event/Insdel.h"
#include "event/TandemDuplication.h"
#include "event/Fake2break.h"

template<class mcolor_t>
struct HistoryGraph { 
  typedef event::TwoBreak<mcolor_t> twobreak_t; 
  typedef event::FakeTwoBreak<mcolor_t> clone_t;            
  typedef event::InsDel<mcolor_t> insertion_t;
  typedef event::TandemDuplication<mcolor_t> tandem_duplication_t;
  typedef std::list<twobreak_t> transform_t;

  inline void save_twobreak(twobreak_t const & twobreak) {
    twobreak_history.push_back(twobreak);
    complete_history.push_back(std::make_pair(twobreak_action, twobreak_history.size() - 1)); 
  }

  inline void save_insertion(insertion_t const & insertion) {
    insertion_history.push_back(insertion);
    complete_history.push_back(std::make_pair(insertion_action, insertion_history.size() - 1));
  }

  inline size_t save_clone(clone_t const & clone) {
    clone_history.push_back(clone);
    complete_history.push_back(std::make_pair(clone_action, clone_history.size() - 1));
    return (clone_history.size() - 1);
  }

  inline void save_tandem_duplication(tandem_duplication_t const & tandem_duplication) {
    tandem_duplication_history.push_back(tandem_duplication);
    complete_history.push_back(std::make_pair(tandem_duplication_action, tandem_duplication_history.size() - 1));
  }

  inline clone_t const & get_clone(size_t index) const { 
    return clone_history[index];
  } 

  void change_history(); 

  typedef typename transform_t::const_iterator twobreak_citer;
  typedef typename transform_t::const_reverse_iterator twobreak_criter;
  DECLARE_CONST_ITERATOR( twobreak_citer, break2_history, cbegin_2break_history, cbegin ) 
  DECLARE_CONST_ITERATOR( twobreak_citer, break2_history, cend_2break_history, cend ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, break2_history, crbegin_2break_history, crbegin ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, break2_history, crend_2break_history, crend ) 

private:
  enum type_action {insertion_action, twobreak_action, clone_action, tandem_duplication_action};

  std::list<std::pair<type_action, size_t> > complete_history;
  std::vector<insertion_t> insertion_history;
  std::vector<twobreak_t> twobreak_history;
  std::vector<clone_t> clone_history;
  std::vector<tandem_duplication_t> tandem_duplication_history;   

  std::list<twobreak_t> break2_history;
}; 

template<class mcolor_t>
void HistoryGraph<mcolor_t>::change_history() {  
  for (auto action = complete_history.rbegin(); action != complete_history.rend(); ++action) { 
    if (action->first == twobreak_action) {
      break2_history.push_front(twobreak_history[action->second]); 
    } else if (action->first == clone_action) { 
      auto const & fakebreak2 = clone_history[action->second];
      auto const & central = fakebreak2.get_central_arc();
      auto const & mother_edge = fakebreak2.get_mother_edge();
      auto const & end_edges = fakebreak2.get_end_edges(); 

      //std::cerr << "Clone " << central.first << " " << central.second << " " << mother_edge.first << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;  

      auto last_twobreak = break2_history.end(); 
      twobreak_t old_two_break;
      for (auto br = break2_history.begin(); br != break2_history.end(); ++br) {
        auto const change_lambda = [&] (size_t ind) -> void { 
          if (br->get_vertex(ind) == mother_edge.first) {
          //std::cerr << "Change " << br->get_vertex(0) << " " << br->get_vertex(1) << " " 
          //<< br->get_vertex(2) << " " << br->get_vertex(3) << " " << genome_match::mcolor_to_name(br->get_mcolor()) << std::endl;
            
            auto const & color = br->get_mcolor();
            if (!color.includes(mother_edge.second)) { 
              //std::cerr << genome_match::mcolor_to_name(color) << " " << genome_match::mcolor_to_name(mother_edge.second) 
              //  << " " << color.includes(mother_edge.second) << std::endl;
              //!mcolor_t(color, mother_edge.second, mcolor_t::Intersection).empty()) { // 
              bool found = false;   
              for (auto arc = end_edges.cbegin(); arc != end_edges.cend() && !found; ++arc) {
                mcolor_t inter_color(color, arc->second, mcolor_t::Intersection);
                if (inter_color.size() == arc->second.size()) {
                  old_two_break = *br;
                  //std::cerr << "Start to change " << br->get_vertex(ind) << " " << central.first << std::endl;
                  br->change_vertex(ind, central.first);
                  last_twobreak = br;
                  found = true;
                }
              }
            }

            if (fakebreak2.is_have_pseudo_vertex() && mother_edge.first == br->get_vertex(ind) && color.includes(mother_edge.second)) {
              //std::cerr << "Start to change on infinity " << br->get_vertex(ind) << " " << Infty << std::endl;
              br->change_vertex(ind, Infty);
            }
          }
        };

        change_lambda(0);
        change_lambda(1);
        change_lambda(2);
        change_lambda(3);
      } 

      if (last_twobreak != break2_history.end()) {
        ++last_twobreak;
        
        vertex_t mother = mother_edge.first; 
        if (fakebreak2.is_have_pseudo_vertex()) { 
          //std::cerr << "Change on infinity vertex" << std::endl;
          mother = Infty;
        } 

        vertex_t where; 
        if (old_two_break.get_vertex(0) == mother_edge.first) { 
          where = old_two_break.get_vertex(2);
        } else if (old_two_break.get_vertex(1) == mother_edge.first) { 
          where = old_two_break.get_vertex(3);
        } else if (old_two_break.get_vertex(2) == mother_edge.first) { 
          where = old_two_break.get_vertex(0);
        } else if (old_two_break.get_vertex(3) == mother_edge.first) { 
          where = old_two_break.get_vertex(1);
        } else { 
          assert(false);
        } 

        twobreak_t two_break = twobreak_t(central.first, where, central.second, mother, mother_edge.second);
        break2_history.insert(last_twobreak, two_break);
        //std::cerr << "Insert two break " << std::endl;
        //std::cerr << central.first << " " << where << " " << central.second << " " 
        //<< mother << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;
      } 
    }
  }
} 

#endif