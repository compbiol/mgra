#ifndef MBGRAPH_HISTORY_H_
#define MBGRAPH_HISTORY_H_

#include "genome_match.h"

#include "graph/mbgraph_colors.h"

#include "event/2break.h"
#include "event/Insdel.h"
#include "event/TandemDuplication.h"
#include "event/Fake2break.h"

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
  void new_change_history(size_t index, vertex_t const & where);

  /*bool check_edge_with_pseudo_vertex() const {
    for (auto const & pseudo_vertex : pseudo_mother_vertex) {
      mcolor_t color = this->get_adjacent_multiedges(pseudo_vertex).get_multicolor(Infty);
      if (color != this->get_complete_color()) {
        return false;
      }  
    } 
    return true;
  }*/
  
  void apply(twobreak_t const & break2, bool record = true); //2-break operations
  void apply(fake_twobreak_t const & fakebreak2, bool record = true); //Fork operations 
  void apply(tandem_duplication_t const & dupl, bool record = true); //(Reverse) tandem duplication operations
  void apply(insertion_t const & insdel) { //Insertion/Deletion operations
    for (auto const &color : insdel) {
      this->add_edge(color.first, insdel.get_edge().first, insdel.get_edge().second); 
    }
  }   

  typedef typename transform_t::const_iterator twobreak_citer;
  typedef typename transform_t::const_reverse_iterator twobreak_criter;
  DECLARE_CONST_ITERATOR( twobreak_citer, break2_history, cbegin_2break_history, cbegin ) 
  DECLARE_CONST_ITERATOR( twobreak_citer, break2_history, cend_2break_history, cend ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, break2_history, crbegin_2break_history, crbegin ) 
  DECLARE_CONST_ITERATOR( twobreak_criter, break2_history, crend_2break_history, crend ) 

private:	
  inline void registrate_real_edge(std::pair<vertex_t, mcolor_t> const & mother_edge, size_t index) { 
    target_real_edges[mother_edge.first].insert(std::make_pair(mother_edge.second, index));
  } 

private:
  enum type_action {twobreak_action, clone_action, tandem_duplication_action};
  std::list<std::pair<type_action, size_t> > complete_history;

  size_t number_events;
  std::list<twobreak_t> break2_history;
  std::map<size_t, twobreak_t> br2_history;
  std::map<size_t, fake_twobreak_t> fbr2_history;
  std::vector<tandem_duplication_t> tandem_duplication_history;		

  //std::unordered_set<vertex_t> pseudo_mother_vertex;
  std::map<vertex_t, std::set<std::pair<mcolor_t, size_t> > > target_real_edges;
};

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply(twobreak_t const & break2, bool record) { 
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
      if (target_real_edges.count(v) != 0) { // && pseudo_mother_vertex.count(v) != 0) { 
        mularcs_t const & mularcs = this->get_adjacent_multiedges(v); 
        auto target_colors = target_real_edges.find(v)->second; 
        for (auto arc = mularcs.cbegin(); arc != mularcs.cend(); ++arc) {
          auto const & color = arc->second;
          for (auto const & target : target_colors) { 
            mcolor_t inter_color(color, target.first, mcolor_t::Intersection);
            if (inter_color.size() == target.first.size()) {
              //std::cerr << "V: " << v << " " << arc->first << " " << target_colors.size() << std::endl;
              //std::cerr << genome_match::mcolor_to_name(color) << " " << genome_match::mcolor_to_name(target.first) << std::endl;
              this->new_change_history(target.second, arc->first);
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
void mbgraph_with_history<mcolor_t>::apply(fake_twobreak_t const & fakebreak2, bool record) { 
  if (record) { 
    fbr2_history.insert(std::make_pair(number_events++, fakebreak2));
  } 

  auto const & central_edge = fakebreak2.get_central_arc(); 
  auto const & mother_edge = fakebreak2.get_mother_edge();

  if (fakebreak2.is_have_pseudo_vertex()) {
    this->vertex_set.insert(mother_edge.first);
    for (auto color = mother_edge.second.cbegin(); color != mother_edge.second.cend(); ++color) {
      this->erase_edge(color->first, central_edge.second, Infty);  
      this->add_edge(color->first, central_edge.second, mother_edge.first);
    }  

    auto compl_color = this->get_complement_color(mother_edge.second);
    for (auto color = compl_color.cbegin(); color != compl_color.cend(); ++color) {
      this->add_edge(color->first, mother_edge.first, Infty);
    }
  }

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
 
  if (fakebreak2.is_have_pseudo_vertex()) {
    registrate_real_edge(mother_edge, number_events - 1);   
    //pseudo_mother_vertex.insert(mother_edge.first);
  }

} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::new_change_history(size_t index, vertex_t const & where) {
  if (fbr2_history.find(index) != fbr2_history.end()) {
    auto const & fakebreak2 = fbr2_history.find(index)->second;
    auto const & mother_edge = fakebreak2.get_mother_edge();
    
    //std::cerr << "start to check" << std::endl;  //pseudo_mother_vertex.count(mother_edge.first) != 0
    if (fakebreak2.is_have_pseudo_vertex()) { //pseudo_mother_vertex.count(mother_edge.first) != 0) { 
      //std::cerr << mother_edge.first << std::endl;
      //std::cerr << "Where = " << where << std::endl;

      //std::cerr << genome_match::mcolor_to_name(mother_edge.second) << std::endl;
      for (auto color = mother_edge.second.cbegin(); color != mother_edge.second.cend(); ++color) {
        this->erase_edge(color->first, where, mother_edge.first); 
        if (where != Infty) {  
          this->add_edge(color->first, where, Infty);  
        }
      } 

      mcolor_t other_color = this->get_adjacent_multiedges(mother_edge.first).get_multicolor(Infty); 
      //std::cerr << genome_match::mcolor_to_name(other_color) << std::endl;
      for (auto color = other_color.cbegin(); color != other_color.cend(); ++color) {
        this->erase_edge(color->first, mother_edge.first, Infty);  
      }
      this->vertex_set.erase(mother_edge.first);
      

      /*if (fakebreak2.is_have_pseudo_vertex()) { 
        mcolor_t other_color = this->get_adjacent_multiedges(mother_edge.first).get_multicolor(Infty); 
        std::cerr << genome_match::mcolor_to_name(other_color) << std::endl;
        for (auto color = other_color.cbegin(); color != other_color.cend(); ++color) {
          this->erase_edge(color->first, mother_edge.first, Infty);  
        }
        this->vertex_set.erase(mother_edge.first);
        std::cerr << "End to clean. remove " << mother_edge.first << std::endl;
      }*/
    } else { 
      ;//std::cerr << mother_edge.first << std::endl;
      //assert(false);
    }    
  } else {
    assert(false);
  }
} 

template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::change_history() {
  //std::cerr << number_events << std::endl;
  for (int index = number_events - 1; index >= 0; --index) { 
    if (br2_history.find(index) != br2_history.end()) {
      break2_history.push_front(br2_history.find(index)->second); 
      //std::cerr << "2-break " << index << std::endl;
    } else if (fbr2_history.find(index) != fbr2_history.end()) {
      auto const & fakebreak2 = fbr2_history.find(index)->second;
      auto const & central = fakebreak2.get_central_arc();
      auto const & mother_edge = fakebreak2.get_mother_edge();
      auto const & colors = fakebreak2.get_end_edges().get_multicolors(); 
      

      //std::cerr << "Clone " << index << " " << central.first << " " << central.second << " " << mother_edge.first << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;  

      auto last_twobreak = break2_history.end(); 
      twobreak_t old_two_break;
      for (auto br = break2_history.begin(); br != break2_history.end(); ++br) {
        
        auto const change_lambda = [&] (size_t ind) -> void { 
          if (br->get_vertex(ind) == mother_edge.first) {
          //std::cerr << "Change " << br->get_arc(0).first << " " << br->get_arc(0).second << " " 
          //<< br->get_arc(1).first << " " << br->get_arc(1).second << " " << genome_match::mcolor_to_name(br->get_mcolor()) << std::endl;
            
            auto const & color = br->get_mcolor();
            if (!color.includes(mother_edge.second)) { // || color == mother_edge.second) { 
              //std::cerr << genome_match::mcolor_to_name(color) << " " << genome_match::mcolor_to_name(mother_edge.second) 
              //  << " " << color.includes(mother_edge.second) << std::endl;
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

        //std::cerr << "Start to change " << br->get_arc(0).first << " " << br->get_arc(0).second << " " 
         //<< br->get_arc(1).first << " " << br->get_arc(1).second << " " << genome_match::mcolor_to_name(br->get_mcolor()) << std::endl; 
      } 

      if (last_twobreak != break2_history.end()) {
        ++last_twobreak;
        
        vertex_t mother = mother_edge.first; 
        if (fakebreak2.is_have_pseudo_vertex()) { //&& pseudo_mother_vertex.count(mother_edge.first) != 0) { 
          //std::cerr << "Change on infinity vertex" << std::endl;
          mother = Infty;
        } 

        //std::cerr << "Insert two break " << std::endl;

        if (old_two_break.get_vertex(0) == mother_edge.first) {
          //std::cerr << "1: " << central.first << " " << old_two_break.get_vertex(2) << " " << central.second << " " 
          //<< mother << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(2), central.second, mother, mother_edge.second));          
        } else if (old_two_break.get_vertex(1) == mother_edge.first) { 
          //std::cerr << "2: " << central.first << " " << old_two_break.get_vertex(3) << " " << central.second << " " 
          //<< mother << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(3), central.second, mother, mother_edge.second));
        } else if (old_two_break.get_vertex(2) == mother_edge.first) { 
          //std::cerr << "3: " << central.first << " " << old_two_break.get_vertex(0) << " " << central.second << " " 
          //<< mother << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(0), central.second, mother, mother_edge.second));
        } else if (old_two_break.get_vertex(3) == mother_edge.first) { 
          //std::cerr << "4: " << central.first << " " << old_two_break.get_vertex(1) << " " << central.second << " " 
          //<< mother << " " << genome_match::mcolor_to_name(mother_edge.second) << std::endl;
          break2_history.insert(last_twobreak, twobreak_t(central.first, old_two_break.get_vertex(1), central.second, mother, mother_edge.second));
        } 
      } else { 
        ;//assert(false);
        //std::cerr << "Problem with " << central.first << " " << central.second << " " << mother_edge.first << std::endl;
      } 
    }
  } 
} 


template<class mcolor_t>
void mbgraph_with_history<mcolor_t>::apply(tandem_duplication_t const & dupl, bool record) {
  if (record) {
    tandem_duplication_history.push_back(dupl);
  }
 
  if (dupl.is_deletion_oper()) {
    for (auto it = dupl.cbegin_edges(); it != (--dupl.cend_edges()); ++it) {
      for (auto const &color : dupl) {  
        this->erase_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (auto const &color : dupl) {
        this->add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (auto const &color : dupl) {
        this->erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } else {
    for (auto it = dupl.cbegin_edges(); it != (--dupl.cend_edges()); ++it) {
      for (auto const &color : dupl) {
        this->add_edge(color.first, it->first, it->second);
      }
    } 

    if (dupl.is_reverse_tandem_duplication()) {
      for (auto const &color : dupl) {
        this->erase_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    } else {
      for (auto const &color : dupl) {
      	this->add_edge(color.first, (--dupl.cend_edges())->first, (--dupl.cend_edges())->second);
      }
    }
  } 
} 

#endif
