#ifndef STAGE7_H_
#define STAGE7_H_

template<class graph_t>
bool Algorithm<graph_t>::stage7() {
  bool isChanged = false; 
  size_t number_rear = 0; // number of rearrangements 

  for (const auto& x: *graph) {  
    const mularcs_t& mularcs = graph->get_adjacent_multiedges(x);
      
    bool found = false;
    for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
      const vertex_t& y = im->first; // Q == im->second - color of central edge
      if (y == Infty) {
        continue;
      }
 
      const auto canformQ_Lambda = [&] (const vertex_t& x, const mcolor_t& Q) -> size_t {
        if (x == Infty) {
          return canformQoo;
        }
        const mularcs_t& mularcs = graph->get_adjacent_multiedges_with_info(x); 
    	bool canform = true;
	for(auto arc = mularcs.cbegin(); (arc != mularcs.cend()) && canform; ++arc) { 
          mcolor_t color(Q, arc->second, mcolor_t::Intersection); 
          if (color.size() > 0 && color.size() < arc->second.size()) { 
            canform = false;
          }
        }
        return canform;
     };

      mularcs_t&& mularcs_x = graph->get_adjacent_multiedges_with_info(x);
      mularcs_x.erase(y);
      mularcs_t&& mularcs_y = graph->get_adjacent_multiedges_with_info(y);
      mularcs_y.erase(x);

      bool sligshot = true;
      for (auto arc = mularcs_x.cbegin(); arc != mularcs_x.cend() && sligshot; ++arc) {
        sligshot = graph->is_vec_T_consistent_color(arc->second);
      } 
       
      if (sligshot && (mularcs_y.size() == 1) && 
           graph->is_vec_T_consistent_color(mularcs_y.cbegin()->second) && (mularcs_y.cbegin()->first != Infty) &&  
           mularcs_y.cbegin()->second == mularcs_x.union_multicolors()
	   && graph->is_T_consistent_color(im->second)) {
         size_t count = 0;
         const vertex_t& mother = mularcs_y.cbegin()->first; 
         mularcs_t&& end_s = graph->get_adjacent_multiedges_with_info(mother);
         end_s.erase(y);
         /*for (auto arc = end_s.cbegin(); arc != end_s.cend(); ++arc) { 
           if (!mularcs_x.defined(arc->first) && canformQ_Lambda(arc->first, mularcs_y.cbegin()->second)) {
             ++count;
           }
         }  
           
         if (count == 0) { */
           fake_twobreak_t ft(arc_t(x, y), mularcs_x, *(mularcs_y.cbegin()));  
           //std::cerr << mother << " " << y << std::endl;
           graph->apply_fake_twobreak(ft);
           graph->registrate_real_edge(mother, mularcs_y.cbegin()->second, arc_t(x, y));
           assert(graph->get_adjacent_multiedges(x).number_unique_edge() == 1 && graph->get_adjacent_multiedges(x).union_multicolors() == graph->get_complete_color()); 
           ++number_rear;
         //} 
      } 
    }
  }  

  return (number_rear != 0);
} 

#endif
