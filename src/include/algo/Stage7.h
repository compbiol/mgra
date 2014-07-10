#ifndef STAGE7_H_
#define STAGE7_H_

template<class graph_t>
bool Algorithm<graph_t>::stage7() {
  size_t number_rear = 0; // number of rearrangements 

  for (auto const & x: *graph) {  
    mularcs_t const & mularcs = graph->get_adjacent_multiedges(x);
      
    bool found = false;
    for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
      vertex_t const & y = im->first; // Q == im->second - color of central edge
      if (y == Infty) {
        continue;
      }

      if (!is_mobility_edge(x, y)) {  
        mularcs_t&& mularcs_y = graph->get_adjacent_multiedges_with_info(y);
        mularcs_y.erase(x);
        mularcs_t&& mularcs_x = graph->get_adjacent_multiedges_with_info(x);
        mularcs_x.erase(y);
       
        vertex_t const & mother = mularcs_y.cbegin()->first;
        bool sligshot = (mularcs_y.size() == 1) && (mularcs_x.size() != 1) 
			&& graph->is_vec_T_consistent_color(mularcs_y.cbegin()->second) && (mother != Infty);
        for (auto arc = mularcs_x.cbegin(); arc != mularcs_x.cend() && sligshot; ++arc) {
          sligshot = graph->is_vec_T_consistent_color(arc->second);
        } 
        sligshot = sligshot && mularcs_y.cbegin()->second == mularcs_x.union_multicolors();
  
       
        if (sligshot) {
          size_t count = 0;

          /*mularcs_t&& end_s = graph->get_adjacent_multiedges_with_info(mother);
          end_s.erase(y);
          for (auto arc = end_s.cbegin(); arc != end_s.cend(); ++arc) { 
            if (!mularcs_x.defined(arc->first) && canformQ(arc->first, mularcs_y.cbegin()->second)) {
              ++count;
            }
          } */ 
          
          if (count == 0) { 
            fake_twobreak_t ft(arc_t(x, y), mularcs_x, *(mularcs_y.cbegin()));  
            //std::cerr << "Create clone " << x << " " << y << " " /*<< *(mularcs_y.cbegin())*/ << std::endl;
            graph->apply_fake_twobreak(ft);
            assert(graph->get_adjacent_multiedges(x).number_unique_edge() == 1 && graph->get_adjacent_multiedges(x).union_multicolors() == graph->get_complete_color()); 
            ++number_rear;
          } 
        } 
      }
    } 
  }  

  return (number_rear != 0);
}

template<class graph_t>
bool Algorithm<graph_t>::stage71() {
  size_t number_rear = 0; // number of rearrangements 

  for (auto const & x: *graph) {  
    mularcs_t const & mularcs = graph->get_adjacent_multiedges(x);
      
    bool found = false;
    for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
      vertex_t const & y = im->first; // Q == im->second - color of central edge
    
      if (y == Infty) {
        continue;
      }

      if (!graph->is_T_consistent_color(im->second)) { 
      mularcs_t mularcs_x = graph->get_adjacent_multiedges(x);
      mularcs_x.erase(y);

      mularcs_t mularcs_y = graph->get_adjacent_multiedges(y);
      mularcs_y.erase(x);
 
      size_t count = 0;
      bool is_y = false;
      std::pair<vertex_t, mcolor_t> pr;
      for (auto left = mularcs_x.begin(); left != mularcs_x.end(); ++left) {  
        if (!graph->is_T_consistent_color(left->second)) {
          pr = *left;
          ++count;
        } 
      }   	

      for (auto right = mularcs_y.begin(); right != mularcs_y.end(); ++right) { 
        if (!graph->is_T_consistent_color(right->second)) {
          pr = *right;
          is_y = true;
          ++count;
        } 
      }

      if (count == 1 && pr.first != Infty) { 
        size_t count_f = graph->split_color(im->second).size(); 
        size_t count_s = graph->split_color(pr.second).size();

        if (count_f == count_s) { 
	  //std::cerr << "Situation " << x << " " << y << " " << pr.first << std::endl;
          mularcs_t mul_x = graph->get_adjacent_multiedges_with_info(x, false);
          mularcs_t mul_y = graph->get_adjacent_multiedges_with_info(y, false);
          size_t count_central = calculate_cost(y, mul_x, mul_y);
          size_t count_another = 0;
          if (is_y) { 
            auto temp = graph->get_adjacent_multiedges_with_info(pr.first, false);
            temp.erase(y);
            count_another = calculate_cost(pr.first, graph->get_adjacent_multiedges_with_info(y, false), temp);
	  } else {
            auto temp = graph->get_adjacent_multiedges_with_info(pr.first, false);
            temp.erase(y);
	    count_another = calculate_cost(pr.first, graph->get_adjacent_multiedges_with_info(x, false), temp);
          } 

          if (count_central == std::min(count_central, count_another)) {
            mul_x.erase(y);
            std::vector<twobreak_t> history;  
            bool good = true;
            for (auto const &arc : mul_x) {
	      vertex_t const & v = mularcs_y.get_vertex(arc.second); 
              if (!v.empty() && graph->is_vec_T_consistent_color(arc.second)) {
                  //std::cerr << "Two_break " << x << " " << arc.first << " " << y << " " << v << " " << genome_match::mcolor_to_name(arc.second) << std::endl;                
                  history.push_back(twobreak_t(x, arc.first, y, v, arc.second));
        	  //graph->apply_two_break(twobreak_t(x, arc.first, y, v, arc.second));
              } else { 
                good = false;
              }  
            }

            if (good) { 
              found = true;
              for (auto const & break2 : history) {
                graph->apply_two_break(break2);
                ++number_rear; 
              }
            } 
          } 
        }
      }
      }  
    } 
  }
 
  return (number_rear != 0);
}
 
/*template<class graph_t>
bool Algorithm<graph_t>::stage7() {
  size_t number_rear = 0; // number of rearrangements 

  for (auto const & x: *graph) {  
    mularcs_t const & mularcs = graph->get_adjacent_multiedges(x);
      
    bool found = false;
    for(auto im = mularcs.cbegin(); (im != mularcs.cend()) && !found; ++im) {
      vertex_t const & y = im->first; // Q == im->second - color of central edge
      if (y == Infty) {
        continue;
      }
 
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
         vertex_t const & mother = mularcs_y.cbegin()->first; 
         mularcs_t&& end_s = graph->get_adjacent_multiedges_with_info(mother);
         end_s.erase(y);
         for (auto arc = end_s.cbegin(); arc != end_s.cend(); ++arc) { 
           if (!mularcs_x.defined(arc->first) && canformQ(arc->first, mularcs_y.cbegin()->second)) {
             ++count;
           }
         }  
           
         if (count == 0) { 
           fake_twobreak_t ft(arc_t(x, y), mularcs_x, *(mularcs_y.cbegin()));  
           //std::cerr << mother << " " << y << std::endl;
           graph->apply_fake_twobreak(ft);
           graph->registrate_real_edge(mother, mularcs_y.cbegin()->second, arc_t(x, y));
           assert(graph->get_adjacent_multiedges(x).number_unique_edge() == 1 && graph->get_adjacent_multiedges(x).union_multicolors() == graph->get_complete_color()); 
           ++number_rear;
         } 
      } 
    }
  }  

  return (number_rear != 0);
}*/ 

#endif
