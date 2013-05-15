#include "estimate.h"

Statistics::Statistics(const MBGraph& gr)
: graph(gr) {  
  count_compl_multiedges();
  count_not_compl_multiedges();	
  count_chromosomes();  
	
#ifdef VERSION2 
  //count_some_statistics();
  count_weak_simple_vertex();
#endif
}	      	

std::vector<std::string> Statistics::get_compl_stat() const { 
  std::multimap<size_t, std::string> answer;

  for(auto im = compl_multiedges_count.cbegin(); im != compl_multiedges_count.cend(); ++im) {
    const Mcolor& current = graph.get_complement_color(im->first);  // complementary multicolor.

    if (im->first < current) {
      continue;
    } 
	 
    size_t m1 = calc_value(compl_multiedges_count, current) / 2;
    size_t m2 = (im->second) / 2;	
    size_t paths = calc_value(simple_vertices_count, current) - (calc_value(simple_multiedges_count, current) + calc_value(simple_multiedges_count, im->first)) - calc_value(simple_vertices_alone_count, current) - calc_value(special_cycle_count, current);
    size_t cycles = calc_value(simple_cycle_count, current) + calc_value(special_cycle_count, current);
 	
    std::ostringstream os;
 
    os << "{";		    
    if (graph.is_T_consistent_color(im->first)) { 
      os << "\\bf ";
    } 

    const Mcolor& first = graph.get_min_complement_color(current); 
    const Mcolor& second = graph.get_complement_color(first);

    os <<  genome_match::mcolor_to_name(first) << " + "  <<  genome_match::mcolor_to_name(second) << "} & " 
      // multiedges
       << m1 << " + " << m2 << " = " << m1 + m2 << " & " 
      // simple vertices
       << calc_value(simple_vertices_count, current) << " & "  
      // simple multiedges
       << calc_value(simple_multiedges_count, current) << " + " << calc_value(simple_multiedges_count, im->first) << " = " << calc_value(simple_multiedges_count, current) + calc_value(simple_multiedges_count, im->first) << " & "
      // simple paths + cycles
       << paths << " + " << cycles << " = " << paths + cycles << " & "
      // irregular multiedges
       << calc_value(good_irrer_multiedges_count, current) << " + " << calc_value(good_irrer_multiedges_count, im->first) << " = " << calc_value(good_irrer_multiedges_count, current) + calc_value(good_irrer_multiedges_count, im->first);
	
    answer.insert(std::make_pair(m1 + m2, os.str()));
  }
	
  std::vector<std::string> output;
  for(auto it = answer.rbegin(); it != answer.rend(); ++it) {
    output.push_back(it->second); 
  }  	
  return output;
} 

std::vector<std::string> Statistics::get_no_compl_stat() const { 
  std::multimap<size_t, std::string> answer;

  for(auto im = not_compl_multiedges_count.cbegin(); im != not_compl_multiedges_count.cend(); ++im) {
    size_t vm1 = im->second;
		
    std::ostringstream os;
  
    os << "{";		    
    if (graph.is_T_consistent_color(im->first)) { 
      os << "\\bf ";
    } 
    
    os <<  genome_match::mcolor_to_name(im->first) << "} & " << vm1;

    answer.insert(std::make_pair(vm1, os.str()));
  } 

  std::vector<std::string> output;
  for(auto it = answer.rbegin(); it != answer.rend(); ++it) {
    output.push_back(it->second); 
  }  	
  return output;	
}

std::vector<Mcolor> Statistics::get_new_color() const { 
  std::vector<Mcolor> output(compl_multiedges_count.size());
  for(auto im = compl_multiedges_count.cbegin(); im != compl_multiedges_count.cend(); ++im) { 
    output.push_back(im->first);
  }  
  return output;	
} 

void Statistics::count_weak_simple_vertex() { 
  size_t count = 0; 

  for(auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
    mularcs_t Mx = graph.get_adjacent_multiedges(*it);

    if (graph.is_duplication_vertice(Mx)) { 
	++count;
    } 
  } 

  std::cerr << "Duplication vertex: " << count << std::endl;
} 

void Statistics::count_some_statistics() { 
  std::map<std::pair<Mcolor, Mcolor>, size_t> count_vertex;
  
  for(auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
    mularcs_t Mx = graph.get_adjacent_multiedges(*it);
    
    if (Mx.size() == 2) { 
	if (Mx.cbegin()->second < Mx.crbegin()->second) {
		++count_vertex[std::make_pair(Mx.cbegin()->second, Mx.crbegin()->second)];
	} else { 
		++count_vertex[std::make_pair(Mx.crbegin()->second, Mx.cbegin()->second)];	
	} 
    }

    Mx.erase(Infty);

    if (Mx.size() == 2) { 
	if (Mx.cbegin()->second < Mx.crbegin()->second) {
		++count_vertex[std::make_pair(Mx.cbegin()->second, Mx.crbegin()->second)];
	} else { 
		++count_vertex[std::make_pair(Mx.crbegin()->second, Mx.cbegin()->second)];	
	} 
    } 
  } 

  std::multimap<size_t, std::pair<Mcolor, Mcolor> > edges; 

  for(auto it = count_vertex.cbegin(); it != count_vertex.cend(); ++it) { 
	edges.insert(std::make_pair(it->second, it->first));
  } 

  std::ofstream ofstat("new_stat.txt");
  for(auto it = edges.crbegin(); it != edges.crend(); ++it) { 
	ofstat << "{" << genome_match::mcolor_to_name(it->second.first) << "+" << genome_match::mcolor_to_name(it->second.second) << "} " << it->first << std::endl; 
  } 
  ofstat.close();
} 

void Statistics::count_compl_multiedges() {
  std::unordered_set<std::string> processed;

  for(auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
    //multimularcs_t current = graph.get_adjacent_multiedges_v2(*it); //current is list with adjacent multiedges
    mularcs_t current = graph.get_adjacent_multiedges(*it); //current is list with adjacent multiedges

    ++multidegree_count[current.size()]; //current.size - is degree vertex *it

    if (graph.is_simple_vertice(current)) {  //we define simple vertices as a regular vertex of multidegree 2. 
      processed.insert(*it);
      ++simple_vertices_count[std::min(current.begin()->second, current.rbegin()->second)]; //simple vertices because degree 2.
    }

    for(auto im = current.cbegin(); im != current.cend(); ++im) {
      if (!im->second.is_good_multiedge()) {  /*|| !graph.is_T_consistent_color(im->second)) { */
	continue;
      } 

      ++compl_multiedges_count[im->second];   // count two times, because same underected edge (u, v) and (v, u)
			
      if (graph.is_simple_vertice(current)) {
	++good_multiedges_count[im->second]; //good if one vertices have degree 2
	
	if (im->first == Infty) { 
	  ++good_multiedges_count[im->second];
	} 
						
	if (processed.find(im->first) != processed.end()) {  
	  ++simple_multiedges_count[im->second]; //if two vertices have degree = 2 - is simple edges
	} 
      } 

      if (im->first == Infty) { 
	++compl_multiedges_count[im->second];
	++good_irrer_multiedges_count[im->second];			
      } 
    }
  }

  // count lonely vertices (short paths) 
  for(auto it = processed.cbegin(); it != processed.cend(); ++it) {
    //multimularcs_t current = graph.get_adjacent_multiedges_v2(*it);
    mularcs_t current = graph.get_adjacent_multiedges(*it);
    if (processed.find(current.begin()->first) == processed.end() && processed.find(current.rbegin()->first) == processed.end()) {
      ++simple_vertices_alone_count[std::min(current.begin()->second, current.rbegin()->second)]; //no good neighbors
    }
  } 
} 

void Statistics::count_not_compl_multiedges() { 
  for(auto it = graph.begin_vertices(); it != graph.end_vertices(); ++it) {
    //multimularcs_t current = graph.get_adjacent_multiedges_v2(*it);
    mularcs_t current = graph.get_adjacent_multiedges(*it);  
    for (auto jt = current.cbegin(); jt != current.cend(); ++jt) {
      if (jt->second.is_good_multiedge()) {
	continue; 
      } 

      ++not_compl_multiedges_count[jt->second]; 

      if (jt->first == Infty) { 
	++not_compl_multiedges_count[jt->second];		
      } 
    } 				 
  } 	
} 

std::map<std::pair<Mcolor, Mcolor>, size_t> Statistics::get_Hsubgraph() { 
  std::unordered_set<vertex_t> processed;
	
  for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {
    mularcs_t Mx = graph.get_adjacent_multiedges(*is);

    if (graph.is_fair_vertice(Mx)) { 
      for(auto im = Mx.cbegin(); im != Mx.cend(); ++im) {
	if (im->first == Infty || processed.find(im->first) != processed.end()) { 
	  continue; 
	} 

	mularcs_t My = graph.get_adjacent_multiedges(im->first);
			
	if (graph.is_fair_vertice(My)) {
	  mularcs_t Mx0 = Mx;
	  
	  Mx0.erase(im->first);
	  My.erase(*is);
	  
	  //Mcolor Q1 = Mx0.begin()->second;
	  //Mcolor Q2 = Mx0.rbegin()->second;
	  if ((Mx0.begin()->second == My.begin()->second) || (Mx0.rbegin()->second == My.begin()->second)) {
	    Mcolor QQ1 = graph.get_min_complement_color(Mx0.begin()->second);
	    Mcolor QQ2 = graph.get_min_complement_color(Mx0.rbegin()->second);	    
	    ++Hcount[std::make_pair(QQ1, QQ2)];
	    ++Hcount[std::make_pair(QQ2, QQ1)];
	    //Hmid[std::make_pair(QQ2, QQ1)] = MBG.is_T_consistent_color(im->second);
	  }
	} 
      } 
    }
    processed.insert(*is);
  }

   return Hcount;	
} 

void Statistics::count_cycles() { 
  std::unordered_set<std::string> processed;

  for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {
    if (processed.find(*is) != processed.end()) { 
      continue; 
    } 

    //multimularcs_t Mx = graph.get_adjacent_multiedges_v2(*is);
    mularcs_t Mx = graph.get_adjacent_multiedges(*is); 
    if (!graph.is_complement_color(Mx)) { 
      continue;
    } 

    std::string current = *is;
    std::string prev = "";
    Mcolor special_Q; 

    do {
      processed.insert(current);
      //multimularcs_t My = graph.get_adjacent_multiedges_v2(current);
      mularcs_t My = graph.get_adjacent_multiedges(current);
      if (!graph.is_simple_vertice(My)) {
	break;
      }

      if (prev == My.begin()->first) {
	prev = current;
	current = My.rbegin()->first;
      } else {
	prev = current;
	current = My.begin()->first;
      }
               
      while (current == Infty) {
	if (special_Q.empty()) {
	  special_Q = My.find(current)->second;
	  prev = *is;
	  current = Mx.begin()->first; 
	} else {
	  if (special_Q != My.find(current)->second) { 
	    ++special_cycle_count[std::min(special_Q, My.find(current)->second)]; 	  
	  }
	  break;
	}
      }
    } while ((current != Infty) && (processed.find(current) == processed.end()));
	
    if (current == *is) { //find cycle. 
      ++simple_cycle_count[std::min(Mx.begin()->second, Mx.rbegin()->second)];
    }
  }
} 

void Statistics::count_chromosomes() { 
  circular_chr.resize(graph.size_graph());
  liniar_chr.resize(graph.size_graph());
  
  for(int i = 0; i < graph.size_graph(); ++i) {
    std::unordered_set<orf_t> processed;		

    for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {		    
      if (processed.find(*is) == processed.end()) { 				  
	processed.insert(*is);
	std::string y = graph.get_adj_vertex(*is);

	while (true) {
	  if (member(processed, y)) {
	    ++circular_chr[i];
	    break;
	  }
	  
	  processed.insert(y);
	  if (!graph.is_there_edge(i, y)) {
	    ++liniar_chr[i];
	    break;
	  }
	  
	  y = graph.get_adj_vertex(i, y);
	  if (member(processed, y)) {
	    ++circular_chr[i];
	    break;
	  }
	  processed.insert(y);
	  y = graph.get_adj_vertex(y);
	}
		
	if (graph.is_there_edge(i, *is)) {
	  std::string y = graph.get_adj_vertex(i, *is);
					
	  while (processed.find(y) == processed.end()) {
	    processed.insert(y);
	    y = graph.get_adj_vertex(y);
	    if (member(processed, y)) { 
	      break;
	    } 

	    processed.insert(y);
	    if (!graph.is_there_edge(i, y))  { 
	      break;
	    } 
	    y = graph.get_adj_vertex(i, y);
	  }
	}		    
      }	
    } 
  }
} 
