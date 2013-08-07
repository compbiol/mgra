#ifndef STAGE_3_H_
#define STAGE_3_H_

// cut the graph into connected components
template<class graph_t>
bool Algorithm<graph_t>::stage5_1() { 
 bool isChanged = false;
 size_t number_rear = 0; // number of rearrangements 

 do {
   number_rear = 0;
   //std::cerr << "Stage 5_1: splitting into connected components" << std::endl;

   // go over all T-consistent multicolors
   for(auto ic = graph.cbegin_T_color(); ic != graph.cend_T_color(); ++ic) {
     const Mcolor& Q = *ic;

     bool repeat = true;
     while(repeat) {
       repeat = false;

       equivalence<vertex_t> CC; // connected components
       std::map<vertex_t, vertex_t> QQ; // multiedges of colors !Q (!*ic)
		    
       for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {    
	 Mularcs<Mcolor> M = graph.get_adjacent_multiedges(*is);

	 if (M.size() == 1 && M.cbegin()->second == graph.get_complete_color()) { 
	   continue; // ignore complete multiedges
	 } 

	 for(auto im = M.cbegin(); im != M.cend(); ++im) {    
	   if (im->second == *ic) { // edges of color Q (*ic)
	     QQ.insert(std::make_pair(*is, im->first));
	   } else if (im->first != Infty) { // reg. edges of color !Q (!*ic)
	     CC.addrel(*is, im->first);
	   }
	 }
       }
		    
       // N.B. for regular edges (x, y) of multicolor Q, QQ[x] = y and QQ[y] = x
       // for irregular edge (x,oo), QQ[x] = oo		    
       CC.update();
       //std::cerr << genome_match::mcolor_to_name(*ic) << " ~ " << CC.classes() << std::endl;
    
       typedef std::string concom_t; //connected components type
       std::map<concom_t, std::set<arc_t> > EC; // reg. edges between diff. connected components of color Q
       std::map<concom_t, std::set<arc_t> > EI; // irreg. edges of color Q
    
       for(auto iq = QQ.begin(); iq != QQ.end(); ++iq) {
	 if (iq->second == Infty) {
	   EI[CC[iq->first]].insert(std::make_pair(iq->first, iq->second));
	 } else if (CC.isequiv(iq->first, iq->second)) { 
	   continue;
	 } else {			
	   EC[CC[iq->first]].insert(std::make_pair(iq->first, iq->second));
	 }
       }
		
       std::set<concom_t> processed;

       // reg. edges between diff. connected components of color Q
       for(auto ie = EC.begin(); ie != EC.end(); ++ie) {
	 if (processed.find(ie->first) != processed.end()) { 
		continue;
	 }

	 // connected component with a single external edge
	 if (ie->second.size() == 1) {
	   const arc_t& p = *(ie->second.begin());
	   arc_t q;

	   if (graph.is_indel_vertex(p.first) || graph.is_indel_vertex(p.second) 
		|| graph.is_duplication_vertex(p.first) || graph.is_duplication_vertex(p.second)) {
		continue;
           }

	   // look for irregular edges
	   bool found = false;
    
	   for(auto ii = EI[ie->first].cbegin(); (ii != EI[ie->first].cend()) && !found; ++ii) {
	     const arc_t& t = *ii; // edge

	     // let check what would happen with (be-)edge e=(p.first,t.first)

	     Mcolor T = Q;
	     for(size_t i = 0; i < graph.count_local_graphs(); ++i) {
	       if (graph.is_exist_edge(i, p.first) && graph.get_adjecent_vertex(i, p.first) == t.first) {
		 T.insert(i);
	       }
	     } 
	     // if e is enriched to T-consistent color, great!
	     if (T.size() > Q.size() && graph.is_T_consistent_color(T) ) {
	       //std::cerr << "perfect edge is found" << std::endl;
	       q = t;
	       found = true;
	       EI[ie->first].erase(ii);
	     }
	   }

	   // we did not find good edge, but there are some candidates - take any
	   if (!found && EI[ie->first].size() == 1) {
	     //std::cerr << "somewhat good but unique edge is found" << std::endl;
	     q = *(EI[ie->first].begin());
	     EI.erase(ie->first);
	     found = true;
	   }

	   if (!found && EI[ie->first].size() == 0) {
	     //std::cerr << "no irregular edges, do fission" << std::endl;
	     q = std::make_pair(Infty, Infty);
	     found = true;
	   }
    
	   if (!found) {
		continue;
           }

	   if ((q.first != Infty && (graph.is_indel_vertex(q.first) || graph.is_duplication_vertex(q.first)))  
		|| (q.second != Infty && (graph.is_indel_vertex(q.second) || graph.is_duplication_vertex(q.second)))) {
		continue;
           }

        
	   graph.apply_two_break(TwoBreak<Mcolor>(p, q, Q));
	   ++number_rear;
			    
	   //std::cerr << "Stage 5_1: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << std::endl;

	   processed.insert(CC[p.first]);
	   if (q.first != Infty) { 
	     processed.insert(CC[q.first]);
	   } 

	   processed.insert(CC[p.second]);
	   if (q.second != Infty) { 
	     processed.insert(CC[q.second]);
	   } 

	   repeat = true;
	 } else if (ie->second.size() == 2 && EI[ie->first].size() == 0) {
	   arc_t p = *(ie->second.begin());
	   arc_t q = *(ie->second.rbegin());

	   if (graph.is_indel_vertex(p.first) || graph.is_indel_vertex(p.second) 
		|| graph.is_duplication_vertex(p.first) || graph.is_duplication_vertex(p.second)) {
		continue;
           }

	   if (graph.is_indel_vertex(q.first) || graph.is_indel_vertex(q.second) 
		|| graph.is_duplication_vertex(q.first) || graph.is_duplication_vertex(q.second)) {
		continue;
           }
			    
	   // N.B. we have CC[p.first] == CC[q.first] == ie->first
	   if ((processed.count(CC[p.second]) != 0) || (processed.count(CC[q.second]) != 0) || CC[p.second] != CC[q.second]) { 
		continue;
	   }

	   graph.apply_two_break(TwoBreak<Mcolor>(p, q, Q));
	   ++number_rear;

	   //std::cerr << "Stage 222.2: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << endl;

	   processed.insert({CC[p.first], CC[q.first], CC[p.second], CC[q.second]});
	   repeat = true;
	 }
       }
     }
   }

   if (number_rear != 0) { 
     isChanged = true;
   } 
 } while (number_rear > 0); 
 
 return isChanged; 
} 

//search for 4-cycles
template<class graph_t>
bool Algorithm<graph_t>::stage5_2() { 
 bool isChanged = false;
 size_t number_rear = 0; // number of rearrangements

 do {
   number_rear = 0;

   for(const auto &x : graph) {
     if (graph.is_indel_vertex(x) || graph.is_duplication_vertex(x)) {
	continue;
     }

     Mularcs<Mcolor> Mx = graph.get_adjacent_multiedges(x);
     bool next = false;

     for(auto im = Mx.cbegin(); (im != Mx.cend()) && (!next); ++im) {
       const vertex_t& y = im->first;
       const Mcolor& Q = im->second;

       if (graph.is_vec_T_color(Q) && y != Infty && !graph.is_indel_vertex(y) && !graph.is_duplication_vertex(y)) { 
         Mularcs<Mcolor> My = graph.get_adjacent_multiedges(y);
         My.erase(x);

	 for(auto jm = My.cbegin(); (jm != My.cend()) && (!next); ++jm) {
	   const vertex_t& z = jm->first;
	  
	   if (z != Infty) { 
	     Mularcs<Mcolor> Cz = graph.get_adjacent_multiedges(z);
	     vertex_t v = Cz.get_vertex(Q);
	     if (!v.empty() && Mx.defined(v)) { 
	       //std::cerr << "Stage 5_2: " << x << " - " << y << "\tX\t" << v << " - " << z << std::endl;
	       graph.apply_two_break(TwoBreak<Mcolor>(x, y, v, z, Q));
	       ++number_rear;
	       next = true;
	     }
	   }
         }
       }
     }
     
   }
 
   if (number_rear != 0) { 
    isChanged = true; 
   } 
 } while (number_rear > 0);

 return isChanged; 
} 
#endif
