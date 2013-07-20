#ifndef STAGE_3_H_
#define STAGE_3_H_

// cut the graph into connected components
template<class graph_t>
bool Algorithm<graph_t>::stage3_1() { 
bool simplified = false;
 size_t nr = 0; // number of rearrangements, 

 do {
   nr = 0;
   outlog << "Stage 222: splitting into connected components" << std::endl;

   // go over all T-consistent multicolors
   for(auto ic = graph.cbegin_T_color(); ic != graph.cend_T_color(); ++ic) {
     if (ic->size() == 0 || ic->size() == graph.size_graph()) { //FIXME. size() == 0
       continue; // except empty and complete multicolor
     } 		
     const Mcolor& Q = *ic;

     bool repeat = true;
     while(repeat) {
       repeat = false;

       equivalence<vertex_t> CC; // connected components
       std::map<vertex_t, vertex_t> QQ; // multiedges of colors !Q (!*ic)
		    
       for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {    
	 Mularcs<Mcolor> M = graph.get_adjacent_multiedges(*is);

	 if (M.size() == 1) { 
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
       outlog << genome_match::mcolor_to_name(*ic) << " ~ " << CC.classes() << std::endl;
    
       typedef std::string concom_t;
       std::map <concom_t, std::set<std::pair<std::string, std::string> > > EC; // reg. edges between diff. connected components of color Q
       std::map <concom_t, std::set<std::pair<std::string, std::string> > > EI; // irreg. edges of color Q
    
       for(auto iq = QQ.begin(); iq != QQ.end(); ++iq) {
	 if (iq->second == Infty) {
	   EI[CC[iq->first]].insert(std::make_pair(iq->first, iq->second));
	   continue;
	 }
			
	 if (CC.isequiv(iq->first, iq->second)) { 
	   continue;
	 } 
			
	 EC[CC[iq->first]].insert(std::make_pair(iq->first, iq->second));
       }
		
       std::map <std::string, std::pair<std::string, std::string> > FE; // free end
       std::set <concom_t> processed;

       // reg. edges between diff. connected components of color Q
       for(auto ie = EC.begin(); ie != EC.end(); ++ie) {
	 if (member(processed, ie->first)) continue;

	 // connected component with a single external edge
	 if (ie->second.size() == 1) {
	   std::pair<std::string, std::string> p = *(ie->second.begin());
	   std::pair<std::string, std::string> q;

	   // look for irregular edges
	   bool found = false;
    
	   for(auto ii = EI[ie->first].cbegin(); ii != EI[ie->first].cend(); ++ii) {
	     const std::pair<std::string, std::string>& t = *ii; // edge

	     // let check what would happen with (be-)edge e=(p.first,t.first)

	     Mcolor T = Q;
	     for(size_t i = 0; i < graph.size_graph(); ++i) {
	       if( graph.get_local_graph(i).defined(p.first) && graph.get_local_graph(i)[p.first]==t.first ) {
		 T.insert(i);
	       }
	     } 
	     // if e is enriched to T-consistent color, great!
	     if( T.size()>Q.size() && graph.is_T_consistent_color(T) ) {
	       outlog << "perfect edge is found" << endl;
	       q = t;
	       found = true;
	       EI[ ie->first ].erase(ii);
	       break;
	     }
	   }

	   // we did not find good edge, but there are some candidates - take any
	   if( !found && EI[ ie->first ].size() == 1 ) {
	     outlog << "somewhat good but unique edge is found" << endl;
	     q = *(EI[ ie->first ].begin());
	     EI.erase( ie->first );
	     found = true;
	   }

	   // save for future
	   if( !found && EI[ ie->first ].size() == 0 ) {

	     /*
	       if( !member(FE,CC[p.second]) ) {
	       FE[ CC[p.second] ] = p;
	       }
	       else {
	       q = FE[ CC[p.second] ];
	       // N.B. we have CC[p.second] == CC[q.second]

	       FE.erase( CC[p.second] );
	       found = true;
	       }
	     */

	     outlog << "no irregular edges, do fission" << endl;
	     q = make_pair(Infty,Infty);
	     found = true;
	   }
    
	   if( !found ) continue;
    
    
	   if( ! TwoBreak<graph_t, Mcolor>(p, q, Q).apply(graph, true) ) continue;
	   ++nr;
			    
	   outlog << "Stage 222.1: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << std::endl;

	   processed.insert(CC[p.first]);
	   if (q.first != Infty) { 
	     processed.insert(CC[q.first]);
	   } 

	   processed.insert(CC[p.second]);
	   if (q.second != Infty) { 
	     processed.insert(CC[q.second]);
	   } 

	   /*
	     EC[ CC[p.second] ].erase( make_pair(p.second,p.first) );
	     if( q.second != Infty ) EC[ CC[q.second] ].erase( make_pair(q.second,q.first) );
	     EC[ CC[p.first] ].erase( make_pair(p.first,p.second) );
	     EC[ CC[q.first] ].erase( make_pair(q.first,q.second) );


	     if( !CC.isequiv(p.first,q.first) ) {
	     EC[ CC[p.first] ].insert( make_pair(p.first,q.first) );
	     EC[ CC[q.first] ].insert( make_pair(q.first,p.first) );
	     }
	   */

	   repeat = true;
	   continue;
	 }

	 if( ie->second.size()==2 && EI[ ie->first ].size() == 0 ) {

	   pair<string,string> p = *(ie->second.begin());
	   pair<string,string> q = *(ie->second.rbegin());
			    
	   // N.B. we have CC[p.first] == CC[q.first] == ie->first

	   if( member(processed,CC[p.second]) || member(processed,CC[q.second]) || CC[p.second]!=CC[q.second] ) continue;

	   if( ! TwoBreak<graph_t, Mcolor>(p,q,Q).apply(graph,true) ) continue;
	   nr++;

	   outlog << "Stage 222.2: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << endl;

	   processed.insert( CC[p.first] );
	   processed.insert( CC[q.first] );
	   processed.insert( CC[p.second] );
	   processed.insert( CC[q.second] );

	   /*
	     EC[ CC[p.second] ].erase( make_pair(p.second,p.first) );
	     EC[ CC[q.second] ].erase( make_pair(q.second,q.first) );
	     EC[ CC[p.first] ].erase( make_pair(p.first,p.second) );
	     EC[ CC[q.first] ].erase( make_pair(q.first,q.second) );
    
	     if( !CC.isequiv(p.second,q.second) ) {
	     EC[ CC[p.second] ].insert( make_pair(p.second,q.second) );
	     EC[ CC[q.second] ].insert( make_pair(q.second,p.second) );
	     }
	   */

	   repeat = true;
	   continue;
	 }
       }
     }
   }
   if (nr != 0) { 
     simplified = true;
   } 
 } while (nr > 0); 
 
 return simplified; 
} 

 // search for 4-cycles
template<class graph_t>
bool Algorithm<graph_t>::stage3_2() { 
 bool simplified = false;
 size_t nr = 0; // number of rearrangements, 

 do {
   nr = 0;

   for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {
     const string& x = *is;
     Mularcs<Mcolor> Mx = graph.get_adjacent_multiedges(x);

     bool next = false;

     for(auto im = Mx.cbegin(); im!=Mx.cend(); ++im) {

       const std::string& y = im->first;
       const Mcolor& Q = im->second;

       if (!graph.is_vec_T_color(Q) || y==Infty) { 
	 continue;
       }

       Mularcs<Mcolor> My = graph.get_adjacent_multiedges(y);
       My.erase(x);

       for(auto jm = My.cbegin(); jm != My.cend(); ++jm) {
	 std::string z = jm->first;
	 if (z == Infty) { 
	   continue;
	 } 
	
	 Mularcs<Mcolor> Cz = graph.get_adjacent_multiedges(z);

	 vertex_t v = "";
	 for (auto ir = Cz.cbegin(); ir != Cz.cend(); ++ir) {
	   if (ir->second == Q) { 
	     v = ir->first;
	   } 
	 }  

	 if ((!v.empty()) && (Mx.find(v) != Mx.cend())) {

	   auto p = std::make_pair(x, y);
	   auto q = std::make_pair(v, z);

	   outlog << "Stage 2222: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << std::endl;

	   if (TwoBreak<graph_t, Mcolor>(p, q, Q).apply(graph, true)) nr++;

	   next = true;

	   break;

	 }
       }
       if(next) break;
     }
   }

   if (nr != 0) { 
     simplified = true; 
   } 
 } while (nr > 0);

 return simplified; 
} 

#endif
