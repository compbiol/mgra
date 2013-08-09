/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008 - 2013 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#include "reader.h"
#include "algo/Algorithm.h"
#include "Wgenome.h"

std::vector<std::string> genome_match::number_to_genome;
genome_match::gen2num genome_match::genome_to_number;   

typedef std::list<TwoBreak<Mcolor> > transform_t;

#include "RecoveredGenomes.h"

/* Given a non-linear genome PG of a multicolor Q and a transformation into a linear genome,
 * find linearizing fissions in the transformation, move them to beginning, and apply to PG
 * i.e., try to reorder transformation into: PG -> PG' -> linear genome, where PG' is linear
 * and obtained from PG with fission.
 * We replace PG with PG' and return the transformation PG -> PG'
 * Transformation may contain only multicolors Q' with Q'\cap Q = 0 or Q.
*/
transform_t decircularize(mbgraph_with_history<Mcolor>& graph, partgraph_t& PG, transform_t& TG, const Mcolor& Q) {
    RecoveredGenomes<mbgraph_with_history<Mcolor> > reductant(graph);	
	
    // decircularizing sub-transform that is removed
    transform_t D;

    size_t CircSize = reductant.numchr(PG).second;
    if (CircSize == 0) {
	return D;
    } 

    //std::cerr << "Eliminating " << CircSize << " circular chromosomes in " << genome_match::mcolor_to_name(Q) << std::endl;

    partgraph_t T = PG; // reconstructed genome ("bad")

    auto start = TG.begin();

    // looking for de-circularizig 2-breaks
    for(auto it = start; it != TG.end();) {
        // check multicolor
	{
	    Mcolor C(it->get_mcolor(), Q, Mcolor::Intersection);
    
	    if (C.empty()) {
		++it;
		continue;
	    }
    
	    if (C != Q) {
		//std::cerr << "Impossible multicolor in the transformation!" << std::endl;
		break;
	    }
	}

	it->apply_single(T);

        size_t ccsize = reductant.numchr(T).second;

	if( ccsize >= CircSize) {
	    ++it;
	    continue;
	}

	//std::cerr << "Found problematic 2-break: ";// << *it << "\t";

	// move t over to beginning
	for(auto jt=it;jt!=TG.begin();) {

	    auto kt = jt--; // jt, kt are successive, *kt == t

	    const TwoBreak<Mcolor>& t = *kt;
	    const TwoBreak<Mcolor>& s = *jt;

//            outlog << "... trying to swap with " << s << endl;

	    std::pair<vertex_t,vertex_t> p1, q1, p2, q2;

	    bool usearc = false;

	    Mcolor C(t.get_mcolor(), s.get_mcolor(), Mcolor::Intersection);
	    if (!C.empty()) {


		/*
			 p1=(x1,x2) x (y1,y2)=q1
			 p2=(x1,y1) x (x3,y3)=q2
    
			 into:
    
			 (x1,x2) x (x3,y3)
			 (y3,x2) x (y1,y2)
		*/
    
		for(int j = 0; j < 2; ++j) {    
		    if (t.get_arc(j) == std::make_pair(jt->get_arc(0).first, jt->get_arc(1).first)) { 
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = jt->get_arc(0);
			q1 = jt->get_arc(1);
		    } else if (t.get_arc(j) == std::make_pair(jt->get_arc(1).first, jt->get_arc(0).first)) {
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = jt->get_arc(1);
			q1 = jt->get_arc(0);
		    } else if (t.get_arc(j) == std::make_pair(jt->get_arc(0).second, jt->get_arc(1).second)) {
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = std::make_pair(jt->get_arc(0).second, jt->get_arc(0).first);
			q1 = std::make_pair(jt->get_arc(1).second, jt->get_arc(1).first);
		    } else if (t.get_arc(j) == std::make_pair(jt->get_arc(1).second, jt->get_arc(0).second)) {
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = std::make_pair(jt->get_arc(1).second, jt->get_arc(1).first);
			q1 = std::make_pair(jt->get_arc(0).second, jt->get_arc(0).first);
		    }
		    if (usearc) break;
		}
	    }

	    // TwoBreak t0 = t;

	    if (usearc) {
		if (t.get_mcolor() != s.get_mcolor()) break;
		*kt = TwoBreak<Mcolor>(q2.second, p1.second, q1.first, q1.second, t.get_mcolor());
		*jt = TwoBreak<Mcolor>(p1.first, p1.second, q2.first, q2.second, t.get_mcolor());
	    } else {
		TwoBreak<Mcolor> temp = *kt;
		*kt = *jt;
                *jt = temp;
	    }

	    {
		Mcolor C(kt->get_mcolor(), Q, Mcolor::Intersection);
    
                // N.B. at this point if C is not empty, then C == Q
		if( !C.empty() ) {
		    kt->inverse().apply_single(T);

		    ccsize = reductant.numchr(T).second;
		}
	    }

	    /*
	    if( CC.size() > ccsize ) {
		outlog << "Cannot pop:" << endl;
		outlog << *jt << " , " << t0 << "  -->  " << t << " , " << *kt << endl;
	    }
	    */

	}


	if( ccsize < CircSize ) {
	    //std::cerr << " SUCCEDED" << std::endl;

	    // move t away from the transformation TG and save it to D
            TG.begin()->apply_single(PG);
	    D.push_back(*TG.begin());

	    TG.erase(TG.begin());

	    CircSize = reductant.numchr(PG).second;

	    if( CircSize==0 ) break;

	    start = TG.begin();
	}
	else {  // did not succeed
	    start++;
	    //std::cerr << " FAILED" << std::endl;
	}

	T = PG;
	for(it = TG.begin(); it != start; ++it) {
	    it->apply_single(T);
	}
    }
   
    //if (CircSize > 0) {
	//std::cerr << "Unable to de-circularize ;-(" << std::endl;
    //}

    return D;
}

void tell_root_besides(const mbgraph_with_history<Mcolor>& graph) {
  // tell where the root resides
  std::clog << "the root resides in between:";

  std::set<Mcolor> T(graph.cbegin_T_consistent_color(), graph.cend_T_consistent_color()); 

  for (auto it = T.begin(); it != T.end(); ++it) {
    for (auto jt = it; ++jt != T.end(); ) {
      Mcolor C(*it, *jt, Mcolor::Intersection);
      if (C.size() == it->size()) {
	T.erase(it++);
	jt = it;
	continue;
      }
      if (C.size() == jt->size()) {
	T.erase(jt++);
	--jt;
      }
    }
    std::clog << " " << genome_match::mcolor_to_name(*it);
  }
  std::clog << std::endl;
}


////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
  std::cout << "MGRA (Multiple Genome Rearrangements & Ancestors) ver. 1.5" << std::endl;
  std::cout << "(c) 2008,12 by Max Alekseyev <maxal@cse.sc.edu>" << std::endl;
  std::cout << "Distributed under GNU GENERAL PUBLIC LICENSE license." << std::endl;
  std::cout << std::endl;

  /*reading flags*/
  std::string name_cfg_file;
  if (argc != 2) {
    std::cerr << "Usage: mgra <ProblemConfiguration>" << std::endl;
    return 1;
  } else { 
    name_cfg_file = argv[1];
  } 

  /*Reading problem configuration*/
  ProblemInstance<Mcolor> PI(reader::read_cfg_file(name_cfg_file)); 
  std::vector<Genome> genomes = reader::read_genomes(PI);
  
  genome_match::init_name_genomes(PI, genomes); //FIXME: IT'S DEBUG

  for(size_t i = 0; i < genomes.size(); ++i) { 
	std::clog << "Genome " << PI.get_priority_name(i) << " blocks: " << genomes[i].size() << std::endl;
  } 

  std::shared_ptr<mbgraph_with_history<Mcolor> > graph(new mbgraph_with_history<Mcolor>(genomes, PI)); 

  std::clog << "vecT-consistent colors: " << graph->count_vec_T_consitent_color() << std::endl;
  for (auto id = graph->cbegin_T_consistent_color(); id != graph->cend_T_consistent_color(); ++id) {
    std::clog << "\t" << genome_match::mcolor_to_name(*id);  ///FIXME: CHANGE
  }
  std::clog << std::endl;

  tell_root_besides(*graph); 	

  Algorithm<mbgraph_with_history<Mcolor> > main_algo(graph);

  main_algo.convert_to_identity_bgraph(PI); 
 
  if (!PI.get_target().empty()) {
	partgraph_t PG;

	for(const auto &x : *graph) {
	    if (PG.defined(x)) { 
		continue;
	    }

	    vertex_t y = Infty;
	    bool good = true;
	    size_t def = 0;
	    const auto &target = PI.get_target();
	    for (const auto &col : target) { 
		size_t numb_color = col.first; 
		if (graph->is_exist_edge(numb_color, x)) {
		    ++def;
		    if (y == Infty) { 
			y = graph->get_adjecent_vertex(numb_color, x);
		    } 

		    if (y != graph->get_adjecent_vertex(numb_color, x)) { 
			good = false;
		    }
		}
	    }

	    if (good && def == target.size() && y != Infty) {
                PG.insert(x, y);
	    }
	}

        RecoveredGenomes<mbgraph_with_history<Mcolor> > reductant(*graph);	
	std::set<std::pair<path_t, bool> > GN;
	reductant.splitchr(PG, GN, reductant.pg_empty);

	writer::Wgenome<std::set<std::pair<path_t, bool> > > writer_genome; 
	writer_genome.save_genom_in_text_format(genome_match::mcolor_to_name(PI.get_target()), GN, PI.get_target().empty());

    } else {  /* empty target */
        for(auto it = graph->cbegin_local_graphs(); it != graph->cend_local_graphs() - 1; ++it) { 
	  if (*it != *(it + 1)) {
	    std::cout << "T-transformation is not complete. Cannot reconstruct genomes." << std::endl; //FIXME clog
	    exit(1);
	  }
        }

        std::vector<partgraph_t> RG(graph->count_vec_T_consitent_color(), *(graph->cbegin_local_graphs())); // recovered genomes
       
        RecoveredGenomes<mbgraph_with_history<Mcolor> > reductant(*graph);
 	reductant.main_algorithm(RG);
    
	// T-transformation complete, we procede with recovering the ancestral genomes
	//std::cerr << "Initial 2-break distances from the root X: " << std::endl;
	// FIXME: check that the order in which circular chromosomes are eliminated    
        std::vector<transform_t> RT = graph->get_vec_TC_events(); // recovered transformations
	size_t i = 0; 
	for (auto im = graph->cbegin_T_consistent_color(); im != graph->cend_T_consistent_color(); ++im, ++i) {    
	    transform_t T = decircularize(*graph, RG[i], RT[i], *im);
    
	    // move to adjacent branches
	    for (const auto &it : T) {
		size_t j = 0; 
		for (auto jt = graph->cbegin_T_consistent_color(); jt != graph->cend_T_consistent_color(); ++jt, ++j) {
		    if ((j != i) && includes(im->cbegin(), im->cend(), jt->cbegin(), jt->cend()) && graph->are_adjacent_branches(*im, *jt)) {
			RT[j].push_back(it);
		    }
		}
	    }
	}

    	//std::cerr << "Final 2-break distances from the root X: " << std::endl;
      
	i = 0; 
	for (auto im = graph->cbegin_T_consistent_color(); im != graph->cend_T_consistent_color(); ++im, ++i) {
    	    std::set<std::pair<path_t, bool> > GN;
	    reductant.splitchr(RG[i], GN, reductant.pg_empty);

	    writer::Wgenome<std::set<std::pair<path_t, bool> > > writer_genome; 
	    writer_genome.save_genom_in_text_format(genome_match::mcolor_to_name(*im), GN, PI.get_target().empty());
        
	    std::ofstream tr((genome_match::mcolor_to_name(*im) + ".trs").c_str());
	    for(const auto &event : RT[i]) {
		tr << event.get_arc(0).first << " " << event.get_arc(0).second << "\t" 
		   << event.get_arc(1).first << " " << event.get_arc(1).second << "\t" 
	  	   << genome_match::mcolor_to_name(event.get_mcolor()) << std::endl;
	    }
	    tr.close(); 
	}
    }

    return 0;
}

