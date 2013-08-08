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

std::vector<std::string> genome_match::number_to_genome;
genome_match::gen2num genome_match::genome_to_number;   

std::ofstream outlog("/dev/null");

typedef std::list<TwoBreak<Mcolor> > transform_t;

std::set<vertex_t> getchrset;
std::pair<path_t, bool> getchr(const mbgraph_with_history<Mcolor>& graph, const partgraph_t& PG, const std::string& x) {
    path_t path;
    bool circular = false;

    getchrset.clear();
    getchrset.insert(x);

    for(vertex_t y = graph.get_obverse_vertex(x); ; ) {
	if (getchrset.find(y) != getchrset.end()) {
	    circular = true;
	    break; // circ
	}
	getchrset.insert(y);

	{
	    vertex_t xx = y;
	    xx.resize(xx.size() - 1);
	    if (*y.rbegin() == 't') {
		path.push_back("-" + xx);
	    } else {
		path.push_back("+" + xx);
	    }
	}

	if (!PG.defined(y)) { 
	  break; // linear
	} 

	y = PG[y];

	if (getchrset.find(y) != getchrset.end()) {
	    circular = true;
	    break; // circ
	}
	getchrset.insert(y);

	if (y == Infty) { 
		break;
	} 
	y = graph.get_obverse_vertex(y);
    }

    if (!circular && PG.defined(x)) {
	for(vertex_t y = x; PG.defined(y);) {
	    y = PG[y];
	    getchrset.insert(y);

	    if (y == Infty) {
		break;
	    }
	    y = graph.get_obverse_vertex(y);
	    getchrset.insert(y);
	    {
		std::string xx = y;
		xx.resize(xx.size() - 1);
		if (*y.rbegin() == 't') {
		    path.push_front("+" + xx);
		} else {
		    path.push_front("-" + xx);
		}
	    }
	}
    }

    return std::make_pair(path, circular);
}


std::list< std::set<vertex_t> > pg_empty;
void splitchr(const mbgraph_with_history<Mcolor>& graph, const partgraph_t& PG, std::set<std::pair<path_t, bool> >& AllChr, std::list<std::set<vertex_t> >& CircChr = pg_empty) {

    if (&CircChr != &pg_empty) { 
	CircChr.clear();
    } 
    AllChr.clear();
    std::unordered_set<vertex_t> processed;

    for(const auto &x : graph) { 
	if (processed.find(x) == processed.end()) { 
	        std::pair<path_t, bool> pathb = getchr(graph, PG, x);
	
		AllChr.insert(pathb);

	        std::copy(getchrset.begin(), getchrset.end(), std::inserter(processed, processed.end()));
	
		if (pathb.second && (&CircChr != &pg_empty)) {
		    CircChr.push_back(getchrset);
		}
	} 
    }
}

std::pair<size_t, size_t> numchr(const mbgraph_with_history<Mcolor>& graph, const partgraph_t& PG) {
    std::set<std::pair<path_t, bool> > AllChr;
    std::list<std::set<vertex_t> > CircChr;
    splitchr(graph, PG, AllChr, CircChr);
    return std::make_pair(AllChr.size(), CircChr.size());
}

//rename and move to namespace writer
void printchr(const std::string& outname, const std::set<std::pair<path_t, bool> >& AllChr, bool isEmptyTarget) { 
	std::ofstream out((outname + ".gen").c_str());

	out << "# Genome " << outname << std::endl;

	std::string ChrTitle;

	if (isEmptyTarget) { 
		ChrTitle = "chromosome"; 
	} else { 
		ChrTitle = "CAR";
	} 

	size_t ncirc = 0; 
	size_t lcirc = 0;

	for(auto ia = AllChr.cbegin(); ia != AllChr.cend(); ++ia) {
		out << std::endl;

		const path_t& path = ia->first;

#ifdef HAVE_BLOCK_LENGTH
		size_t len = 0; // compute length in bp
#endif
		for(auto ip = path.begin(); ip != path.end(); ++ip) {
			std::string t = *ip;
			if (t[0] == '-' || t[0] == '+' ) { 
				t = t.substr(1);
			} 

#ifdef HAVE_BLOCK_LENGTH
			len += BL[t + 't'];
#endif
		}

		if (ia->second) {
			++ncirc;
			lcirc += path.size();
			out << "# circular ";
		} else {
			out << "# linear ";
		}

		out << ChrTitle << " of length " << path.size();

#ifdef HAVE_BLOCK_LENGTH
		out << " (" << len << " bp)";
#endif
		out << " follows:" << std::endl;

		if( (*path.begin())[0] == '+' || (*path.rbegin())[0] == '+') {
			for(auto ip = path.cbegin(); ip != path.cend(); ++ip) {
				out << *ip << " ";
			}
		} else {
			for(auto ip = path.rbegin(); ip != path.rend();++ip) {
				std::string e = *ip;
				if (e[0] == '-') { 
					e[0] = '+'; 
				} else if (e[0] == '+') { 
					e[0] = '-';
				} 
				out << e << " ";
			}
		}

		out << "$" << std::endl;
	}

	out << std::endl << "# Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << std::endl;
	std::cout << std::endl << "Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << std::endl;

	if (ncirc) {
		out << "#\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl;
		std::cout << "\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl;
	}

	out.close();
}

// fill in OP container with endpoints of q-obverse paths,
// starting and ending at OP
void get_obverse_paths(const mbgraph_with_history<Mcolor>& graph, std::map<vertex_t, std::set<vertex_t> >& OP, const Mcolor Q) {
    std::map<vertex_t, std::set<int> > processed;

    for(auto iq = Q.cbegin(); iq != Q.cend(); ++iq) {
        const partgraph_t& PG = graph.get_local_graph(iq->first);

	for(auto ip = OP.begin(); ip != OP.end(); ++ip) {
            const vertex_t& x = ip->first;

	    if (x == Infty || (processed[x].find(iq->first) != processed[x].end())) {
		continue;
	    }

	    for(vertex_t y = graph.get_obverse_vertex(ip->first); PG.defined(y);) {
		if (OP.find(y) != OP.end()) {
		    ip->second.insert(y);
		    OP[y].insert(x);
                    processed[y].insert(iq->first);
		    break;
		}
		y = PG[y];

		if( y==x ) {
		    //already is circular, we don't care
		    break;
		}

		y = graph.get_obverse_vertex(y);
	    }
	    processed[x].insert(iq->first);
	}
    }
}

/* Given a non-linear genome PG of a multicolor Q and a transformation into a linear genome,
 * find linearizing fissions in the transformation, move them to beginning, and apply to PG
 * i.e., try to reorder transformation into: PG -> PG' -> linear genome, where PG' is linear
 * and obtained from PG with fission.
 * We replace PG with PG' and return the transformation PG -> PG'
 * Transformation may contain only multicolors Q' with Q'\cap Q = 0 or Q.
*/
transform_t decircularize(mbgraph_with_history<Mcolor>& graph, partgraph_t& PG, transform_t& TG, const Mcolor& Q) {

    // decircularizing sub-transform that is removed
    transform_t D;

    size_t CircSize = numchr(graph, PG).second;
    if( CircSize == 0 ) return D;

    outlog << "Eliminating " << CircSize << " circular chromosomes in " << genome_match::mcolor_to_name(Q) << std::endl;

    /*
    transform_t p1,p2;
    for(transform_t::iterator it=TG.begin();it!=TG.end();++it) {
	if( it->MultiColor==Q ) {
	    p1.push_back(*it);
	}
	else {
	    p2.push_back(*it);
	}
    }
    p1.insert(p1.end(),p2.begin(),p2.end());
    TG = p1;
    outlog << "Transformation reordered" << endl;
    */

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
		outlog << "Impossible multicolor in the transformation!" << std::endl;
		break;
	    }
	}

	it->apply_single(T);

        size_t ccsize = numchr(graph, T).second;

	if( ccsize >= CircSize) {
	    ++it;
	    continue;
	}

	outlog << "Found problematic 2-break: ";// << *it << "\t";

	// move t over to beginning
	for(auto jt=it;jt!=TG.begin();) {

	    auto kt = jt--; // jt, kt are successive, *kt == t

	    const TwoBreak<Mcolor>& t = *kt;
	    const TwoBreak<Mcolor>& s = *jt;

//            outlog << "... trying to swap with " << s << endl;

	    std::pair<vertex_t,vertex_t> p1, q1, p2, q2;

	    bool usearc = false;

	    Mcolor C(t.get_mcolor(), s.get_mcolor(), Mcolor::Intersection);
	    if( !C.empty() ) {


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
    
			p1 = make_pair(jt->get_arc(0).second, jt->get_arc(0).first);
			q1 = make_pair(jt->get_arc(1).second, jt->get_arc(1).first);
		    } else if (t.get_arc(j) == std::make_pair(jt->get_arc(1).second, jt->get_arc(0).second)) {
			usearc = true;
    
			p2 = t.get_arc(j);
			q2 = t.get_arc(1 - j);
    
			p1 = make_pair(jt->get_arc(1).second, jt->get_arc(1).first);
			q1 = make_pair(jt->get_arc(0).second, jt->get_arc(0).first);
		    }
		    if(usearc) break;
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

		    ccsize = numchr(graph, T).second;
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
	    outlog << " SUCCEDED" << std::endl;

	    // move t away from the transformation TG and save it to D
            TG.begin()->apply_single(PG);
	    D.push_back(*TG.begin());

	    TG.erase(TG.begin());

	    CircSize = numchr(graph, PG).second;

	    if( CircSize==0 ) break;

	    start = TG.begin();
	}
	else {  // did not succeed
	    start++;
	    outlog << " FAILED" << std::endl;
	}

	T = PG;
	for(it = TG.begin(); it != start; ++it) {
	    it->apply_single(T);
	}
    }
    //if( start == TG.end() ) {
    if( CircSize>0 ) {
	outlog << "Unable to de-circularize ;-(" << std::endl;
    }

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

std::vector<partgraph_t> RG; // recovered genomes
std::vector<transform_t> RT; // and transformations

///////////////////////////////////////////////////////////////////////////
void RecoverGenomes(const mbgraph_with_history<Mcolor>& graph) {
    size_t NC = graph.count_vec_T_consitent_color();
    RG.resize(NC);
    RT.resize(NC);

    for(size_t i = 0; i < NC; ++i) { 
	RG[i] = *(graph.cbegin_local_graphs());
    } 

    // track changes in the number of chromosomes

    // number of reversals, interchromosomal translocations, and fissions/fusions
    /*std::vector< std::vector<size_t> > RTF(NC);
    for(size_t j = 0; j < NC; ++j) { 
	RTF[j].resize(3);
    } */

    for(auto it = graph.crbegin_2break_history(); it != graph.crend_2break_history(); ++it) {
	const Mcolor& Q = it->get_mcolor();

	//std::cerr << "Reverting (" << it->get_arc(0).first << "," << it->get_arc(0).second << ")x(" << it->get_arc(1).first << "," << it->get_arc(1).second << "):{" << genome_match::mcolor_to_name(Q) << "} " << " in";

	size_t i = 0;
	for(auto im = graph.cbegin_T_consistent_color(); im != graph.cend_T_consistent_color(); ++im) {
    	
	    if (!Q.includes(*im)) { 
	        ++i;
		continue;
	    } 

            // TColor[i] is subset of Q

            size_t nchr_old = 0;
	    if (Q == *im) {
		nchr_old = numchr(graph, RG[i]).first;
	    }

	    it->inverse().apply_single(RG[i]);

	    if (Q == *im) {
		//std::cerr << " " << genome_match::mcolor_to_name(*im);
		RT[i].push_front(*it);
	    }

	    if (Q == *im) {
		bool samechr = true;

		std::unordered_set<vertex_t> vert;
		if (it->get_arc(0).first != Infty) vert.insert(it->get_arc(0).first);
		if (it->get_arc(0).second != Infty) vert.insert(it->get_arc(0).second);
		if (it->get_arc(1).first != Infty) vert.insert(it->get_arc(1).first);
		if (it->get_arc(1).second != Infty) vert.insert(it->get_arc(1).second);
    
		getchr(graph, RG[i], *vert.begin());
    
		vert.erase(vert.begin());
		for(const auto &iv : vert) {
		    if (getchrset.find(iv) == getchrset.end()) {
			samechr = false;
			break;
		    }
		}

		size_t nchr_new = numchr(graph, RG[i]).first;
		/*if (nchr_new != nchr_old) {
		    ++RTF[i][2];
		} else {
		    if (samechr) {
			++RTF[i][0];
		    } else { 
			++RTF[i][1];
		    } 
		}*/
	    }
	    ++i;
	}
	
	//std::cerr << std::endl;
    }

    /*std::vector<size_t> tot(3);
    std::cerr << "% Number of reversals / translocations / fissions+fusions: " << std::endl;
    for(size_t j = 0; j < NC; ++j) {
	tot[0] += RTF[j][0];
	tot[1] += RTF[j][1];
	tot[2] += RTF[j][2];
    }
    std::cerr << "Total\t&\t" << tot[0] << " & " << tot[1] << " & " << tot[2] << " &\t" << tot[0]+tot[1]+tot[2] << " \\\\" << std::endl;*/
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
 
#ifndef VERSION2  
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
	
	std::set< std::pair<path_t, bool> > GN;
	splitchr(*graph, PG, GN);
	printchr(genome_match::mcolor_to_name(PI.get_target()), GN, PI.get_target().empty());

    } else {  /* empty target */
        for(auto it = graph->cbegin_local_graphs(); it != graph->cend_local_graphs() - 1; ++it) { 
	  if (*it != *(it + 1)) {
	    std::cout << "T-transformation is not complete. Cannot reconstruct genomes." << std::endl; //FIXME clog
	    exit(1);
	  }
        }

	RecoverGenomes(*graph);
     
	// T-transformation complete, we procede with recovering the ancestral genomes
	//std::cerr << "Initial 2-break distances from the root X: " << std::endl;

	// FIXME: check that the order in which circular chromosomes are eliminated    
	size_t i = 0; 
	for (auto im = graph->cbegin_T_consistent_color(); im != graph->cend_T_consistent_color(); ++im) {    
	    transform_t T = decircularize(*graph, RG[i], RT[i], *im);
    
	    // move to adjacent branches
	    for (const auto &it : T) {
		size_t j = 0; 
		for (auto jt = graph->cbegin_T_consistent_color(); jt != graph->cend_T_consistent_color(); ++jt) {
		    if ((j != i) && includes(im->cbegin(), im->cend(), jt->cbegin(), jt->cend()) && graph->are_adjacent_branches(*im, *jt)) {
			RT[j].push_back(it);
		    }
		    ++j; 
		}
	    }
 	    ++i; 
	}
    
	//std::cerr << "Final 2-break distances from the root X: " << std::endl;
	
	i = 0; 
	for(auto im = graph->cbegin_T_consistent_color(); im != graph->cend_T_consistent_color(); ++im) {
    	    std::set<std::pair<path_t, bool> > GN;
	    splitchr(*graph, RG[i], GN);
	    printchr(genome_match::mcolor_to_name(*im), GN, PI.get_target().empty());
        
	    std::ofstream tr((genome_match::mcolor_to_name(*im) + ".trs").c_str());
	    for(const auto &event : RT[i]) {
		tr << event.get_arc(0).first << " " << event.get_arc(0).second << "\t" 
		   << event.get_arc(1).first << " " << event.get_arc(1).second << "\t" 
	  	   << genome_match::mcolor_to_name(event.get_mcolor()) << std::endl;
	    }
	    tr.close();
     	    ++i; 
	}
    }
#endif
    return 0;
}

