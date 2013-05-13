/* 
** Module: MGRA main body
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008,12 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <list>
#include <vector>
#include <array> 

#include <iterator>   // ostream_iterator etc.
using namespace std;

#include "reader.h"
#include "algo/Algorithm.h"
#include "writer/Wstats.h"
#include "writer/Wdots.h"

bool ConvPhylTreeAll(MBGraph&, int Stage);

#ifdef REMOVE_SHORT_BLOCKS
#define HAVE_BLOCK_LENGTH
#include "remove_short_blocks.h"
#endif

#ifdef BPG_COMPRESSED
#define HAVE_BLOCK_LENGTH
#endif

#ifdef HAVE_BLOCK_LENGTH
std::map<std::string, size_t> BL;  // block -> length FIXME
#endif

vector<partgraph_t> RG; // recovered genomes
vector<transform_t> RT; // and transformations

bool RecoverGenomes(MBGraph&, const transform_t&);
set <vertex_t> getchrset;

std::pair<path_t, bool> getchr(const MBGraph& graph, const partgraph_t& PG, const std::string& x) {
    path_t path;
    bool circular = false;

    getchrset.clear();
    getchrset.insert(x);

    for(vertex_t y = graph.get_adj_vertex(x); ; ) {
	if( member(getchrset,y) ) {
	    circular = true;
	    break; // circ
	}
	getchrset.insert(y);

	{
	    string xx = y;
	    xx.resize(xx.size()-1);
	    if( *y.rbegin()=='t' ) {
		path.push_back("-"+xx);
	    }
	    else {
		path.push_back("+"+xx);
	    }
	}

	if( !PG.defined(y) ) break; // linear

	y = PG[y];

	if( member(getchrset,y) ) {
	    circular = true;
	    break; // circ
	}
	getchrset.insert(y);

	y = graph.get_adj_vertex(y);
    }

    if( !circular && PG.defined(x) ) {
	for(string y = x;PG.defined(y);) {
	    y = PG[y];
	    getchrset.insert(y);

	    y = graph.get_adj_vertex(y);
	    getchrset.insert(y);
	    {
		string xx = y;
		xx.resize(xx.size()-1);
		if( *y.rbegin()=='t' ) {
		    path.push_front("+"+xx);
		}
		else {
		    path.push_front("-"+xx);
		}
	    }
	}
    }

    return make_pair( path, circular );
}


list< set<vertex_t> > pg_empty;
void splitchr(const MBGraph& graph, const partgraph_t& PG, set< pair<path_t,bool> >& AllChr, const bool Xonly = false, list< set<vertex_t> >& CircChr = pg_empty) {

    if (&CircChr != &pg_empty) { 
	CircChr.clear();
    } 
    AllChr.clear();
    std::set<orf_t> processed;

    for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {
	const string& x = *is;
	
	if( member(processed,x) ) continue;
       
        pair< path_t, bool > pathb = getchr(graph, PG, x);

	AllChr.insert( pathb );

        copy(getchrset.begin(),getchrset.end(),inserter(processed,processed.end()));

	if( pathb.second && (&CircChr != &pg_empty) ) {
	    CircChr.push_back(getchrset);
	}
    }
}

std::pair<size_t, size_t> numchr(const MBGraph& graph, const partgraph_t& PG) {
    set< pair<path_t,bool> > AllChr;
    list< set<vertex_t> > CircChr;
    splitchr(graph, PG, AllChr, false, CircChr);
    return make_pair(AllChr.size(),CircChr.size());
}

//rename and move to namespace writer
void printchr(const std::string& outname, const std::set<std::pair<path_t, bool> >& AllChr, bool isEmptyTarget) { 
	ofstream out((outname + ".gen").c_str());

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
		out << endl;

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

		out << "$" << endl;
	}

	out << std::endl << "# Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << std::endl;
	cout << std::endl << "Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << std::endl;

	if (ncirc) {
		out << "#\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl;
		std::cout << "\tout of which " << ncirc << " are circular of total length " << lcirc << std::endl;
	}

	out.close();
}

// fill in OP container with endpoints of q-obverse paths,
// starting and ending at OP
void get_obverse_paths(const MBGraph& graph, map< vertex_t, set<vertex_t> >& OP, const Mcolor Q) {
    map< vertex_t, set<int> > processed;

    for(auto iq = Q.cbegin(); iq != Q.cend(); ++iq) {
        const partgraph_t& PG = graph.get_local_graph(iq->first);

	for(auto ip = OP.begin(); ip != OP.end(); ++ip) {

            const vertex_t& x = ip->first;

	    if( x==Infty || member(processed[x], iq->first) ) continue;

	    for(vertex_t y = graph.get_adj_vertex(ip->first); PG.defined(y);) {
		if( member(OP,y) ) {
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

		y = graph.get_adj_vertex(y);
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
transform_t decircularize(const MBGraph& graph, partgraph_t& PG, transform_t& TG, const Mcolor& Q) {

    // decircularizing sub-transform that is removed
    transform_t D;

    size_t CircSize = numchr(graph, PG).second;
    if( CircSize == 0 ) return D;

    outlog << "Eliminating " << CircSize << " circular chromosomes in " << genome_match::mcolor_to_name(Q) << endl;

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

    transform_t::iterator start = TG.begin();

    // looking for de-circularizig 2-breaks
    for(auto it = start; it != TG.end();) {
        // check multicolor
	{
	    Mcolor C(it->MultiColor, Q, Mcolor::Intersection);
    
	    if( C.empty() ) {
		++it;
		continue;
	    }
    
	    if( C != Q ) {
		outlog << "Impossible multicolor in the transformation!" << endl;
		break;
	    }
	}

	it->applySingle(T);

        size_t ccsize = numchr(graph, T).second;

	//bool hotfix = ( it->OldArc[0] == arc_t("770h","770t") );

	if( ccsize >= CircSize /* && !hotfix */ ) {
	    ++it;
	    continue;
	}


	//TwoBreak t = *it;
	//t.normalize();

	outlog << "Found problematic 2-break: " << *it << "\t";

	// move t over to beginning
	for(transform_t::iterator jt=it;jt!=TG.begin();) {

	    transform_t::iterator kt = jt--; // jt, kt are successive, *kt == t

	    const TwoBreak& t = *kt;
	    const TwoBreak& s = *jt;
	    //s.normalize();

//            outlog << "... trying to swap with " << s << endl;

	    pair<vertex_t,vertex_t> p1, q1, p2, q2;

	    bool usearc = false;

	    Mcolor C(t.MultiColor, s.MultiColor, Mcolor::Intersection);
	    if( !C.empty() ) {


		/*
			 p1=(x1,x2) x (y1,y2)=q1
			 p2=(x1,y1) x (x3,y3)=q2
    
			 into:
    
			 (x1,x2) x (x3,y3)
			 (y3,x2) x (y1,y2)
		*/
    
		for(int j=0;j<2;++j) {
    
		    if( t.OldArc[j] == make_pair(jt->OldArc[0].first,jt->OldArc[1].first) ) {
			usearc = true;
    
			p2 = t.OldArc[j];
			q2 = t.OldArc[1-j];
    
			p1 = jt->OldArc[0];
			q1 = jt->OldArc[1];
		    }
		    else if( t.OldArc[j] == make_pair(jt->OldArc[1].first,jt->OldArc[0].first) ) {
			usearc = true;
    
			p2 = t.OldArc[j];
			q2 = t.OldArc[1-j];
    
			p1 = jt->OldArc[1];
			q1 = jt->OldArc[0];
		    }
		    else if( t.OldArc[j] == make_pair(jt->OldArc[0].second,jt->OldArc[1].second) ) {
			usearc = true;
    
			p2 = t.OldArc[j];
			q2 = t.OldArc[1-j];
    
			p1 = make_pair(jt->OldArc[0].second, jt->OldArc[0].first);
			q1 = make_pair(jt->OldArc[1].second, jt->OldArc[1].first);
		    }
		    else if( t.OldArc[j] == make_pair(jt->OldArc[1].second,jt->OldArc[0].second) ) {
			usearc = true;
    
			p2 = t.OldArc[j];
			q2 = t.OldArc[1-j];
    
			p1 = make_pair(jt->OldArc[1].second, jt->OldArc[1].first);
			q1 = make_pair(jt->OldArc[0].second, jt->OldArc[0].first);
		    }
		    if(usearc) break;
		}
	    }

	    // TwoBreak t0 = t;

	    if( usearc ) {
		if( t.MultiColor != s.MultiColor ) break;
		*kt = TwoBreak(q2.second,p1.second,q1.first,q1.second,t.MultiColor);
		*jt = TwoBreak(p1.first,p1.second,q2.first,q2.second,t.MultiColor);
	    }
	    else {
		TwoBreak temp = *kt;
		*kt = *jt;
                *jt = temp;
	    }

	    {
		Mcolor C(kt->MultiColor, Q, Mcolor::Intersection);
    
                // N.B. at this point if C is not empty, then C == Q
		if( !C.empty() ) {
		    kt->revertSingle(T);

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
	    outlog << " SUCCEDED" << endl;

	    // move t away from the transformation TG and save it to D
            TG.begin()->applySingle(PG);
	    D.push_back(*TG.begin());

	    TG.erase(TG.begin());

	    CircSize = numchr(graph, PG).second;

	    if( CircSize==0 ) break;

	    start = TG.begin();
	}
	else {  // did not succeed
	    start++;
	    outlog << " FAILED" << endl;
	}

	T = PG;
	for(it = TG.begin();it!=start;++it) {
	    it->applySingle(T);
	}
    }
    //if( start == TG.end() ) {
    if( CircSize>0 ) {
	outlog << "Unable to de-circularize ;-(" << endl;
    }

    return D;
}

/*
Метод собирает статистику после шага. 
Обновляет complement цвета в графе. 
Распечатывает эту статистику в файл stat.txt
*/
void save_information(writer::Wstats& wstats, size_t stage, const ProblemInstance& cfg, MBGraph& graph, bool flag = false) { 
  Statistics st(graph); 
  graph.update_complement_color(st.get_new_color());
  auto p = st.get_compl_stat(graph);
  wstats.print_all_statistics(stage, st, cfg, graph);
  writer::Wdots dots;
  dots.save_dot(graph, cfg, stage);
} 

void main_algorithm(const ProblemInstance& cfg, MBGraph& graph) {
  writer::Wstats write_stats("stats.txt");
  save_information(write_stats, 0, cfg, graph, false);
	
  std::array<bool, 5> print_dots;
  print_dots.fill(true);
  bool process_compl = true; 
  bool isChanged = true; 

  while(isChanged) {
    isChanged = false; 

    if ((cfg.get_stages() >= 1) && !isChanged) {
      Stage1 st1(graph);
      isChanged = st1.stage1();
      graph = st1.get_graph();	

      //if (!isChanged) { isChanged = ConvPhylTreeAll(MBG,26); }  // cut hanging free ends  

      if (print_dots[1] && !isChanged) {
	print_dots[1] = false;    	
	save_information(write_stats, 1, cfg, graph);		
      }

    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      Stage2 st2;//(canformQoo);
      isChanged = st2.stage2(graph);

      if (print_dots[2] && !isChanged) {
	print_dots[2] = false;
	save_information(write_stats, 2, cfg, graph);		    
      }
    }

    if ((cfg.get_stages() >= 3) && !isChanged) {     // STAGE 3, somewhat unreliable
      outlog << "Stage: 3" << std::endl;

      isChanged = stage3(graph);
    
      if (canformQoo && !isChanged) {
	isChanged = true;
	canformQoo = false; // more flexible
      }    

      if (print_dots[3] && !isChanged) {
	print_dots[3] = false;
	save_information(write_stats, 3, cfg, graph, true);
      }
    }

    if ((cfg.get_stages() >= 4) && !isChanged) {
      outlog << "Stage: 4" << std::endl;

      isChanged = stage4(graph);

      if (print_dots[4] && !isChanged) {
	print_dots[4] = false;
	save_information(write_stats, 4, cfg, graph, true);
      }
    }

#ifndef VERSION2
    if (process_compl && !cfg.get_completion().empty() && !isChanged) {     
      outlog << "Manual Completion Stage" << std::endl;

      auto completion = cfg.get_completion();
      for(auto il = completion.begin(); il != completion.end(); ++il) {
	TwoBreak t;
	t.OldArc[0].first = (*il)[0];
	t.OldArc[0].second = (*il)[1];
	t.OldArc[1].first = (*il)[2];
	t.OldArc[1].second = (*il)[3];
	t.MultiColor = genome_match::name_to_mcolor((*il)[4]);

	t.apply(graph, true);
      }
      process_compl = false;
      isChanged = true;
    } 
#endif

  }	

  writer::Wdots dot; 
  dot.save_dot(graph, cfg, 99);

#ifndef VERSION2
  write_stats.print_fair_edges(graph);
  write_stats.histStat();
#else 
  writer::Wdots components;
  components.save_components(graph, cfg, 5);
#endif
} 


////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
  std::cout << "MGRA (Multiple Genome Rearrangements & Ancestors) ver. 1.1" << std::endl;
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
  ProblemInstance PI(reader::read_cfg_file(name_cfg_file)); 

  std::vector<Genome> genomes = reader::read_genomes(PI);
  genome_match::init_name_genomes(genomes);

  MBGraph graph; 
  graph.init(genomes, PI); //create constructor and not global variable

  main_algorithm(PI, graph);

#ifndef VERSION2  
  if (!PI.get_target().empty()) {

#if 0
	// EXPERIMENTAL DECIRCULARIZATION of PI.target
	transform_t H;
	for(transform_t::const_iterator ih=TwoBreak::History.begin();ih!=TwoBreak::History.end();++ih) {
	    H.push_front(ih->inverse());
	}
	for(int i=0;i<N;++i) {
	    transform_t T = decircularize(graph, graph.LG[i],H,TColor[i]); // assume that TColor[i] = { i }
    
	    // move to adjacent branches
	    for(transform_t::const_iterator it = T.begin(); it!=T.end(); ++it) {
		for(int j=0;j<N;++j) if( j!=i && member(it->MultiColor,j) ) {
                    it->applySingle(graph.LG[j]);
		}
	    }
	}
#endif

	partgraph_t PG;

	//ofstream cf("st1comp.res");
	for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {
	    const string& x = *is;
	    if( PG.defined(x) ) continue;

	    string y = Infty;
	    bool good = true;
	    int def = 0;
		std::string target = PI.get_target();
	    for(int i = 0; i < target.size(); ++i) {
		int j = genome_match::get_number(target.substr(i, 1));//PI.get_number_genome_to_name(target.substr(i, 1));
		if (graph.is_there_edge(j, x)) {
		    def++;
		    if( y==Infty ) y = graph.get_adj_vertex(j, x);
		    if( y != graph.get_adj_vertex(j, x) ) good = false;
		}
	    }
	    if( good && def == target.size() && y!=Infty ) {
                PG.insert(x,y);
		//cf << x << "\t" << y << endl;
	    }
	}
	//cf.close();

	set< pair<path_t,bool> > GN;
	splitchr(graph, PG, GN);
	printchr(PI.get_target(), GN, PI.get_target().empty());

	//for(int i=0;i<N;++i) {
	//    set< pair<path_t,bool> > GN;
	//    splitchr(graph.LG[i], GN);
	//    printchr(sname[i].substr(0,1) + "_m",GN, PI.get_target().empty());
	//}


    } else {  /* empty target */

	const size_t NC = graph.DiColor.size();
    
	RG.resize(NC);
	RT.resize(NC);
    
	if( !RecoverGenomes(graph, TwoBreak::History) ) exit(1);
    
	// T-transformation complete, we procede with recovering the ancestral genomes
    
	outlog << "Initial 2-break distances from the root X: " << std::endl;
	for(int i = 0; i < NC; ++i) {
	    outlog << genome_match::mcolor_to_name(graph.TColor[i]) << ":\t" << RT[i].size() << std::endl;
	}
    
	// FIXME: check that the order in which circular chromosomes are eliminated
    
	for(int i = 0; i < NC; ++i) {
    
	    transform_t T = decircularize(graph, RG[i], RT[i], graph.TColor[i]);
    
	    // move to adjacent branches
	    for(transform_t::const_iterator it = T.begin(); it!=T.end(); ++it) {
		for(int j=0;j<NC;++j) {
		    if( j!=i && includes( graph.TColor[i].begin(), graph.TColor[i].end(), graph.TColor[j].begin(), graph.TColor[j].end() ) 
			&& graph.AreAdjacentBranches(graph.TColor[i],graph.TColor[j]) ) {
			RT[j].push_back(*it);
		    }
		}
	    }
	}
    
	outlog << "Final 2-break distances from the root X: " << endl;
	for(int i = 0; i < NC; ++i) {
	    outlog << genome_match::mcolor_to_name(graph.TColor[i]) << ":\t" << RT[i].size() << endl;
	}
    
	for(int i = 0;i < NC; ++i) {
    	    std::set<std::pair<path_t, bool> > GN;
	    splitchr(graph, RG[i], GN);
	    printchr(genome_match::mcolor_to_name(graph.TColor[i]),GN, PI.get_target().empty());
    
	    //splitchr(RG[i], GN, true);
	    //printchr(genome_match::mcolor_to_name(graph.TColor[i]) + "_x",GN, PI.get_target().empty());
    
		std::ofstream tr( (genome_match::mcolor_to_name(graph.TColor[i]) + ".trs").c_str() );
		for(transform_t::const_iterator it=RT[i].begin();it!=RT[i].end();++it) {
			tr << it->OldArc[0].first << " " << it->OldArc[0].second << "\t" << it->OldArc[1].first << " " << it->OldArc[1].second << "\t" << genome_match::mcolor_to_name(it->MultiColor) << endl;
		}
	    tr.close();
	}

    }
#endif
    return 0;
}



///////////////////////////////////////////////////////////////////////////
bool RecoverGenomes(MBGraph& graph, const transform_t& tr) {

    /*
    for(int i=0;i<N;++i) {
      for( partgraph_t::const_iterator il=LG[i].begin();il!=LG[i].end();++il) {
	 if( il->first > il->second) continue;
         cout << il->first << "-" << il->second << " ";
      }
      cout << endl;
    }
    */

    size_t NC = graph.DiColor.size();

    for(int i=0; i < graph.size_graph() - 1; ++i) {
	if( graph.get_local_graph(i) != graph.get_local_graph(i + 1)) {//FIXME
	    std::cout << "T-transformation is not complete. Cannot reconstruct genomes." << std::endl;
	    return false;
	}
    }

    RG.clear(); RG.resize(NC);
    RT.clear(); RT.resize(NC);

    for(int i=0; i < NC; ++i) { 
	RG[i] = graph.get_local_graph(0);
    } 


    // track changes in the number of chromosomes

    // number of reversals, interchromosomal translocations, and fissions/fusions
    std::vector< std::vector<size_t> > RTF(NC);
    for(size_t j=0; j < NC; ++j) RTF[j].resize(3);

    for(auto it = tr.rbegin(); it != tr.rend(); ++it) {
	const Mcolor& Q = it->MultiColor;

	outlog << "Reverting (" << it->OldArc[0].first << "," << it->OldArc[0].second << ")x(" << it->OldArc[1].first << "," << it->OldArc[1].second << "):{" << genome_match::mcolor_to_name(it->MultiColor) << "} " << " in";

	for(size_t i = 0; i < NC; ++i) {

	    if (!Q.includes(graph.TColor[i])) { 
		continue;
	    } 

            // TColor[i] is subset of Q

            size_t nchr_old = 0;
	    if (Q == graph.TColor[i]) {
		nchr_old = numchr(graph, RG[i]).first;
	    }


	    it->revertSingle(RG[i]);

	    if( Q==graph.TColor[i] ) {
		outlog << " " << genome_match::mcolor_to_name(graph.TColor[i]);
		RT[i].push_front(*it);
	    }

	    if( Q == graph.TColor[i] ) {

		bool samechr = true;

		set< string > Vert;
		if( it->OldArc[0].first != Infty ) Vert.insert( it->OldArc[0].first );
		if( it->OldArc[0].second != Infty ) Vert.insert( it->OldArc[0].second );
		if( it->OldArc[1].first != Infty ) Vert.insert( it->OldArc[1].first );
		if( it->OldArc[1].second != Infty ) Vert.insert( it->OldArc[1].second );
    
		getchr(graph, RG[i],*Vert.begin());
    
		for(set<string>::const_iterator iv=++Vert.begin();iv!=Vert.end();++iv) {
		    if( !member(getchrset,*iv) ) {
			samechr = false;
			break;
		    }
		}

		size_t nchr_new = numchr(graph, RG[i]).first;
		if( nchr_new != nchr_old ) {
		    ++RTF[i][2];
		}
		else {
		    if( samechr ) ++RTF[i][0];
		    else ++RTF[i][1];
		}
	    }

	}
	
	outlog << endl;
    }

    vector< size_t > tot(3);
    outlog << "% Number of reversals / translocations / fissions+fusions: " << endl;
    for(size_t j = 0; j < NC; ++j) {
	outlog << genome_match::mcolor_to_name(graph.TColor[j]) << "+" << genome_match::mcolor_to_name(graph.CColor(graph.TColor[j])) << "\t&\t" << RTF[j][0] << " & " << RTF[j][1] << " & " << RTF[j][2]
	    << " &\t" << RTF[j][0]+RTF[j][1]+RTF[j][2] << " \\\\" << endl;
        outlog << "\\hline" << endl;
	tot[0] += RTF[j][0];
	tot[1] += RTF[j][1];
	tot[2] += RTF[j][2];
    }
    outlog << "Total\t&\t" << tot[0] << " & " << tot[1] << " & " << tot[2] << " &\t" << tot[0]+tot[1]+tot[2] << " \\\\" << endl;

    return true;
}

// Given a set of partial BP-graphs (genomes), convolve them
// IG gives a list of genomes for which we need to reconstruct a common ancestor
// return true if LG was modified
bool ConvPhylTreeAll(MBGraph& graph, int stage) {
    bool simplified = false;
    size_t nr = 0; // number of rearrangements, 
    size_t nf = 0; // number of fussions/fissions

    // simplifying graph
    outlog << "Stage: " << stage << endl;

    do {
	nr = 0; 
	nf = 0;

	for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {  
	     const std::string& x = *is;
	     mularcs_t M = graph.get_adjacent_multiedges(x);
	    // generalized reliable simple path
		    if( stage==22 && M.size()>=2 ) {
			for(auto im = M.begin();im!=M.end();++im) {
			  const std::string& y = im->first;
   		          const Mcolor& Q = im->second;
	
	                  if( y==Infty ) continue;
	                    
			  mularcs_t My = graph.get_adjacent_multiedges(y);
			  My.erase(x);
     			  mularcs_t Mx = M;
		    Mx.erase(y);

		    if( member(Mx,Infty) && member(My,Infty) && Mx[Infty]==My[Infty] && graph.is_T_consistent_color(Mx[Infty]) ) {
			Mcolor C(Q, Mx[Infty], Mcolor::Union);
			if (!graph.is_T_consistent_color(C)) continue;
			for(auto iq = Mx[Infty].begin(); iq!=Mx[Infty].end(); ++iq) {
			    graph.add_edge(iq->first, x, y);
			}
			outlog << "Stage 22: fusion " << x << " + " << y << std::endl;
			nf++;
			break;
		    }
		    
		    auto Cx = Mx.cend();
		    for(auto jm = Mx.cbegin();jm!=Mx.cend();++jm) {
			if( !graph.is_T_consistent_color(jm->second) ) continue;
			Mcolor C(Q, jm->second, Mcolor::Union);
			if (graph.is_T_consistent_color(C)) {
			    if( Cx!=Mx.end() ) { Cx=Mx.end(); break; }
			    Cx = jm;
			}
	    	    }
	    	    if( Cx == Mx.end() ) continue;
		    
		    auto Cy = My.cend();
		    for(auto jm = My.cbegin(); jm != My.cend(); ++jm) {
			if( !graph.is_T_consistent_color(jm->second) ) continue;
			Mcolor C(Q, jm->second, Mcolor::Union);
			if( graph.is_T_consistent_color(C) ) {
			    if( Cy!=My.end() ) { Cy=My.end(); break; }
			    Cy = jm;
			}
	    	    }
	    	    if( Cy == My.end() ) continue;
	    	    
	    	    if( Cx->second == Cy->second ) {
                        if( TwoBreak(x,Cx->first,y,Cy->first,Cx->second).apply(graph,true) ) nr++;
	    		outlog << "Stage 22: fusion " << x << " + " << y << endl;
	    		break;
	    	    }
		}
	    }

            // cut hanging free ends
	    if (stage == 26 && M.size() == 2 && member(M,Infty) ) {
		const Mcolor& Q1 = M[Infty];
                vertex_t y;
		Mcolor Q2;
		if( M.begin()->first == Infty) {
                    y = M.rbegin()->first;
		    Q2 = M.rbegin()->second;
		}
		else {
                    y = M.begin()->first;
		    Q2 = M.begin()->second;
		}
		if( !member(graph.DiColor,Q1) && member(graph.DiColor,Q2) /* && !member(graph.get_adjacent_multiedges(y),Infty) */ ) {
		//if( member(DiColor,Q2) && !member(graph.get_adjacent_multiedges(y),Infty) ) {
		    outlog << "Unhanging fission:" << endl;
		    if( TwoBreak(x,y,Infty,Infty,Q2).apply(graph,true) ) nf++;
		}
	    }
	}
	outlog << nr << " 2-breaks performed" << std::endl;
	outlog << nf << " untangling 2-breaks performed" << std::endl;
	if (nr > 0 || nf > 0) simplified = true;
    } while (nr > 0 || nf > 0);
    outlog << std::endl;

    return simplified;
}

