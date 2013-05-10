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

bool RecoverGenomes(const transform_t&);
set <vertex_t> getchrset;

std::pair<path_t, bool> getchr(const partgraph_t& PG, const std::string& x) {
    path_t path;
    bool circular = false;

    getchrset.clear();
    getchrset.insert(x);

    for(vertex_t y = MBG.get_adj_vertex(x); ; ) {
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

	y = MBG.get_adj_vertex(y);
    }

    if( !circular && PG.defined(x) ) {
	for(string y = x;PG.defined(y);) {
	    y = PG[y];
	    getchrset.insert(y);

	    y = MBG.get_adj_vertex(y);
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
void splitchr(const partgraph_t& PG, set< pair<path_t,bool> >& AllChr, const bool Xonly = false, list< set<vertex_t> >& CircChr = pg_empty) {

    if( &CircChr != &pg_empty ) CircChr.clear();
    AllChr.clear();

    set<orf_t> processed;

    for(auto is=MBG.begin_vertices();is!=MBG.end_vertices();++is) {

	const string& x = *is;
	
	if( member(processed,x) ) continue;
       
        pair< path_t, bool > pathb = getchr(PG,x);

	AllChr.insert( pathb );

        copy(getchrset.begin(),getchrset.end(),inserter(processed,processed.end()));

	if( pathb.second && (&CircChr != &pg_empty) ) {
	    CircChr.push_back(getchrset);
	}
    }
}

pair<size_t,size_t> numchr(const partgraph_t& PG) {
    set< pair<path_t,bool> > AllChr;
    list< set<vertex_t> > CircChr;
    splitchr(PG,AllChr,false,CircChr);
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
void get_obverse_paths(map< vertex_t, set<vertex_t> >& OP, const Mcolor Q) {
    map< vertex_t, set<int> > processed;

    for(auto iq = Q.cbegin(); iq != Q.cend(); ++iq) {
        const partgraph_t& PG = MBG.LG[iq->first];

	for(auto ip = OP.begin(); ip != OP.end(); ++ip) {

            const vertex_t& x = ip->first;

	    if( x==Infty || member(processed[x], iq->first) ) continue;

	    for(vertex_t y = MBG.get_adj_vertex(ip->first); PG.defined(y);) {
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

		y = MBG.get_adj_vertex(y);
	    }
	    processed[x].insert(iq->first);
	}
    }
}

// Save .dot file and output statistics of synteny blocks representing breakpoints
void save_dot(const ProblemInstance& cfg, size_t stage, bool outbad = false) {
  std::string dotname = cfg.get_graphname() + toString(stage) + ".dot";
  std::ofstream dot(dotname.c_str());

  dot << "graph {" << std::endl;
  if (!cfg.get_colorscheme().empty()) { 
    dot << "edge [colorscheme=" << cfg.get_colorscheme() << "];" << std::endl;
  } 

  int infv = 0;
  size_t ncls = 0;
  std::unordered_set<vertex_t> mark; // vertex set
  for(auto is = MBG.begin_vertices(); is != MBG.end_vertices(); ++is) {
    const std::string& x = *is;

    //multimularcs_t Mx = MBG.get_adjacent_multiedges_v2(x);
    mularcs_t Mx = MBG.get_adjacent_multiedges(x);

    if (Mx.size() == 1) { 
      continue; // trivial cycle
    } 

    for(auto im = Mx.cbegin(); im != Mx.cend(); ++im) {
      const std::string& y = im->first;

      if (mark.find(y) != mark.end()) { 
	continue; // already output
      } 

      const Mcolor& C = im->second;
      for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
       	/* ************** output edge (x,y) **************** */
#ifndef VERSION2
	dot << "\t\"" << x << "\"\t--\t\"";
	if (y == Infty) {
	  if (ic == C.cbegin()) { 
	    --infv;
	  } 
	  dot << infv << "\"\t[len=0.75,";
	} else { 
	  dot << y << "\"\t[";
	} 
	dot << "color=" <<  cfg.get_RGBcolor(cfg.get_RGBcoeff() * (ic->first)) << "];" << std::endl;
#else 
	if (y != Infty) { 
	  dot << "\t\"" << x << "\"\t--\t\"";
	  dot << y << "\"\t[";
	  dot << "color=" <<  cfg.get_RGBcolor(cfg.get_RGBcoeff() * (ic->first)) << "];" << std::endl;
	} 
#endif
      }
    }
    mark.insert(x);
  }

#ifndef VERSION2
  for(int i = infv; i < 0; ++i) {
    dot << "\t\"" << i << "\"\t[shape=point,color=black];" << std::endl;
  }
#endif

  dot << "}" << std::endl;
  dot.close();
}

void save_components(const ProblemInstance& cfg, size_t stage) { 
  std::string dotname = cfg.get_graphname() + toString(stage);

  equivalence<vertex_t> CC; // connected components
  		
  for(auto is = MBG.begin_vertices(); is != MBG.end_vertices(); ++is) {
    CC.addrel(*is, *is);
  } 
  
  for(size_t i = 0; i < MBG.size_graph(); ++i) {
    for(auto il = MBG.LG[i].cbegin(); il != MBG.LG[i].cend(); ++il) {
      CC.addrel(il->first, il->second);
    }
  }

  CC.update();    
  /*for(auto is = MBG.begin_vertices(); is != MBG.end_vertices(); ++is) {    
    mularcs_t M = MBG.get_adjacent_multiedges(*is);

    for(auto im = M.cbegin(); im != M.cend(); ++im) {    
	CC.addrel(*is, im->first);
    }
  }
  
  CC.update();*/
  
  std::map<std::string, std::set<std::string> > components;
  CC.get_eclasses(components);
    
  std::cerr << components.size() << std::endl;

  size_t i = 0; 
  for(auto it = components.cbegin(); it != components.cend(); ++it) { 
    const std::set<std::string>& current = it->second;
    if (current.size() <= 7) {
      continue;
    }

    std::cerr << current.size() << std::endl; 
  
    std::string namefile = dotname + "_" + toString(++i) + ".dot"; 
    std::ofstream dot(namefile.c_str());

    dot << "graph {" << std::endl;
    if (!cfg.get_colorscheme().empty()) { 
      dot << "edge [colorscheme=" << cfg.get_colorscheme() << "];" << std::endl;
    } 

    std::unordered_set<vertex_t> mark; // vertex set
    for(auto is = current.cbegin(); is != current.cend(); ++is) {
      const std::string& x = *is;
  
      if (x == Infty) { 
	continue;
      } 

      mularcs_t Mx = MBG.get_adjacent_multiedges(x);

      if (Mx.size() == 1) { 
	continue; // trivial cycle
      } 
      
      for(auto im = Mx.cbegin(); im != Mx.cend(); ++im) {
	const std::string& y = im->first;
	
	if (mark.find(y) != mark.end()) { 
	  continue; // already output
	} 

	const Mcolor& C = im->second;
	for(auto ic = C.cbegin(); ic != C.cend(); ++ic) {
	  if (y != Infty) { 
	    dot << "\t\"" << x << "\"\t--\t\"";
	    dot << y << "\"\t[";
	    dot << "color=" <<  cfg.get_RGBcolor(cfg.get_RGBcoeff() * (ic->first)) << "];" << std::endl;
	  } 
	}
      }
      mark.insert(x);
    }

    dot << "}" << std::endl;
    dot.close();
  } 

} 

/*
canformQ(x,Q) говорит, можно ли из мультицветов мультиребер инцидентных вершине x образовать мультицвет Q.
*/
/* 
can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
if return false, then Q cannot be formed
if true - who knows...
*/
bool canformQoo = true; // safe choice, at later stages may change to false

bool canformQ(const std::string& x, const Mcolor& Q) {
    if (x == Infty) {
	return canformQoo;
    }

    // color Q can be formed if some it adjacent multicolors form a partition of Q
    // OR 
    // if every intersection Q \cap QQ = \emptyset or QQ.

    multimularcs_t M = MBG.get_adjacent_multiedges_with_split(x);
    for(auto im = M.cbegin(); im != M.cend(); ++im) { 
	Mcolor C(Q, im->second, Mcolor::Intersection); 
	if (C.size() > 0 && C.size() < im->second.size()) { 
		return false;
	} 
    }
    return true;
}


/* Given a non-linear genome PG of a multicolor Q and a transformation into a linear genome,
 * find linearizing fissions in the transformation, move them to beginning, and apply to PG
 * i.e., try to reorder transformation into: PG -> PG' -> linear genome, where PG' is linear
 * and obtained from PG with fission.
 * We replace PG with PG' and return the transformation PG -> PG'
 * Transformation may contain only multicolors Q' with Q'\cap Q = 0 or Q.
*/
transform_t decircularize(partgraph_t& PG, transform_t& TG, const Mcolor& Q) {

    // decircularizing sub-transform that is removed
    transform_t D;

    size_t CircSize = numchr(PG).second;
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

        size_t ccsize = numchr(T).second;

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

		    ccsize = numchr(T).second;
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

	    CircSize = numchr(PG).second;

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
  MBG.update_complement_color(st.get_new_color());
  auto p = st.get_compl_stat(graph);
  wstats.print_all_statistics(stage, st, cfg, graph);
  save_dot(cfg, stage, flag);
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
      isChanged = stage1(graph);
      //if (!isChanged) { isChanged = ConvPhylTreeAll(MBG,26); }  // cut hanging free ends  

      if (print_dots[1] && !isChanged) {
	print_dots[1] = false;    		
	save_information(write_stats, 1, cfg, graph);		
      }
    }

    if ((cfg.get_stages() >= 2) && !isChanged) {
      isChanged = ConvPhylTreeAll(graph, 2); 

      if (print_dots[2] && !isChanged) {
	print_dots[2] = false;
	save_information(write_stats, 2, cfg, graph);		    
      }
    }

    if ((cfg.get_stages() >= 3) && !isChanged) {     // STAGE 3, somewhat unreliable
      outlog << "Stage: 3" << std::endl;
      isChanged = ConvPhylTreeAll(graph, 222); // cut the graph into connected components

      if (!isChanged) { 
	isChanged = ConvPhylTreeAll(graph, 2222); // process 4-cycles
      } 
    
      if (canformQoo && !isChanged) {
	isChanged = true;
	canformQoo = false; // more flexible    FIXME
      }    

      if (print_dots[3] && !isChanged) {
	print_dots[3] = false;
	save_information(write_stats, 3, cfg, graph, true);
      }
    }

    if ((cfg.get_stages() >= 4) && !isChanged) {
      outlog << "Stage: 4" << std::endl;

      outlog << "SplitBadColors is ON" << std::endl;
      graph.SplitBadColors = true; //FIXME

      isChanged = ConvPhylTreeAll(graph, 2);

      graph.SplitBadColors = false;
      outlog << "SplitBadColors is OFF" << std::endl;

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


  save_dot(cfg, 99, true);

#ifndef VERSION2
  write_stats.print_fair_edges(graph);
  write_stats.histStat();
#else 
  save_components(cfg, 5);
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

  MBG.init(genomes, PI); //create constructor and not global variable

  main_algorithm(PI, MBG);

#ifndef VERSION2  
  if (!PI.get_target().empty()) {

#if 0
	// EXPERIMENTAL DECIRCULARIZATION of PI.target
	transform_t H;
	for(transform_t::const_iterator ih=TwoBreak::History.begin();ih!=TwoBreak::History.end();++ih) {
	    H.push_front(ih->inverse());
	}
	for(int i=0;i<N;++i) {
	    transform_t T = decircularize(MBG.LG[i],H,TColor[i]); // assume that TColor[i] = { i }
    
	    // move to adjacent branches
	    for(transform_t::const_iterator it = T.begin(); it!=T.end(); ++it) {
		for(int j=0;j<N;++j) if( j!=i && member(it->MultiColor,j) ) {
                    it->applySingle(MBG.LG[j]);
		}
	    }
	}
#endif

	partgraph_t PG;

	//ofstream cf("st1comp.res");
	for(auto is=MBG.begin_vertices();is!=MBG.end_vertices();++is) {
	    const string& x = *is;
	    if( PG.defined(x) ) continue;

	    string y = Infty;
	    bool good = true;
	    int def = 0;
		std::string target = PI.get_target();
	    for(int i = 0; i < target.size(); ++i) {
		int j = genome_match::get_number(target.substr(i, 1));//PI.get_number_genome_to_name(target.substr(i, 1));
		if( MBG.LG[j].defined(x) ) {
		    def++;
		    if( y==Infty ) y=MBG.LG[j][x];
		    if( y!=MBG.LG[j][x] ) good = false;
		}
	    }
	    if( good && def == target.size() && y!=Infty ) {
                PG.insert(x,y);
		//cf << x << "\t" << y << endl;
	    }
	}
	//cf.close();

	set< pair<path_t,bool> > GN;
	splitchr(PG, GN);
	printchr(PI.get_target(), GN, PI.get_target().empty());

	//for(int i=0;i<N;++i) {
	//    set< pair<path_t,bool> > GN;
	//    splitchr(MBG.LG[i], GN);
	//    printchr(sname[i].substr(0,1) + "_m",GN, PI.get_target().empty());
	//}


    } else {  /* empty target */

	const size_t NC = MBG.DiColor.size();
    
	RG.resize(NC);
	RT.resize(NC);
    
	if( !RecoverGenomes(TwoBreak::History) ) exit(1);
    
	// T-transformation complete, we procede with recovering the ancestral genomes
    
	outlog << "Initial 2-break distances from the root X: " << std::endl;
	for(int i = 0; i < NC; ++i) {
	    outlog << genome_match::mcolor_to_name(MBG.TColor[i]) << ":\t" << RT[i].size() << std::endl;
	}
    
	// FIXME: check that the order in which circular chromosomes are eliminated
    
	for(int i = 0; i < NC; ++i) {
    
	    transform_t T = decircularize(RG[i],RT[i],MBG.TColor[i]);
    
	    // move to adjacent branches
	    for(transform_t::const_iterator it = T.begin(); it!=T.end(); ++it) {
		for(int j=0;j<NC;++j) {
		    if( j!=i && includes( MBG.TColor[i].begin(), MBG.TColor[i].end(), MBG.TColor[j].begin(), MBG.TColor[j].end() ) 
			&& MBG.AreAdjacentBranches(MBG.TColor[i],MBG.TColor[j]) ) {
			RT[j].push_back(*it);
		    }
		}
	    }
	}
    
	outlog << "Final 2-break distances from the root X: " << endl;
	for(int i = 0; i < NC; ++i) {
	    outlog << genome_match::mcolor_to_name(MBG.TColor[i]) << ":\t" << RT[i].size() << endl;
	}
    
	for(int i = 0;i < NC; ++i) {
    	    std::set<std::pair<path_t, bool> > GN;
	    splitchr(RG[i], GN);
	    printchr(genome_match::mcolor_to_name(MBG.TColor[i]),GN, PI.get_target().empty());
    
	    //splitchr(RG[i], GN, true);
	    //printchr(genome_match::mcolor_to_name(MBG.TColor[i]) + "_x",GN, PI.get_target().empty());
    
		std::ofstream tr( (genome_match::mcolor_to_name(MBG.TColor[i]) + ".trs").c_str() );
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
bool RecoverGenomes(const transform_t& tr) {

    /*
    for(int i=0;i<N;++i) {
      for( partgraph_t::const_iterator il=LG[i].begin();il!=LG[i].end();++il) {
	 if( il->first > il->second) continue;
         cout << il->first << "-" << il->second << " ";
      }
      cout << endl;
    }
    */

    size_t NC = MBG.DiColor.size();

    for(int i=0;i<MBG.size_graph()-1;++i) {
	if( MBG.LG[i]!=MBG.LG[i+1] ) {
	    cout << "T-transformation is not complete. Cannot reconstruct genomes." << endl;
	    return false;
	}
    }

    RG.clear(); RG.resize(NC);
    RT.clear(); RT.resize(NC);

    for(int i=0;i<NC;++i) RG[i]=MBG.LG[0];


    // track changes in the number of chromosomes

    // number of reversals, interchromosomal translocations, and fissions/fusions
    std::vector< std::vector<size_t> > RTF(NC);
    for(size_t j=0; j < NC; ++j) RTF[j].resize(3);

    for(auto it = tr.rbegin(); it != tr.rend(); ++it) {
	const Mcolor& Q = it->MultiColor;

	outlog << "Reverting (" << it->OldArc[0].first << "," << it->OldArc[0].second << ")x(" << it->OldArc[1].first << "," << it->OldArc[1].second << "):{" << genome_match::mcolor_to_name(it->MultiColor) << "} " << " in";

	for(size_t i = 0; i < NC; ++i) {

	    if (!Q.includes(MBG.TColor[i])) { 
		continue;
	    } 

            // TColor[i] is subset of Q

            size_t nchr_old = 0;
	    if (Q == MBG.TColor[i]) {
		nchr_old = numchr(RG[i]).first;
	    }


	    it->revertSingle(RG[i]);

	    if( Q==MBG.TColor[i] ) {
		outlog << " " << genome_match::mcolor_to_name(MBG.TColor[i]);
		RT[i].push_front(*it);
	    }

	    if( Q == MBG.TColor[i] ) {

		bool samechr = true;

		set< string > Vert;
		if( it->OldArc[0].first != Infty ) Vert.insert( it->OldArc[0].first );
		if( it->OldArc[0].second != Infty ) Vert.insert( it->OldArc[0].second );
		if( it->OldArc[1].first != Infty ) Vert.insert( it->OldArc[1].first );
		if( it->OldArc[1].second != Infty ) Vert.insert( it->OldArc[1].second );
    
		getchr(RG[i],*Vert.begin());
    
		for(set<string>::const_iterator iv=++Vert.begin();iv!=Vert.end();++iv) {
		    if( !member(getchrset,*iv) ) {
			samechr = false;
			break;
		    }
		}

		size_t nchr_new = numchr(RG[i]).first;
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
	outlog << genome_match::mcolor_to_name(MBG.TColor[j]) << "+" << genome_match::mcolor_to_name(MBG.CColor(MBG.TColor[j])) << "\t&\t" << RTF[j][0] << " & " << RTF[j][1] << " & " << RTF[j][2]
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
bool ConvPhylTreeAll(MBGraph& MBG, int stage) {
    bool simplified = false;
    size_t nr = 0; // number of rearrangements, 
    size_t nf = 0; // number of fussions/fissions

    // simplifying graph
    outlog << "Stage: " << stage << endl;

    do {
	nr = 0; 
	nf = 0;

	// cut the graph into connected components
	if (stage == 222) {

	    outlog << "Stage 222: splitting into connected components" << std::endl;

            // go over all T-consistent multicolors
	    for(auto ic = MBG.DiColor.begin(); ic != MBG.DiColor.end(); ++ic) {
		if (ic->size() == 0 || ic->size() == MBG.size_graph()) { 
			continue; // except empty and complete multicolor
		} 		
		const Mcolor& Q = *ic;

		bool repeat = true;
		while(repeat) {
		    repeat = false;

		    equivalence<vertex_t> CC; // connected components
		    std::map<vertex_t, vertex_t> QQ; // multiedges of colors !Q (!*ic)
		    
		    for(auto is = MBG.begin_vertices(); is != MBG.end_vertices(); ++is) {    
			mularcs_t M = MBG.get_adjacent_multiedges(*is);

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
				for(size_t i = 0; i < MBG.size_graph(); ++i) {
					if( MBG.LG[i].defined(p.first) && MBG.LG[i][p.first]==t.first ) {
				   	 T.insert(i);
					}
				} 
                                // if e is enriched to T-consistent color, great!
				if( T.size()>Q.size() && MBG.is_T_consistent_color(T) ) {
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
    
    
			    if( ! TwoBreak(p, q, Q).apply(MBG, true) ) continue;
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

			    if( ! TwoBreak(p,q,Q).apply(MBG,true) ) continue;
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
	}

	if( stage == 2222 ) {

            // search for 4-cycles

	    for(auto is = MBG.begin_vertices(); is != MBG.end_vertices(); ++is) {
		const string& x = *is;
		mularcs_t Mx = MBG.get_adjacent_multiedges(x);

                bool next = false;

		for(auto im = Mx.begin(); im!=Mx.end(); ++im) {

		    const std::string& y = im->first;
		    const Mcolor& Q = im->second;

		    if (!member(MBG.DiColor, Q) || y==Infty) { 
			continue;
		    }

		    mularcs_t My = MBG.get_adjacent_multiedges(y);
		    My.erase(x);

		    for(auto jm = My.begin(); jm != My.end(); ++jm) {
			std::string z = jm->first;
			if (z == Infty) { 
				continue;
			} 
	
			mularcs_t Cz = MBG.get_adjacent_multiedges(z);

			vertex_t v = "";
			for (auto ir = Cz.cbegin(); ir != Cz.cend(); ++ir) {
				if (ir->second == Q) { 
					v = ir->first;
				} 
			}  

			if ((!v.empty()) && member(Mx, v) ) {

			    auto p = std::make_pair(x, y);
			    auto q = std::make_pair(v, z);

			    outlog << "Stage 2222: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << std::endl;

			    if (TwoBreak(p, q, Q).apply(MBG, true)) nr++;

                            next = true;

			    break;

			}
		    }
                    if(next) break;
		}
	    }
	}


   for(auto is = MBG.begin_vertices(); is!=MBG.end_vertices(); ++is) {  
     const std::string& x = *is;
     mularcs_t M = MBG.get_adjacent_multiedges(x);
    // generalized reliable simple path
	    if( stage==22 && M.size()>=2 ) {
		for(auto im = M.begin();im!=M.end();++im) {
		  const std::string& y = im->first;
		    const Mcolor& Q = im->second;

                    if( y==Infty ) continue;
                    
		    mularcs_t My = MBG.get_adjacent_multiedges(y);
		    My.erase(x);
		    mularcs_t Mx = M;
		    Mx.erase(y);

		    if( member(Mx,Infty) && member(My,Infty) && Mx[Infty]==My[Infty] && MBG.is_T_consistent_color(Mx[Infty]) ) {
			Mcolor C(Q, Mx[Infty], Mcolor::Union);
			if (!MBG.is_T_consistent_color(C)) continue;
			for(auto iq = Mx[Infty].begin(); iq!=Mx[Infty].end(); ++iq) {
			    MBG.LG[iq->first].insert(x,y);
			}
			outlog << "Stage 22: fusion " << x << " + " << y << std::endl;
			nf++;
			break;
		    }
		    
		    
		    

		    auto Cx = Mx.cend();
		    for(auto jm = Mx.cbegin();jm!=Mx.cend();++jm) {
			if( !MBG.is_T_consistent_color(jm->second) ) continue;
			Mcolor C(Q, jm->second, Mcolor::Union);
			if (MBG.is_T_consistent_color(C)) {
			    if( Cx!=Mx.end() ) { Cx=Mx.end(); break; }
			    Cx = jm;
			}
	    	    }
	    	    if( Cx == Mx.end() ) continue;
		    
		    auto Cy = My.cend();
		    for(auto jm = My.cbegin(); jm != My.cend(); ++jm) {
			if( !MBG.is_T_consistent_color(jm->second) ) continue;
			Mcolor C(Q, jm->second, Mcolor::Union);
			if( MBG.is_T_consistent_color(C) ) {
			    if( Cy!=My.end() ) { Cy=My.end(); break; }
			    Cy = jm;
			}
	    	    }
	    	    if( Cy == My.end() ) continue;
	    	    
	    	    if( Cx->second == Cy->second ) {
                        if( TwoBreak(x,Cx->first,y,Cy->first,Cx->second).apply(MBG,true) ) nr++;
	    		outlog << "Stage 22: fusion " << x << " + " << y << endl;
	    		break;
	    	    }
		}
	    }

		    


	    if( stage==2 ) {
/*
		bool xgood = true;
		for(map< string, set<int> >::const_iterator im = M.begin();im!=M.end();++im) {
		    if( !member(Color, im->second) ) xgood = false;
		}
		if( !xgood ) continue;
*/

#ifdef VERSION2
if (!MBG.is_fair_vertice(M)) { 
	continue; 
} 
#endif
		for(auto im = M.begin(); im != M.end(); ++im) {
			const std::string& y = im->first;

			const Mcolor& Q = im->second; // color of central edge

			if (y == Infty || Q.size() == MBG.size_graph()) { 
				continue;
			} 

			multimularcs_t Cx = MBG.get_adjacent_multiedges_with_split(x);
			multimularcs_t Cy = MBG.get_adjacent_multiedges_with_split(y);
#ifdef VERSION2
if (!MBG.is_fair_vertice(MBG.get_adjacent_multiedges(y))) { 
	continue;
} 
#endif
//if( MBG.get_adjacent_multiedges(y).size()!=3 ) continue;

        		outlog << "Testing mobility of edge " << x << "-" << y << " " << genome_match::mcolor_to_name(Q) << " ";

			// test "mobility" of central edge
			// can it be ever find neighboring edge of the same multicolor
			bool mobilQ = false;

//here 
			for(auto jc = Cx.cbegin(); jc!= Cx.cend(); ++jc) {
				if (jc->first != y) { continue; } // not a cental sub-edge

				const Mcolor& QQ = jc->second; // color of central sub-edge (QQ is sub-multicolor of Q)
	                        if (!member(MBG.DiColor, QQ)) { continue; } 


				for(auto ix = Cx.begin(); ix != Cx.end(); ++ix) { 
					if (ix->first == y  /*|| ix->first == Infty*/ ) { continue; } 

					if (canformQ(ix->first, QQ)) {
						outlog << "MOBIL: " << x << "-" << ix->first << " canForm: " << genome_match::mcolor_to_name(QQ) << std::endl;
						mobilQ = true;
						break;
					}
				}
	
				if (mobilQ) { 
					break;
				} 
    
				for(auto iy = Cy.cbegin(); iy != Cy.cend(); ++iy) { 
					if (iy->first == x  /*|| iy->first == Infty*/ ) { continue; } 
			    		if (canformQ(iy->first, QQ)) {
						outlog << "MOBIL: " << y << "-" << iy->first << " canForm: " << genome_match::mcolor_to_name(QQ) << std::endl;
						mobilQ = true;
						break;
			    		}
				}

				if (mobilQ) { 
					break;
				} 
		    	}

			if (mobilQ) continue;

			outlog << "NOT MOBIL" << std::endl;

			bool found = false;

			for (auto ix = Cx.cbegin(); ix != Cx.cend(); ++ix) { 
				if (ix->first == y) continue;
				const Mcolor& QQ = ix->second;

				outlog << " Sub-multiedge " << genome_match::mcolor_to_name(ix->second) << std::endl;

				vertex_t temp = "";   
				for(auto iy = Cy.cbegin(); iy != Cy.cend(); ++iy) { 
					if (iy->second == ix->second) { 	
						temp = iy->first;
						break; 
					} 
				} 

				if (!member(MBG.DiColor, QQ) || temp.empty()) continue; 

				/*
				Mcolor C(Q, QQ, Mcolor::Union);
				// TODO: graph on central edges
				//if( !MBG.is_T_consistent_color(C) ) continue; // do not create T-consistent color
	                        */

				if(TwoBreak(x, ix->first, y, temp, QQ).apply(MBG, true)) {
	                        	found = true;
					++nf;
				}
			}

			if (found) { break; } 
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
		if( !member(MBG.DiColor,Q1) && member(MBG.DiColor,Q2) /* && !member(MBG.get_adjacent_multiedges(y),Infty) */ ) {
		//if( member(DiColor,Q2) && !member(MBG.get_adjacent_multiedges(y),Infty) ) {
		    outlog << "Unhanging fission:" << endl;
		    if( TwoBreak(x,y,Infty,Infty,Q2).apply(MBG,true) ) nf++;
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

