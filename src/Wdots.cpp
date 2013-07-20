#include "Wdots.h"

void writer::Wdots::save_dot(const mbgraph_with_history<Mcolor>& graph, const ProblemInstance& cfg, size_t stage) { 
  std::string dotname = cfg.get_graphname() + toString(stage) + ".dot";
  std::ofstream dot(dotname.c_str());

  dot << "graph {" << std::endl;
  if (!cfg.get_colorscheme().empty()) { 
    dot << "edge [colorscheme=" << cfg.get_colorscheme() << "];" << std::endl;
  } 

  int infv = 0;
  size_t ncls = 0;
  std::unordered_set<vertex_t> mark; // vertex set
  for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {
    const std::string& x = *is;

    Mularcs<Mcolor> Mx = graph.get_adjacent_multiedges(x);

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
       	/*************** output edge (x,y) **************** */
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

void writer::Wdots::save_components(const mbgraph_with_history<Mcolor>& graph, const ProblemInstance& cfg, size_t stage) { 
  std::string dotname = cfg.get_graphname() + toString(stage);

  equivalence<vertex_t> CC; // connected components
  		
  for(auto is = graph.begin_vertices(); is != graph.end_vertices(); ++is) {
    CC.addrel(*is, *is);
  } 
  
  for(auto lc = graph.begin_local_graphs(); lc != graph.end_local_graphs(); ++lc) { 
    for(auto il = lc->cbegin(); il != lc->cend(); ++il) {
      CC.addrel(il->first, il->second);
    }
  }

  CC.update();    
  
  std::map<std::string, std::set<std::string> > components;
  CC.get_eclasses(components);
  
  size_t i = 0; 
  for(auto it = components.cbegin(); it != components.cend(); ++it) { 
    const std::set<std::string>& current = it->second;
    if (current.size() <= 7) {
      continue;
    }

    //std::cerr << current.size() << std::endl; 
  
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

      Mularcs<Mcolor> Mx = graph.get_adjacent_multiedges(x);

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

/*void writer::Wdots::write_legend_dot(size_t size, const std::vector<std::string>& output, const ProblemInstance& cfg) { 
    	std::ofstream output("legend.dot");

	output << "digraph Legend {" << std::endl;
	output << "\tnode [style=filled];" << std::endl;

	/*for (size_t j = 0; j < size; ++j) {
	    flegend << "\t\"" << cfg.get_name(j) << "\"\t[fillcolor=" <<  cfg.get_RGBcolor(cfg.get_RGBcoeff() * j)  << "];" << std::endl;
	} 

  	for(auto it = info.cbegin(); it != info.cend(); ++it) {
		output << *it << std::endl;
	} 

	output << "}" << std::endl;
	output.close();
} */


/*writer::Wdots::Wdots(std::string name_file) { 
	output.open(name_file); 
} */

#if 0 
// Save .dot file and output statistics of synteny blocks representing breakpoints
const set<string> EmptySS;
//void save_dot(const string& dotname, const set<string>& VV = EmptySS) {
void writer::Wdots::save_dot(bool outbad = false, const set<string>& VV = EmptySS) { 
    ofstream dot(dotname.c_str());
    dot << "graph {" << endl;
    dot << "edge [colorscheme=" << PI.get_colorscheme() << "];" << endl;
    //dot << "node [fontname=\"freefont/FreeSans\",fontsize=8];" << endl;

#ifdef BPG_COLORED

    map<string,int> gen2chr;
    equivalence<string> SameColor;
    set<string> Biff = graph.VertexSet;
    set<string> noniso; // non-isolated vertices
    
    for(set<string>::const_iterator is=graph.VertexSet.begin();is!=graph.VertexSet.end();++is) {
	const string& x = *is;
	if( !member(Biff,x) ) continue;

	string xx = x;
	xx.resize(xx.size()-1);
	string chr = G[REFGENOME].orf2cpan[xx].first;
	//if(chr=="chrX") chr = "0";
	for(int i=0;i<chr.size();++i) if(chr[i]<'0'||chr[i]>'9') chr[i]=' ';
	istringstream istr(chr); int cn; istr >> cn;
	gen2chr[x] = cn;
        gen2chr[graph.OB[x]] = cn;

#ifndef BPG_COMPRESSED
	Mularcs Mx = graph.mularcs(x);
	if( Mx.size() == 1 ) Biff.erase(x);
    }
#else

	SameColor.addrel(x,graph.OB[x]);

	Mularcs Mx = graph.mularcs(x);
#ifndef UNICOLCOMP
	if( Mx.size() == 1 ) {
	    const string& y = Mx.begin()->first;
#else
	for(Mularcs::const_iterator im=Mx.begin();im!=Mx.end();++im) {
	    const string& y = im->first;
#endif
	    if( y==Infty || (member(gen2chr,y) && gen2chr[y]==cn) ) {
		Biff.erase(x);
		if( y!=Infty ) {
		    Biff.erase(y);
		    SameColor.addrel(x,y);
		}
	    }
	}
    }
    SameColor.update();
#endif

#endif /* BPG_COLORED */




    int infv = 0;
    size_t ncls = 0;
    set< orf_t > V; // vertex set
    for(set<string>::const_iterator is=graph.VertexSet.begin();is!=graph.VertexSet.end();++is) {
	const string& x = *is;

	Mularcs Mx = graph.mularcs(x);
#ifndef FULL_BPG
	if( Mx.size() == 1 ) continue; // trivial cycle
#endif	
	//if( !VV.empty() && !member(VV,x) ) continue;



#ifdef FULL_BPG
        // obverse edges
	if( x[x.size()-1]=='t' ) {
	    string xt = x, xh = x;
	    xh[xh.size()-1]='h';
	    // dot << "\t\"" << x << "\"\t--\t\"" << xh << "\"\t[style=dashed,color=9];" << endl;

#ifndef UNICOLCOMP
#ifdef BPG_COLORED
	    if( member(Biff,xt) || member(Biff,xh) ) {
	        if( !member(Biff,xt) ) {
		    noniso.insert( SameColor[xt] );
		    xt = SameColor[xt] + "CC";
		}
		if( !member(Biff,xh) ) {
		    noniso.insert( SameColor[xh] );
		    xh = SameColor[xh] + "CC";
		}
#endif

		dot << "\t\"" << xt << "\"\t--\t\"" << xh << "\"\t[style=dashed,color=" << graph.size_graph()+1 << "];" << endl;
	    }
            // EQ.addrel(x,xh);
#endif  // UNICOLCOMP

	}
#endif  // FULL_BPG



	//if( Mx.size()==N ) continue;

	for(Mularcs::const_iterator im=Mx.begin();im!=Mx.end();++im) {
	    const string& y = im->first;

	    if( member(V,y) ) continue; // already output

#if 0
	    Mularcs Mx1 = graph.mularcs(x);
	    Mx1.erase(Infty);
	    if( Mx1.size() >= 3 ) continue;
	    if( y!=Infty ) {
		Mularcs My = graph.mularcs(y);
		My.erase(Infty);
		if( My.size() >= 3 ) continue;
		EQ.addrel(x,y);
	    }
#endif


	    const set<int>& C = im->second;
	    for(set<int>::const_iterator ic=C.begin();ic!=C.end();++ic) {
       	// ************** output edge (x,y) ****************

#ifdef BPG_COMPRESSED

#ifndef UNICOLCOMP
		if( member(Biff,x) || member(Biff,y) ) {
#else
		if( y!=Infty && !SameColor.isequiv(x,y) ) {
#endif

		    string u = x;
		    string v = y;
#ifndef UNICOLCOMP
		    if( !member(Biff,u) )
#endif
		    {
			noniso.insert( SameColor[u] );
			u = SameColor[u] + "CC";
		    }
#ifndef UNICOLCOMP
		    if( !member(Biff,v) && v!=Infty )
#endif
		    {
			noniso.insert( SameColor[v] );
			v = SameColor[v] + "CC";
		    }

		    dot << "\t\"" << u << "\"\t--\t\"";
		    if( v==Infty ) {
			if(ic==C.begin()) --infv;
			dot << infv;
		    }
                    else dot << v;
		    dot << "\"\t[color=";  // << 1 + *ic;

#ifdef RESTRICT2MDH
		    switch( *ic ) {
		    case 0: dot << 1; break;
		    case 1: dot << 3; break;
		    case 2: dot << 4; break;
		    }
#else
		    dot << RGBcols[RGBcoeff * (*ic)];
#endif
		    
		    if( Mx.size()==1 ) {
			dot << ",penwidth=3";
		    }
		    dot << "];" << endl;

		}

#else /* ! BPG_COMPRESSED */

		dot << "\t\"" << x << "\"\t--\t\"";
		if( y==Infty ) {
                    if(ic==C.begin()) --infv;
		    dot << infv << "\"\t[len=0.75,";
		}
		else dot << y << "\"\t[";
		//dot << "color=" << color[*ic] << "];" << endl;
		dot << "color=" << RGBcols[RGBcoeff * (*ic)] << "];" << endl;
#endif


	    }
	    //if( y!=Infty && !VV.empty() && !member(VV,y) ) {
	    //    dot << "\t\"" << y << "\"\t[shape=diamond];" << endl;
	    //}
	}

	//if( member(chrX,x) ) {
	//    dot << "\t\"" << x << "\"\t[style=filled,fillcolor=8];" << endl;
	//}
	//else if( Mx.size()>3 ) {
	//    dot << "\t\"" << x << "\"\t[style=filled,fillcolor=9];" << endl;
	//}
	
	V.insert(x);
    }

#ifdef BPG_COMPRESSED
    map< string, set<string> > CC;
    SameColor.get_eclasses(CC);
    // output boxes
    for(map< string, set<string> >::const_iterator ic=CC.begin();ic!=CC.end();++ic) {
    
	int cn = gen2chr[ ic->first ];
    
    //    bool bf = false;
	size_t sz = 0;
	set<string> processed; // to account possible obverses in the same chain
	for(set<string>::const_iterator is=ic->second.begin();is!=ic->second.end();++is) {
	    string gen = *is;
	    gen.resize(gen.size()-1);
	    if( member(processed,gen) ) continue;
	    //pair<int,int> p = G[REFGENOME].orf2cpan[gen].second;
	    //sz += p.second - p.first;
            sz += BL[*is];
	    processed.insert(gen);
    
    //        if( member(Biff,gen) ) bf = true;
	}
    
	//if( sz<8000000 && !bf ) continue;
	if( sz<15000000 && !member(noniso,ic->first) ) continue;
    
	sz /= 1000;
	dot << "\t\"" << ic->first + "CC" << "\"\t[label=\"";
	if( sz<1000 ) dot << sz << "K / ";
	else dot << sz/1000 << "." << (sz/100)%10 << "M / ";
	dot << ic->second.size() << " (";
	if( cn ) dot << cn; else dot << "X";
	dot << ")\",style=filled,fillcolor=\"" << RGBcols[(cn*6)%COL] << "\",shape=box];" << endl;
    }
#endif

#if defined(BPG_COLORED) && !defined(UNICOLCOMP)
    // output Biffurcation vertices
    for(set<string>::const_iterator ib=Biff.begin();ib!=Biff.end();++ib) {
	int cn = gen2chr[ *ib ];
    
	dot << "\t\"" << *ib << "\"\t[label=\"" << *ib;
#ifdef BPG_COMPRESSED
	dot << " (";
	if( cn ) dot << cn; else dot << "X";
	dot << ")";
#endif
	dot << "\"";

#ifdef BPG_COLORED
	dot << ",style=filled,fillcolor=\"" << RGBcols[(cn*6)%COL] << "\"";
#endif

	dot << "];" << endl;
    //    dot << "\t\"" << *ib << "\"\t[label=\"\",style=filled,fillcolor=\"" << RGBcols[cn*6] << "\"];" << endl;
    }
#endif


    // output infinity vertices
    for(int i=infv;i<0;++i) {
	dot << "\t\"" << i << "\"\t[shape=point,color=black];" << endl;
    }

//outlog << "# closing edges: " << ncls/2 << endl;
    dot << "}" << endl;
    dot.close();

#if 0
    // output table of bad blocks
    if( outbad ) {
	outlog << endl << " Bad blocks: " << endl;

	map<int,string> B;
	for(set<string>::const_iterator iv=V.begin();iv!=V.end();++iv) {
	    string x = *iv;
	    x.resize(x.size()-1);
	    istringstream is( x );
	    int n;
	    is >> n;
	    B[n] = x;
	}

	for(map<int,string>::const_iterator ib=B.begin();ib!=B.end();++ib) {
	    const string& x = ib->second;
	    outlog << x;
	    for(int i=0;i<graph.size_graph();++i) {
		outlog << " & " << (G[i].sign[x]>0?"+ ":"- ") << G[i].orf2cpan[x].first << ":" 
		//<< G[i].orf2cpan[x].second.first << "-" << G[i].orf2cpan[x].second.second << " (" 
		<< G[i].orf2cpan[x].second.second - G[i].orf2cpan[x].second.first; 
		//<< ")";
	    }
	    outlog << " \\\\" << endl;
	}
    }
#endif

#if 0
    {
        // output statistics on colored connected components
	EQ.update();
	map< string, set<string> > CC;
	EQ.get_eclasses(CC);
	cout << "Number of CCC: " << EQ.classes() << endl;
	for(map< string, set<string> >::const_iterator ic=CC.begin();ic!=CC.end();++ic) {
	    const set<string>& S = ic->second;
	    map<int,int> chr2count;
	    for(set<string>::const_iterator is=S.begin();is!=S.end();++is) {
		++chr2count[ gen2chr[*is] ];
	    }
	    for(map<int,int>::const_iterator jc=chr2count.begin();jc!=chr2count.end();++jc) {
		cout << " + " << jc->first;
	    }
	    cout << " (";
	    for(map<int,int>::const_iterator jc=chr2count.begin();jc!=chr2count.end();++jc) {
		cout << " + " << jc->second / 2;
	    }
	    cout << " blocks )" << endl;
	}
    }
#endif
} 
#endif
