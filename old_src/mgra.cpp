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
#include <set>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>   // ostream_iterator etc.
using namespace std;


const string Infty = "oo";

#define member(S,x) ((S).find(x)!=(S).end())


ofstream outlog("/dev/null");


#ifdef RGBCOLORS
#include "colors.h"
#else
vector<int> RGBcols(136);
#endif
int RGBcoeff = 1;

#include "mcolor.h"
#include "pconf.h"
#include "bpgraph.h"


///////////////////////////////////////////////

ProblemInstance PI;
MBGraph MBG;

bool ConvPhylTreeAll(MBGraph&, int Stage);

set< string > chrX;          // chromosome X (vertices)
//set< string > chrXfilt;      // filtered chromosome X genes (no microinv)

#ifdef REMOVE_SHORT_BLOCKS
#define HAVE_BLOCK_LENGTH
#include "remove_short_blocks.h"
#endif

#ifdef BPG_COMPRESSED
#define HAVE_BLOCK_LENGTH
#endif

#ifdef HAVE_BLOCK_LENGTH
map<string,size_t> BL;  // block -> length
#endif

#include "mcolor.h"
#include "2break.h"

vector<partgraph_t> RG; // recovered genomes
vector<transform_t> RT; // and transformations

bool RecoverGenomes(const transform_t&);


vector<size_t> LinChr;
vector<size_t> CirChr;


ofstream ofstat;


map< size_t, size_t > MDcount;  // MDcount[n] = # vertices of multidegree n
map< set<int>, size_t > MEcount; // ME[S] = # multiedges of multicolor S
map< set<int>, size_t > GoodMEcount; // ME[S] = # good multiedges of multicolor S
map< set<int>, size_t > ChrMEcount; // ME[S] = # good irregular multiedges of multicolor S
map< set<int>, size_t > MD2count;  // MDcount[min(S,!S)] = # simple vertices incident to S-colored
map< set<int>, size_t > MD2alone;  // MDcount[min(S,!S)] = # simple vertices incident to S-colored, with no good neighbors
map< set<int>, size_t > SimpleMEcount; // ME[S] = # simple multiedges of multicolor S
map< set<int>, size_t > SimpleCycleCount; // cycle of simple vertices
map< set<int>, size_t > SpecialCycleCount; // cycle of simple vertices and oo, of even length


void estimate(const set<string>& S) {

    MDcount.clear();
    MD2count.clear();
    MD2alone.clear();
    MEcount.clear();
    GoodMEcount.clear();
    SimpleMEcount.clear();
    ChrMEcount.clear();
    SimpleCycleCount.clear();
    SpecialCycleCount.clear();
    
    set< string > Good;

    for(set<string>::const_iterator is=S.begin();is!=S.end();++is) {
    
	const string& x = *is;
        mularcs_t M = MBG.mularcs(x);
	
	if( member(M,Infty) ) {
	    MEcount[ M[Infty] ]+=2;
	    ChrMEcount[ M[Infty] ]++;
	}
	
	if( M.size() == 2 ) {
	    Good.insert(x);
	    if( member(M,Infty) ) GoodMEcount[ M[Infty] ]+=2;
	    MD2count[ min(M.begin()->second, M.rbegin()->second) ]++;
	}
	MDcount[M.size()]++;
	
	for(map<string, set<int> >::const_iterator im=M.begin();im!=M.end();++im) {
	    if( im->first != Infty ) {
		MEcount[ im->second ]++;   // count two times
		if( M.size() == 2 ) {
		    GoodMEcount[ im->second ]++;
		    if( member(Good,im->first) ) SimpleMEcount[ im->second ]++;
		}
	    }
	}
    }

    // count lonely vertices (short paths)
    for(set<string>::const_iterator is=S.begin();is!=S.end();++is) {
        const string& x = *is;

	if( !member(Good,x) ) continue;

        map<string, set<int> > M = MBG.mularcs(x);

	if( !member(Good,M.begin()->first) && !member(Good,M.rbegin()->first) ) {
	    MD2alone[ min(M.begin()->second, M.rbegin()->second) ]++;
	}
    }

    // count cycles
    {
	set<string> processed;
	for(set<string>::const_iterator is=S.begin();is!=S.end();++is) {
    	    const string& x = *is;
	    if( member(processed, x) ) continue;


	    //if( !member(Good,x) ) continue;
	    mularcs_t Mx = MBG.mularcs(x);
	    if( Mx.size()!=2 ) continue;

	    
	    string y = x;
	    string z = "";
	    set<int> SpecialQ;
	    
	    do {
		processed.insert(y);
		mularcs_t My = MBG.mularcs(y);
		if( My.size()!=2 ) {
		    break;
		}
		if( z==My.begin()->first ) {
		    z = y;
		    y = My.rbegin()->first;
		}
		else {
		    z = y;
		    y = My.begin()->first;
		}
                
		while( y==Infty ) {
		    if( SpecialQ.empty() ) {
			SpecialQ = My[y];
			z = x;
			y = Mx.rbegin()->first;
		    }
		    else {
			if( SpecialQ != My[y] ) SpecialCycleCount[min(SpecialQ,My[y])]++;
			break;
		    }
		}
		if( y==Infty ) break;
	    } while( !member(processed,y) );
	
	    if( y==x ) SimpleCycleCount[min(Mx.begin()->second,Mx.rbegin()->second)]++;
	}
    }	    


    CirChr.clear();
    LinChr.clear();
    CirChr.resize(MBG.NumGen());
    LinChr.resize(MBG.NumGen());
    // counting chromosomes
    for(int i=0;i<MBG.NumGen();++i) {

	set<orf_t> processed;
	
	for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {

    	    const string& x = *is;
	    
	    if( member(processed,x) ) continue;
	    
	    processed.insert(x);

	    string y = MBG.OB[x];
	    while( true ) {
		if( member(processed,y) ) {
		  CirChr[i]++;
		  break;
		}
		processed.insert(y);
		if( !MBG.LG[i].defined(y) ) {
		  LinChr[i]++;
		  break;
		}
		y = MBG.LG[i][y];
		if( member(processed,y) ) {
		  CirChr[i]++;
		  break;
		}
		processed.insert(y);
		y = MBG.OB[y];
	    }
		
	    if( MBG.LG[i].defined(x) ) {
		string y = MBG.LG[i][x];
		while( !member(processed,y) ) {
		    processed.insert(y);
		    y = MBG.OB[y];
		    if( member(processed,y) ) break;
		    processed.insert(y);
		    if( !MBG.LG[i].defined(y) ) break;
		    y = MBG.LG[i][y];
		}
	    }
		    
	}
    }
}	      	






set< vertex_t > getchrset;
pair<path_t,bool> getchr(const partgraph_t& PG, const string& x) {

    path_t path;
    bool circular = false;

    getchrset.clear();
    getchrset.insert(x);

    for(vertex_t y = MBG.OB[x];;) {
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

	y = MBG.OB[y];
    }

    if( !circular && PG.defined(x) ) {
	for(string y = x;PG.defined(y);) {
	    y = PG[y];
	    getchrset.insert(y);

	    y = MBG.OB[y];
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

    for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {

	const string& x = *is;
	
	if( member(processed,x) ) continue;
        if( Xonly && !member(chrX,x) ) continue;

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

void printchr(const string& outname, const set< pair<path_t,bool> >& AllChr ) {

    ofstream out((outname+".gen").c_str());

    out << "# Genome " << outname << endl;

    string ChrTitle;
    if( PI.target.empty() ) ChrTitle = "chromosome"; else ChrTitle = "CAR";

    size_t ncirc = 0, lcirc = 0;

    for(set< pair<path_t,bool> >::const_iterator ia=AllChr.begin();ia!=AllChr.end();++ia) {

        out << endl;

        const path_t& path = ia->first;

#ifdef HAVE_BLOCK_LENGTH
        // compute length in bp
	size_t len = 0;
#endif

	for(list<orf_t>::const_iterator ip=path.begin();ip!=path.end();++ip) {
	    string t = *ip;
	    if( t[0]=='-' || t[0]=='+' ) t = t.substr(1);

#ifdef HAVE_BLOCK_LENGTH
	    len += BL[t+'t'];
#endif
	}

	if( ia->second ) {
	    ncirc++;
	    lcirc += path.size();

	    out << "# circular ";
	}
	else {
	    out << "# linear ";
	}
	out << ChrTitle << " of length " << path.size();
#ifdef HAVE_BLOCK_LENGTH
	out << " (" << len << " bp)";
#endif
	out << " follows:" << endl;

	if( (*path.begin())[0]=='+' || (*path.rbegin())[0]=='+' ) {
	    for(list<orf_t>::const_iterator ip=path.begin();ip!=path.end();++ip) {
		out << *ip << " ";
	    }
	}
	else {
	    for(list<orf_t>::const_reverse_iterator ip=path.rbegin();ip!=path.rend();++ip) {
                string e = *ip;
		if( e[0]=='-' ) e[0]='+'; else if( e[0]=='+' ) e[0]='-';
		out << e << " ";
	    }
	}

        out << "$" << endl;
    }
    out << endl << "# Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << endl;
    cout << endl << "Reconstructed genome " << outname << " has " << AllChr.size() << " " << ChrTitle << "s" << endl;

    if( ncirc ) {
	out << "#\tout of which " << ncirc << " are circular of total length " << lcirc << endl;
	cout << "\tout of which " << ncirc << " are circular of total length " << lcirc << endl;
    }

    out.close();
}



// fill in OP container with endpoints of q-obverse paths,
// starting and ending at OP

void get_obverse_paths(map< vertex_t, set<vertex_t> >& OP, const mcolor_t Q) {
    map< vertex_t, set<int> > processed;

    for(mcolor_t::const_iterator iq=Q.begin();iq!=Q.end();++iq) {

        const partgraph_t& PG = MBG.LG[*iq];

	for(map< vertex_t, set<vertex_t> >::iterator ip=OP.begin();ip!=OP.end();++ip) {

            const vertex_t& x = ip->first;

	    if( x==Infty || member(processed[x],*iq) ) continue;

	    for(vertex_t y = MBG.OB[ip->first];PG.defined(y);) {
		if( member(OP,y) ) {
		    ip->second.insert(y);
		    OP[y].insert(x);
                    processed[y].insert(*iq);
		    break;
		}
		y = PG[y];

		if( y==x ) {
		    //already is circular, we don't care
		    break;
		}

		y = MBG.OB[y];
	    }
	    processed[x].insert(*iq);
	}
    }
}





set<string> ToLabel;





// Save .dot file and output statistics of synteny blocks representing breakpoints
const set<string> EmptySS;
//void save_dot(const string& dotname, const set<string>& VV = EmptySS) {
void save_dot(const string& dotname, const bool outbad = false, const set<string>& VV = EmptySS) {

    ofstream dot(dotname.c_str());
    dot << "graph {" << endl;
    dot << "edge [colorscheme=" << PI.colorscheme << "];" << endl;
    //dot << "node [fontname=\"freefont/FreeSans\",fontsize=8];" << endl;

#ifdef BPG_COLORED

    map<string,int> gen2chr;
    equivalence<string> SameColor;
    set<string> Biff = MBG.VertexSet;
    set<string> noniso; // non-isolated vertices
    
    for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
	const string& x = *is;
	if( !member(Biff,x) ) continue;

	string xx = x;
	xx.resize(xx.size()-1);
	string chr = G[REFGENOME].orf2cpan[xx].first;
	if(chr=="chrX") chr = "0";
	for(int i=0;i<chr.size();++i) if(chr[i]<'0'||chr[i]>'9') chr[i]=' ';
	istringstream istr(chr); int cn; istr >> cn;
	gen2chr[x] = cn;
        gen2chr[MBG.OB[x]] = cn;

#ifndef BPG_COMPRESSED
	mularcs_t Mx = MBG.mularcs(x);
	if( Mx.size() == 1 ) Biff.erase(x);
    }
#else

	SameColor.addrel(x,MBG.OB[x]);

	mularcs_t Mx = MBG.mularcs(x);
#ifndef UNICOLCOMP
	if( Mx.size() == 1 ) {
	    const string& y = Mx.begin()->first;
#else
	for(mularcs_t::const_iterator im=Mx.begin();im!=Mx.end();++im) {
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
    for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
	const string& x = *is;

	mularcs_t Mx = MBG.mularcs(x);
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

		dot << "\t\"" << xt << "\"\t--\t\"" << xh << "\"\t[style=dashed,color=" << MBG.NumGen()+1 << "];" << endl;
	    }
            // EQ.addrel(x,xh);
#endif  // UNICOLCOMP

	}
#endif  // FULL_BPG



	//if( Mx.size()==N ) continue;

	for(mularcs_t::const_iterator im=Mx.begin();im!=Mx.end();++im) {
	    const string& y = im->first;

	    if( member(V,y) ) continue; // already output

#if 0
	    mularcs_t Mx1 = MBG.mularcs(x);
	    Mx1.erase(Infty);
	    if( Mx1.size() >= 3 ) continue;
	    if( y!=Infty ) {
		mularcs_t My = MBG.mularcs(y);
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

#if 0
    dot << "subgraph clusterX {" << endl;
    for(set<string>::const_iterator is=S.begin();is!=S.end();++is) {
	const string& x = *is;
        if( MBG.mularcs(x).size()==1 ) continue;
	if( member(chrX,x) ) {
	    if( member(VV,x) ) {
		dot << "\t\"" << x << "\"\t[style=filled,fillcolor=7];" << endl;
	    }
	    else {
		dot << "\t\"" << x << "\"\t[style=filled,fillcolor=8];" << endl;
	    }
	}
	else if( member(VV,x) ) {
		dot << "\t\"" << x << "\"\t[style=filled,fillcolor=9];" << endl;
	}
	if( !VV.empty() && member(ToLabel,x) ) {
	    dot << "\t\"" << x << "\"\t[style=filled,fillcolor=7];" << endl;
	}
    }
    dot << "}" << endl;
#endif

//	    outlog << "# closing edges: " << ncls/2 << endl;
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
	    for(int i=0;i<MBG.NumGen();++i) {
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



void save_dot_path(const string& dotname) {

    // count H-subgraphs
    map< pair<mcolor_t,mcolor_t>, size_t > Hcount;
    // middle edge is T-consistent?
    map< pair<mcolor_t,mcolor_t>, bool > Hmid;


    //ofstream dot(dotname.c_str());
    //dot << "graph {" << endl;
    //dot << "graph [overlap=false,splines=true];" << endl;
    //dot << "node [fontname=\"freefont/FreeSans\",fontsize=8];" << endl;
    //dot << "edge [fontname=\"freefont/FreeSans\",fontsize=6,colorscheme=" << PI.colorscheme << "];" << endl;

    // go over multiedges

    int infv = 0;
    set< orf_t > V; // vertex set
    for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
	const string& x = *is;

	const mularcs_t Mx = MBG.mularcs(x);
	if( Mx.size() == 1 ) continue; // trivial cycle
    
	for(mularcs_t::const_iterator im=Mx.begin();im!=Mx.end();++im) {
	    const string& y = im->first;
	    const set<int>& C = im->second;

	    if( member(V,y) ) continue; // already output

	    if( y==Infty ) {
		--infv;
		//dot << "\t\"" << infv << "\"\t[shape=point,color=black];" << endl;
	    }

	    bool hgraph = false;
	    while(true) {
		if( Mx.size()!=3 ) break;
		if( y==Infty ) break;
		mularcs_t My = MBG.mularcs(y);
		if( My.size()!=3 ) break;
		mularcs_t Mx0 = Mx;
		Mx0.erase(y);
		My.erase(x);
		const mcolor_t& Q1 = Mx0.begin()->second;
		const mcolor_t& Q2 = Mx0.rbegin()->second;
		if( (Q1 == My.begin()->second) || (Q2 == My.begin()->second) ) {

		    mcolor_t QQ1 = MBG.CColorRep(Q1);
		    mcolor_t QQ2 = MBG.CColorRep(Q2);

		    Hcount[make_pair(QQ1,QQ2)]++;
		    Hcount[make_pair(QQ2,QQ1)]++;
		    Hmid[make_pair(QQ2,QQ1)] = member(MBG.Color,im->second);

		    if( !member(MBG.Color,Q1) || !member(MBG.Color,Q2) ) break;

		    hgraph = true;
		    //dot << "\t\"" << x << "\"\t--\t\"" << y << "\"\t[" << "color=" << MBG.OColor[Q1] << "];" << endl;
		    //dot << "\t\"" << x << "\"\t--\t\"" << y << "\"\t[" << "color=" << MBG.OColor[Q2] << "];" << endl;

		}
		break;
	    }
#if 0
	    if( member(MBG.Color,C) ) {
		dot << "\t\"" << x << "\"\t--\t\"";
		if( y==Infty ) dot << infv << "\"\t[len=0.5,";
		else dot << y << "\"\t[";
		if( hgraph ) dot << "penwidth=2,bold,";
		dot << "color=" << MBG.OColor[C] << "];" << endl;
	    }
	    else if( !hgraph ) {    
		// output individual edges
		for(set<int>::const_iterator ic=C.begin();ic!=C.end();++ic) {
		    // output edge (x,y)
		    dot << "\t\"" << x << "\"\t--\t\"";
		    if( y==Infty ) dot << infv << "\"\t[len=0.5,";
		    else dot << y << "\"\t[";
		    dot << "style=dashed,color=" << 1 + *ic << "];" << endl;
		}
	    }
#endif
	}
	V.insert(x);
    }
    //dot << "}" << endl;
    //dot.close();

    // output H-subgraphs count
    ofstat << endl << "% Fair multi-edges count: " << endl << endl;

    /*
    // the set of multicolors that appear in H-subraphs
    set< mcolor_t > HCrow,HCcol;
    */

    list< mcolor_t > HCrow;
    set< mcolor_t > proc;

#if defined(MRDQHC) && !defined(RESTRICT2MDH)

    string myord[] = { "0", "1", "01", "2", "012", "3", "45", "4", "5" };
    for(int i=0;i<9;++i) {
	mcolor_t Q;
	for(size_t k=0;k<myord[i].size();++k) {
	    int g = myord[i][k] - '0';
	    Q.insert( g );
	}
	mcolor_t QQ = MBG.CColorRep(Q);
	HCrow.push_back(QQ);
        proc.insert(QQ);
    }


#endif

#ifdef MRODQHC
    string myord[] = { "0", "1", "01", "6", "016", "2", "0162", "3", "45", "4", "5" };
    for(int i=0;i<11;++i) {
	mcolor_t Q;
	for(size_t k=0;k<myord[i].size();++k) {
	    int g = myord[i][k] - '0';
	    Q.insert( g );
	}
	mcolor_t QQ = MBG.CColorRep(Q);
	HCrow.push_back(QQ);
        proc.insert(QQ);
    }
#endif




#ifndef HG_TONLY
    // adding other colors
    for(map< pair<mcolor_t,mcolor_t>, size_t >::const_iterator ih=Hcount.begin();ih!=Hcount.end();++ih) {
	if(!member(proc,ih->first.first)) {
	    HCrow.push_back(ih->first.first);
	    proc.insert(ih->first.first);
	}
	if(!member(proc,ih->first.second)) {
	    HCrow.push_back(ih->first.second);
	    proc.insert(ih->first.second);
	}
    }
#endif


    const list< mcolor_t >& HCcol = HCrow;



    ofstat << "\\begin{table}[h]" << endl;
    //ofstat << "\\begin{footnotesize}" << endl;
    ofstat << "\\centering \\begin{tabular}{|c||";
    for(int i=0;i<HCcol.size();++i) ofstat << "c|";
    ofstat << "}" << endl;
    ofstat << "\\hline" << endl;

    for(list< mcolor_t >::const_iterator ic=HCcol.begin();ic!=HCcol.end();ic++) {
	const mcolor_t& Q1 = *ic;
	ofstat << " & ${";
	if( member(MBG.Color,Q1) ) ofstat << "\\bf ";
//	ofstat << set2name(Q1) << "+" << set2name(CColor(Q1)) << "}";
	ofstat << MBG.set2name(Q1) << "+}$";
    }
    ofstat << "\\\\" << endl;
    ofstat << "\\hline \\hline" << endl;

    for(list< mcolor_t >::const_iterator ic=HCrow.begin();ic!=HCrow.end();ic++) {
	const mcolor_t& Q1 = *ic;

	ofstat << "${";
	if( member(MBG.Color,Q1) ) ofstat << "\\bf ";
	ofstat <<  MBG.set2name(Q1) << "+}$"; // << "+" << set2name(CColor(Q1)) << "}";

	for(list< mcolor_t >::const_iterator jc=HCcol.begin();jc!=HCcol.end();jc++) {
	    const mcolor_t& Q2 = *jc;
	    ofstat << " & ";
	    //if( Q2 <= Q1 ) continue;

	    if( member(Hcount,make_pair(Q1,Q2)) ) {
		ofstat << " "; //"${";
		if( MBG.AreAdjacentBranches(Q1,Q2) )
		    //if( Hmid[make_pair(Q1,Q2)] )
		    //ofstat << "\\bf ";
		    ofstat << "{\\cellcolor[gray]{.9}}";
#ifdef HG_TONLY
		if( member(BranchLen,set2name(Q1)) && member(BranchLen,MBG.set2name(Q2)) ) {
		    long t = (long)(1e5*Hcount[make_pair(Q1,Q2)]/(double)BranchLen[set2name(Q1)]/(double)BranchLen[set2name(Q2)]);
		    ofstat << "$" << t << "\\cdot 10^{-5}$ ";
		}
#else
		ofstat << Hcount[make_pair(Q1,Q2)] << " "; // << "}$";
#endif
	    }
	    if( Q1==Q2 ) {
		ofstat << "$\\star$";
	    }
	}
	ofstat << " \\\\" << endl;
	ofstat << "\\hline" << endl;
    }
    //ofstat << "\\hline" << endl;
    ofstat << "\\end{tabular}" << endl;
    //ofstat << "\\end{footnotesize}" << endl;
    ofstat << "\\end{table}" << endl;
    ofstat << endl;
}







// can incident multiedges of x form multicolor Q (current don't check T-consistent formation)
// if return false, then Q cannot be formed
// if true - who knows...

bool canformQoo = true; // safe choice, at later stages may change to false
//bool canformQoo = false;

bool canformQ(const string& x, const mcolor_t& Q) {
    if( x==Infty ) {
	//cout << "ERROR: canformQ(oo)" << endl;
	//exit(1);
	return canformQoo;
    }

    // color Q can be formed if some it adjacent multicolors form a partition of Q
    // OR 
    // if every intersection Q \cap QQ = \emptyset or QQ.

    mulcols_t M = MBG.mulcols(x);
    for(mulcols_t::const_iterator im=M.begin();im!=M.end();++im) {
	mcolor_t C;
	set_intersection(Q.begin(),Q.end(),im->first.begin(),im->first.end(),inserter(C,C.begin()));
	if( C.size()>0 && C.size()<im->first.size() ) return false;
    }
    return true;
}



void print_cc_xlen() {
    // count connected components in the graph
    equivalence<vertex_t> C;

    size_t nc = 0;

    ofstat << "... complete multiedges:";
    for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
	const string& x = *is;
	C.addrel(x,x);
	mularcs_t M = MBG.mularcs(x);
	if( M.size()==1 && M.begin()->second.size()==MBG.NumGen() && (x < M.begin()->first || M.begin()->first==Infty) ) {
	    ofstat << " " << x << "~" << M.begin()->first;
	    ++nc;
	}
    }
    ofstat << "\t(total: " << nc << ")" << endl;

    for(int i=0;i<MBG.NumGen();++i) {
	for(partgraph_t::const_iterator il=MBG.LG[i].begin();il!=MBG.LG[i].end();++il) {
	    C.addrel(il->first,il->second);
	}
    }

    C.update();

    map< vertex_t,set<vertex_t> > cls;
    C.get_eclasses(cls);

    map< size_t, size_t > stx;
    for(map< vertex_t,set<vertex_t> >::const_iterator ic=cls.begin();ic!=cls.end();++ic) {
	stx[ ic->second.size() ]++;
    }

    ofstat << "... connected components:";
    for(map< size_t,size_t >::const_iterator is=stx.begin();is!=stx.end();++is) {
	ofstat << " " << is->first << "^" << is->second;
    }
    ofstat << endl;

#ifdef HAVE_BLOCK_LENGTH
    // filter out microiversions in X chromosome
    //set< string > chrXfilt;
    ofstream out("Xlen.res");
    for(map< vertex_t,set<vertex_t> >::const_iterator ic=cls.begin();ic!=cls.end();++ic) {
	const set<vertex_t>& V = ic->second;
	if( member(chrX,*V.begin()) ) {
	    if( V.size()<=4 ) {
		set<string> T;
		for(set<vertex_t>::const_iterator iv=V.begin();iv!=V.end();++iv) {
		    string x = *iv;
		    x.resize(x.size()-1);
		    T.insert(x);
		}
		if( T.size()<V.size() ) continue;
	    }
	    for(set<vertex_t>::const_iterator iv=V.begin();iv!=V.end();++iv) {
		string x = *iv;
		x.resize(x.size()-1);
		//chrXfilt.insert(x);
		out << x << "\t" << BL[*iv] << endl;
	    }
	}
    }
    out.close();
#endif
}


void histStat() {
    ofstat << endl << "Total number of 2-breaks: " << TwoBreak::History.size() << endl;
    { // output statistics
	map< mcolor_t, size_t > n2br;
	for(list<TwoBreak>::const_iterator il=TwoBreak::History.begin();il!=TwoBreak::History.end();++il) {
	    n2br[il->MultiColor]++;
	}
	for(map< mcolor_t, size_t >::const_iterator im=n2br.begin();im!=n2br.end();++im) {
	    ofstat << MBG.set2name(im->first) << "\t" << im->second << endl;
	}
    }
    ofstat << endl;
}




/* Given a non-linear genome PG of a multicolor Q and a transformation into a linear genome,
 * find linearizing fissions in the transformation, move them to beginning, and apply to PG
 * i.e., try to reorder transformation into: PG -> PG' -> linear genome, where PG' is linear
 * and obtained from PG with fission.
 * We replace PG with PG' and return the transformation PG -> PG'
 * Transformation may contain only multicolors Q' with Q'\cap Q = 0 or Q.
*/

transform_t decircularize(partgraph_t& PG, transform_t& TG, const mcolor_t& Q) {

    // decircularizing sub-transform that is removed
    transform_t D;

    size_t CircSize = numchr(PG).second;
    if( CircSize == 0 ) return D;

    outlog << "Eliminating " << CircSize << " circular chromosomes in " << MBG.set2name(Q) << endl;

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
    for(transform_t::iterator it=start;it!=TG.end();) {

        // check multicolor
	{
	    mcolor_t C;
	    set_intersection( it->MultiColor.begin(), it->MultiColor.end(), Q.begin(), Q.end(), inserter(C,C.begin()) );
    
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

	    mcolor_t C;
	    set_intersection( t.MultiColor.begin(), t.MultiColor.end(), s.MultiColor.begin(), s.MultiColor.end(), inserter(C,C.begin()) );
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
		mcolor_t C;
		set_intersection( kt->MultiColor.begin(), kt->MultiColor.end(), Q.begin(), Q.end(), inserter(C,C.begin()) );
    
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



void printstat(int st) {

    //ofstat << "% == After Stage " << st << ":" << endl << endl;

    estimate(MBG.VertexSet);
 
    size_t other = 0;
    
    multimap< size_t, string > MMM;
    for(map< set<int>, size_t >::const_iterator im=MEcount.begin();im!=MEcount.end();++im) {
 
 /*
	if( !member(Color,im->first) ) {
	    other += m1;
	    continue;
	}
 */
	set<int> C = MBG.CColor(im->first);  // complementary multicolor
	if(C>im->first) continue;
	
	size_t m1 = MEcount[C]/2;
	size_t m2 = (im->second)/2;
 
	size_t paths = MD2count[C] - (SimpleMEcount[C] + SimpleMEcount[im->first]) - MD2alone[C] - SpecialCycleCount[C];
	size_t cycles = SimpleCycleCount[C] + SpecialCycleCount[C];
 
	ostringstream os;
 
	os << "{";		    
	if( member( MBG.Color, im->first ) ) os << "\\bf ";
	os << MBG.set2name(MBG.CColorRep(C)) << " + " << MBG.set2name(MBG.CColor(MBG.CColorRep(C))) << "} & "
	// multiedges
	<< m1 << " + " << m2 << " = " << m1+m2 << " & " 
	// simple vertices
	<< MD2count[C] << " & " 
	// simple multiedges
	<< SimpleMEcount[C] << " + " << SimpleMEcount[im->first] << " = " << SimpleMEcount[C] + SimpleMEcount[im->first] << " & "
	// simple paths + cycles
	<< paths << " + " << cycles << " = " << paths+cycles << " & "
	// irregular multiedges
	<< ChrMEcount[C] << " + " << ChrMEcount[im->first] << " = " << ChrMEcount[C] + ChrMEcount[im->first];
	
	MMM.insert( make_pair( m1+m2, os.str() ) );
    }
 
    /*
    {
	ostringstream os;
	os << other;
	MMM.insert( make_pair( other, make_pair( "others",  os.str() ) ) );
    }
    */

    ofstat << endl;
    ofstat << "% Rearrangement characters:" << endl << endl;
    ofstat << "\\begin{table}[h]" << endl;
    //ofstat << "\\begin{footnotesize}" << endl;
    ofstat << "\\centering \\begin{tabular}{|c|c|c|c|c|c|}" << endl;
    ofstat << "\\hline" << endl;
    ofstat << "Multicolors & multiedges & simple vertices & simple multiedges & simple paths+cycles & irreg. multiedges\\\\" << endl;
    ofstat << "\\hline" << endl;
    for(multimap< size_t, string >::reverse_iterator im=MMM.rbegin();im!=MMM.rend();++im) {
	ofstat << im->second << "\\\\" << endl;
    }
    ofstat << "\\hline" << endl;
    ofstat << "\\end{tabular}" << endl;
    //ofstat << "\\end{footnotesize}" << endl;
    ofstat << "\\end{table}" << endl;
    ofstat << endl;
 
 
 
    ofstat << "% Estimated distances:" << endl << endl;
    ofstat << "\\begin{table}[h]" << endl;
    ofstat << "\\centering \\begin{tabular}{|c|";
    for(size_t i=0;i<MBG.NumGen();++i) ofstat << "c|";
    ofstat << "}" << endl;
    ofstat << "\\hline" << endl;
    ofstat << "Stage " << st;
    for(size_t i=0;i<MBG.NumGen();++i) ofstat << " & " << MBG.num2gen[i];
    ofstat << " \\\\" << endl << "\\hline" << endl;
    for(size_t i=0;i<MBG.NumGen();++i) {
	ofstat << MBG.num2gen[i] << " & ";
	for(size_t j=0;j<MBG.NumGen();++j) {
	    if( j>i ) {
		ofstat << genome_dist(MBG.LG[i],MBG.LG[j],MBG.OB)[2];
	    }
	    if(j<MBG.NumGen()-1) ofstat << " & ";
	    else ofstat << "\\\\";
	}
	ofstat << endl;
    }
    ofstat << "\\hline" << endl;
    ofstat << "\\end{tabular}" << endl;
    ofstat << "\\end{table}" << endl;
    ofstat << endl;

}

















////////////////////////////////////////////////////////////////////////




int main(int argc, char* argv[]) {


    cout << "MGRA (Multiple Genome Rearrangements & Ancestors) ver. 1.1" << endl;
    cout << "(c) 2008,12 by Max Alekseyev <maxal@cse.sc.edu>" << endl;
    cout << "Distributed under GNU GENERAL PUBLIC LICENSE license." << endl;
    cout << endl;

    if( argc != 2 ) {
	cout << endl << "Usage: mgra <ProblemConfiguration>" << endl;
        return 1;
    }

    // reding problem configuration
    if( !ReadProblem( argv[1], PI ) ) return 1;


#ifdef RGBCOLORS
    for(size_t i=0;i<COL;++i) RGBcols[i] = "\"" + RGBcols[i] + "\"";
    RGBcoeff = (COL-1)/(PI.NumGen-1);
#else
    for(size_t i=0;i<RGBcols.size();++i) RGBcols[i] = i+1;
    RGBcoeff = 1;
#endif


    vector<Genome> G(PI.NumGen);
    if( PI.blkformat == "infercars" ) {
	read_infercars(PI,G);
    }
    else if( PI.blkformat == "grimm" ) {
	read_grimm(PI,G);
    }
    else {
	cerr << "ERROR: unknown synteny blocks format" << endl;
	return 2;
    }

    if( PI.NumGen < 2 ) {
	cerr << "ERROR: at least two input genomes required" << endl;
	return 3;
    }

    MBG.init(G,PI.trees,PI.target);



///////////////////////////////////////////////////////////////////////////////

#ifdef REMOVE_SHORT_BLOCKS
    {
	set<vertex_t> ToDel;
	get_short_blocks(MBG,RG,BL,ToDel);

	while( !ToDel.empty() ) {
	    const vertex_t& x = *ToDel.begin();
            const vertex_t& y = MBG.OB[x];

	    for(size_t i=0;i<MBG.NumGen();++i) {
		string b0, b1;
		if( MBG.LG[i].defined(x) ) {
		    b0 = MBG.LG[i][x];
		    MBG.LG[i].erase(x);
		}
		if( MBG.LG[i].defined(y) ) {
		    b1 = MBG.LG[i][y];
		    MBG.LG[i].erase(y);
		}
		if( !b0.empty() && !b1.empty() ) {
		    MBG.LG[i].insert(b0,b1);
		}
	    }
	    MBG.VertexSet.erase(x);
	    MBG.VertexSet.erase(y);
	    MBG.OB.erase(x);
	    
	    ToDel.erase(x);
	    ToDel.erase(y);
	}
	cout << "Blocks left: " << MBG.VertexSet.size()/2 << endl;
    }
#endif

    ofstat.open("stats.txt");


    ofstat << "Initial graph:" << endl;
    print_cc_xlen();
    printstat(0);


    save_dot(PI.graphfname+"0.dot");
    save_dot_path(PI.graphfname+"0p.dot");


    bool dots1 = true;
    bool dots2 = true;
    bool dots3 = true;
    bool dots4 = true;


#if 0
    ofstream ftable("table.res");


    ftable << "\\begin{table}[h]" << endl;
    ftable << "\\centering \\begin{tabular}{|c|c|";
    for(size_t i=1;i<=MBG.NumGen();++i) ftable << "c|";
    ftable << "}" << endl;
    ftable << "\\hline" << endl;
    ftable << "Stage & 2-Breaks ";
    for(size_t i=1;i<=MBG.NumGen();++i) ftable << "& $v_{" << i << "}(G)$ ";
    ftable << "\\\\" << endl;
    ftable << "\\hline" << endl;
#endif


    while(true) {

	static size_t st = 0;


#if 0
	//if( !NumRea.empty() ) {
	if( false ) {
	    ftable << st << " & ";
	    for(int i=0;i<4;++i) {
		//ftable << NumRea[i];
		switch(i) {
		case 0: ftable << "("; break;
		case 1: ftable << ","; break;
		case 2: ftable << "):"; break;
		}
	    }
	}
	else ftable << " & ";

//	    ftable << " & ";
//	    for(int i=0;i<N;++i) {
//	      ftable << LinChr[i] << "/" << CirChr[i] << " ";
//	    }
//	    ftable << " & " << MDcount[1];
	
	for(int i=1;i<=MBG.NumGen();++i) ftable << " & " << MDcount[i];
	ftable << " \\\\" << endl;
#endif





#if 0
	for(int i=0;i<MBG.NumGen();++i) {
	    set< pair<path_t,bool> > GN;
	    splitchr(MBG.LG[i], GN);
	    printchr(sname[i].substr(0,1) + "_m",GN);
	}
	break;
#endif



	if( PI.stages >= 1 ) {

	    st = 1;
	    if( ConvPhylTreeAll(MBG,1) ) continue;
	    //if( ConvPhylTreeAll(MBG,26) ) continue;  // cut hanging free ends
    
	    if( dots1 ) {
		dots1=false;
    
		ofstat << "After Stage 1 graph:" << endl;
		print_cc_xlen();
                printstat(1);

		save_dot(PI.graphfname+"1.dot");
		save_dot_path(PI.graphfname+"1p.dot");


		/*
		size_t subtle = 0;
    
		partgraph_t PG;
		ofstream cf("st1comp.res");
		for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
		    const string& x = *is;
		    if( PG.defined(x) ) continue;
    
		    string y = Infty;
		    bool good = true;
		    int def = 0;
		    for(int i=0;i<PI.target.size();++i) {
			int j = // PI.target[i] - '0';
			if( MBG.LG[j].defined(x) ) {
			    def++;
			    if( y==Infty ) y=MBG.LG[j][x];
			    if(y!=MBG.LG[j][x] ) good = false;
			}
		    }
		    if( good && def==PI.target.size() && y!=Infty ) {
			PG.insert(x,y);
			cf << x << "\t" << y << endl;
    
			if( !member(MBG.LeafAdj,make_pair(x,y)) ) subtle++;
		    }
		}
		cf.close();
    
		cout << "Subtle adjacencies: " << subtle << " / " << PG.size()/2 << endl;
		*/
    
	    }
	}


/////////////////////////////////// STAGE 2 /////////////////////////////////////////

	if( PI.stages >= 2 ) {

	    st = 2;

	    if( ConvPhylTreeAll(MBG,2) ) continue;

	    /* if( FairPathSearch(2) ) continue; */
    
	    if( dots2 ) {
		dots2=false;
    
		ofstat << "After Stage 2 graph:" << endl;
		print_cc_xlen();
		printstat(2);

		save_dot(PI.graphfname+"2.dot");
		save_dot_path(PI.graphfname+"2p.dot");
    
	    }
	}

/////////////////////////////////// STAGE 3 /////////////////////////////////////////
/// STAGE 3, somewhat unreliable


	if( PI.stages >= 3 ) {

	    outlog << "Stage: 3" << endl;

	    st = 222;
	    if( ConvPhylTreeAll(MBG,222) ) continue;   // cut the graph into connected components
	    st = 2222;
	    if( ConvPhylTreeAll(MBG,2222) ) continue;  // process 4-cycles
    
	    if( canformQoo ) {
		canformQoo = false; // more flexible
		continue;
	    }
    
	    /*
	    static bool hotfix = true;
	    if( hotfix ) {
		hotfix = false;
		set<int> M; M.insert(0);
		TwoBreak("127h","792h","671h","127t",M).apply(MBG,true);
		continue;
	    }
	    */
    
	    if( dots3 ) {
		dots3=false;

		ofstat << "After Stage 3 graph:" << endl;
		print_cc_xlen();
		printstat(3);

		save_dot(PI.graphfname+"3.dot",true);
		save_dot_path(PI.graphfname+"3p.dot");
	    }
	}

#if defined(LEAVE4GENOMES) || defined(RESTRICT2MDH) || defined(STAGES12)
	break;
#endif


	if( PI.stages >= 4 ) {

	    outlog << "Stage: 4" << endl;

	    st = 4;
	    MBG.SplitBadColors = true;
	    outlog << "SplitBadColors is ON" << endl;
	    if( ConvPhylTreeAll(MBG,2) ) {
		MBG.SplitBadColors = false;
		outlog << "SplitBadColors is OFF" << endl;
		continue;
	    }
	    MBG.SplitBadColors = false;
	    outlog << "SplitBadColors is OFF" << endl;

	    if( dots4 ) {
		dots4=false;

		ofstat << "After Stage 4 graph:" << endl;
		print_cc_xlen();
		printstat(4);

		save_dot(PI.graphfname+"4.dot",true);
		save_dot_path(PI.graphfname+"4p.dot");

	    }
	}



	static bool stage4 = true;

	if( stage4 && !PI.completion.empty() ) {

	    stage4 = false;

	    outlog << "Manual Completion Stage" << endl;

	    for(list< vector<string> >::const_iterator il=PI.completion.begin(); il!=PI.completion.end(); ++il) {

		TwoBreak t;
		t.OldArc[0].first = (*il)[0];
		t.OldArc[0].second = (*il)[1];
		t.OldArc[1].first = (*il)[2];
		t.OldArc[1].second = (*il)[3];
		t.MultiColor = MBG.name2set( (*il)[4] );

                t.apply(MBG,true);

	    }

            continue;
	}

	break;
    }	

#if 0
    ftable << "\\hline" << endl;
    ftable << "\\end{tabular}" << endl;
    ftable << "\\end{table}" << endl;
    ftable << endl;
    ftable.close();
#endif


    save_dot(PI.graphfname+"99.dot",true);
    save_dot_path(PI.graphfname+"99p.dot");



#ifndef OPOSSUM
    histStat();

    /*
    cout << "Chromosome X 2-breaks: " << TwoBreak::HistoryX.size() << endl;
    { // output statistics
	map< mcolor_t, size_t > n2br;
	for(list<TwoBreak>::const_iterator il=TwoBreak::HistoryX.begin();il!=TwoBreak::HistoryX.end();++il) {
	    n2br[CColorRep(il->MultiColor)]++;
	}
	for(map< mcolor_t, size_t >::const_iterator im=n2br.begin();im!=n2br.end();++im) {
	    cout << set2name(im->first) << "+" << set2name(CColor(im->first)) << "\t" << im->second << endl;
	}
    }
    */

    //save_dotX();
#endif



    if( !PI.target.empty() ) {
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
	for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
	    const string& x = *is;
	    if( PG.defined(x) ) continue;

	    string y = Infty;
	    bool good = true;
	    int def = 0;
	    for(int i=0;i<PI.target.size();++i) {
		int j = MBG.gen2num[PI.target.substr(i,1)];
		if( MBG.LG[j].defined(x) ) {
		    def++;
		    if( y==Infty ) y=MBG.LG[j][x];
		    if( y!=MBG.LG[j][x] ) good = false;
		}
	    }
	    if( good && def==PI.target.size() && y!=Infty ) {
                PG.insert(x,y);
		//cf << x << "\t" << y << endl;
	    }
	}
	//cf.close();

	set< pair<path_t,bool> > GN;
	splitchr(PG, GN);
	printchr(PI.target, GN);

	//for(int i=0;i<N;++i) {
	//    set< pair<path_t,bool> > GN;
	//    splitchr(MBG.LG[i], GN);
	//    printchr(sname[i].substr(0,1) + "_m",GN);
	//}


    }
    else {  /* empty target */

	const size_t NC = MBG.DiColor.size();
    
	RG.resize(NC);
	RT.resize(NC);
    
	if( !RecoverGenomes(TwoBreak::History) ) exit(1);
    
	// T-transformation complete, we procede with recovering the ancestral genomes
    
	outlog << "Initial 2-break distances from the root X: " << endl;
	for(int i=0;i<NC;++i) {
	    outlog << MBG.set2name(MBG.TColor[i]) << ":\t" << RT[i].size() << endl;
	}
    
	// FIXME: check that the order in which circular chromosomes are eliminated
    
	for(int i=0;i<NC;++i) {
    
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
	for(int i=0;i<NC;++i) {
	    outlog << MBG.set2name(MBG.TColor[i]) << ":\t" << RT[i].size() << endl;
	}
    
	for(int i=0;i<NC;++i) {
    
	    set< pair<path_t,bool> > GN;
	    splitchr(RG[i], GN);
	    printchr(MBG.set2name(MBG.TColor[i]),GN);
    
	    //splitchr(RG[i], GN, true);
	    //printchr(MBG.set2name(MBG.TColor[i]) + "_x",GN);
    
	    ofstream tr( (MBG.set2name(MBG.TColor[i]) + ".trs").c_str() );
	    //ofstream trx( (MBG.set2name(MBG.TColor[i]) + "_x.trs").c_str() );
	    for(transform_t::const_iterator it=RT[i].begin();it!=RT[i].end();++it) {
		tr << it->OldArc[0].first << " " << it->OldArc[0].second << "\t" << it->OldArc[1].first << " " << it->OldArc[1].second << "\t" << MBG.set2name(it->MultiColor) << endl;
	    //   if( member(chrX,it->OldArc[0].first) || member(chrX,it->OldArc[0].second) || member(chrX,it->OldArc[1].first) || member(chrX,it->OldArc[1].second) ) {
	    //	trx << "(" << it->OldArc[0].first << "," << it->OldArc[0].second << ")x(" << it->OldArc[1].first << "," << it->OldArc[1].second << "):{" << MBG.set2name(it->MultiColor) << "} " << endl;
	    //    }
	    }
	    tr.close();
	    //trx.close();
	}

    }


    ofstat.close();

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

    for(int i=0;i<MBG.NumGen()-1;++i) {
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
    vector< vector<size_t> > RTF(NC);
    for(size_t j=0;j<NC;++j) RTF[j].resize(3);

    for(transform_t::const_reverse_iterator it=tr.rbegin();it!=tr.rend();++it) {
	const mcolor_t& Q = it->MultiColor;

	outlog << "Reverting (" << it->OldArc[0].first << "," << it->OldArc[0].second << ")x(" << it->OldArc[1].first << "," << it->OldArc[1].second << "):{" << MBG.set2name(it->MultiColor) << "} " << " in";

	for(int i=0;i<NC;++i) {

	    if( !includes( Q.begin(),Q.end(),MBG.TColor[i].begin(),MBG.TColor[i].end() ) ) continue;

            // TColor[i] is subset of Q

            size_t nchr_old = 0;
	    if( Q == MBG.TColor[i] ) {
		nchr_old = numchr(RG[i]).first;
	    }


	    it->revertSingle(RG[i]);

	    if( Q==MBG.TColor[i] ) {
		outlog << " " << MBG.set2name(MBG.TColor[i]);
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
    for(size_t j=0;j<NC;++j) {
	outlog << MBG.set2name(MBG.TColor[j]) << "+" << MBG.set2name(MBG.CColor(MBG.TColor[j])) << "\t&\t" << RTF[j][0] << " & " << RTF[j][1] << " & " << RTF[j][2]
	    << " &\t" << RTF[j][0]+RTF[j][1]+RTF[j][2] << " \\\\" << endl;
        outlog << "\\hline" << endl;
	tot[0] += RTF[j][0];
	tot[1] += RTF[j][1];
	tot[2] += RTF[j][2];
    }
    outlog << "Total\t&\t" << tot[0] << " & " << tot[1] << " & " << tot[2] << " &\t" << tot[0]+tot[1]+tot[2] << " \\\\" << endl;

    return true;
}


/*
typedef couple<vertex_t> edge2_t;
map< pair< mcolor_t, edge2_t >, size_t > EdgeW;

size_t PathLength(const path_t& path, const mcolor_t& QQ) {

    size_t w = 0;

    path_t::const_iterator z0 = path.begin();
    path_t::const_iterator z1 = ++path.begin();

    while( z1!=path.end() ) {
	const pair< mcolor_t, edge2_t > p = make_pair(QQ,edge2_t(*z0,*z1));
	if( member(EdgeW,p) ) {
	    w += EdgeW[p];
	}
	else w++;
	z0 = z1++;
    }

    return w;

    //return p.size();
}
*/


int ProcessSimplePath(path_t path) {
    int nr = 0;
    if( path.size() >= 4 || (path.size()==3 && *path.begin()==*path.rbegin()) ) {
	outlog << endl << "Processing a path of length " << path.size()-1 << endl;

	outlog << "path:\t" << *path.begin();
	for(list<string>::const_iterator ip=++path.begin();ip!=path.end();++ip) {
	    outlog << " -- " << *ip;
	}
	outlog << endl;

	if( member(chrX,*(++path.begin())) ) {
	    outlog << "... on X chromosome" << endl;
	}

	if( path.size()%2 && (*path.begin()!=*path.rbegin()) ) {
            outlog << "... ";
	    if( !member(MBG.DiColor,MBG.mularcs(*(++path.begin()))[*path.begin()]) ) {
	       //*path.begin()==Infty ) {
		path.erase(path.begin());
		outlog << "left";
	    }
	    else {
		path.erase(--path.end());
		outlog << "right";
	    }
            outlog << " end removed" << endl;
	}

        /*
	map< vertex_t, set<vertex_t> > OP; // obverse paths
	{
            set<vertex_t> empty;
	    for(path_t::const_iterator ip=path.begin();ip!=path.end();++ip) {
		OP[*ip] = empty;
	    }

	    mcolor_t Q = MBG.mularcs(*(++path.begin()))[*path.begin()];
	    if( member(DiColor,Q) ) {
		get_obverse_paths(OP,Q);
	    }
	    else {
		get_obverse_paths(OP,CColor(Q));
	    }
	}
        */

	if( *path.begin()==Infty && *path.rbegin()==Infty ) {
	    outlog << "... affecting two chromosome ends" << endl;
	}
	else if( *path.begin()==Infty || *path.rbegin()==Infty ) {
	    outlog << "... affecting a chromosome end" << endl;
	}

       if( *path.begin()==*path.rbegin() ) {
	    if( path.size()%2==0 ) {
		if( *path.begin()!=Infty ) {
		    outlog << "ERROR: Semi-cycle w/o infinity!" << endl;
		    exit(1);
		}
		if( member(MBG.DiColor,MBG.mularcs(*(++path.begin()))[*path.begin()]) ) {
		    outlog << "... semi-cycle, fusion applied" << endl;
		    const string x0 = *(++path.begin());
		    const string y0 = *(++path.rbegin());
    
		    TwoBreak t(Infty,x0,Infty,y0,MBG.mularcs(x0)[Infty]);
		    if( true ) {
			if( t.apply(MBG,true) ) {
			    path.erase(--path.end());
			    *path.begin() = y0;
			    nr++;
			}
		    }
		    else {
			outlog << "... OOPS! Circularity detected! Fusion reverted, semi-cycle truncated" << endl;
			path.erase(--path.end());
		    }
		}
		else {
		    outlog << "... semi-cycle, fission applied" << endl;
		    const string y0 = *(++path.rbegin());
		    const string y1 = *(++++path.rbegin());

		    if( TwoBreak(y0,y1,Infty,Infty,MBG.mularcs(y0)[y1]).apply(MBG,true) ) {
			nr++;

			path.erase(--path.end());
			*path.rbegin() = Infty;
		    }
		}
		if( path.size() < 4 ) return nr;
	    }
	    else outlog << "... cycle" << endl;
	}

        mcolor_t Q;

	while( true ) {
	    // multicolor of (z1,z2). N.B.: x2 is NOT oo
	    outlog << "... multicolors of first and second multiedges: ";
    
	    Q = MBG.mularcs(*(++path.begin()))[*path.begin()];
    
	    outlog << MBG.set2name(Q) << " + " << MBG.set2name(MBG.CColor(Q)) << endl;

	    if( member(MBG.DiColor,Q) ) break;

	    if( *path.begin()==*path.rbegin() ) {
		outlog << "... rotating" << endl;
		path.push_back(*path.begin());
		path.erase(path.begin());
	    }
	    else {
		if( *path.begin()==Infty && *path.rbegin()!=Infty ) {
		    outlog << "... flipping" << endl;
		    for(path_t::iterator ip=++path.begin();ip!=path.end();) {
			path.push_front(*ip);
			path.erase(ip++);
		    }
		}
		if( *path.rbegin()==Infty ) {
		    outlog << "... extending beyond oo" << endl;
		    path.push_back(Infty);
		    path.erase(path.begin());
		}
		else {
		    outlog << "... truncating ??" << endl;
		    path.erase(path.begin());
		    path.erase(--path.end());
		    if( path.size()<4 ) return nr;
		}
	    }
	}

	mcolor_t Qrep = MBG.CColorRep(Q);

	// x1 -- x2 -- x3 -- ... -- x2k
	// results in:
	// x2 == x3   x4 == x5 ... x(2k-2) == x(2k-1)   x1 -- x2k

	path_t::const_iterator z3 = path.begin();
	path_t::const_iterator z0 = z3++;
	path_t::const_iterator z1 = z3++;
	path_t::const_iterator z2 = z3++;

	while(z3!=path.end()) {
	    if( TwoBreak(*z0,*z1,*z3,*z2,Q).apply(MBG,true) ) {
        	nr++;
	    }
	    else {
		z0 = z2;
	    }
	    z1 = z3++;
    	    if( z3==path.end() ) break;
	    z2 = z3++;
	}
	
	//if( *path.begin()!=*path.rbegin() && *path.begin()!=Infty && *path.rbegin()!=Infty ) {
	//    for(set<int>::const_iterator iq=Q.begin();iq!=Q.end();++iq) {
	//	LG[*iq].insert(*path.begin(),*path.rbegin());
	//    }
	//}


	//if( *path.begin()!=*path.rbegin() && *path.begin()!=Infty && *path.rbegin()!=Infty ) {
	//    EdgeW[make_pair(Qrep,edge2_t(*path.begin(),*path.rbegin()))] += PathLength(path,Qrep);
	//}

	//nr += (path.size()-2)/2;
	outlog << "... resolved with " << nr << " 2-breaks" << endl;
      	return nr;
    }
    return 0;
}



















// Given a set of partial BP-graphs (genomes), convolve them
// IG gives a list of genomes for which we need to reconstruct a common ancestor
// return true if LG was modified
bool ConvPhylTreeAll(MBGraph& MBG, int Stage) {

    /*
    // create set of vertices
//    const set<string>& S = B.CG; // without "h" and "t"
    set<string> S;
    for(partgraph_t::const_iterator ig=LG[0].begin();ig!=LG[0].end();++ig) {
        S.insert(ig->first);
        S.insert(ig->second);
    }
    */

    size_t nr, nf;  // number of rearrangements, fussions/fissions
    bool simplified = false;

    // simplifying graph
    outlog << "Stage: " << Stage << endl;
//    outlog << "Simplifying BPG-graphs:" << endl;
    do {
	nr = nf = 0;

	// cut the graph into connected components
	if( Stage==222 ) {

	    outlog << "Stage 222: splitting into connected components" << endl;

            // go over all T-consistent multicolors
	    for(set<mcolor_t>::const_iterator ic=MBG.DiColor.begin();ic!=MBG.DiColor.end();++ic) {
		if( ic->size()==0 || ic->size()==MBG.NumGen() ) continue; // except empty and complete multicolor
		
		const mcolor_t& Q = *ic;
		


		bool repeat = true;
		while(repeat) {
		    repeat = false;


		    equivalence<string> CC;     // connected components
		    map<string,string> QQ;      // multiedges of colors !Q
		    
		    for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
    
			const string& x = *is;
    
			mularcs_t M = MBG.mularcs(x);

                        if( M.size()==1 ) continue; // ignore complete multiedges

			for(mularcs_t::const_iterator im=M.begin();im!=M.end();++im) {
			    //if( im->first == Infty ) continue;
    
			    // edges of color Q
			    if( im->second == Q ) {
				//if( !QQ.defined(x) ) QQ.insert(x,im->first);
				QQ[x] = im->first;
			    }
			    else if( im->first!=Infty ) {  // reg. edges of color !Q
				CC.addrel(x,im->first);
			    }
			}
		    }
		    
		    // N.B. for regular edges (x,y) of m-c Q, QQ[x]=y and QQ[y]=x
		    // for irregular edge (x,oo), QQ[x] = oo
		    
		    CC.update();
		    outlog << MBG.set2name(*ic) << " ~ " << CC.classes() << endl;
    
		    typedef string concom_t;
		    map< concom_t, set< pair<string,string> > > EC; // reg. edges between diff. connected components of color Q
		    map< concom_t, set< pair<string,string> > > EI; // irreg. edges of color Q
    
		    for(map<string,string>::const_iterator iq=QQ.begin();iq!=QQ.end();++iq) {
			//if(iq->first>iq->second) continue;
			if( iq->second == Infty ) {
			    EI[ CC[iq->first] ].insert( make_pair(iq->first,iq->second) );
			    continue;
			}
			
			if( CC.isequiv(iq->first,iq->second) ) continue;
			
			EC[ CC[iq->first] ].insert( make_pair(iq->first,iq->second) );
			//EC[ CC[iq->second] ].insert( make_pair(iq->second,iq->first) );
		    }
		
		


		    map< string, pair<string,string> > FE; // free end

                    set< concom_t > processed;

                    // reg. edges between diff. connected components of color Q
		    for(map< concom_t, set< pair<string,string> > >::const_iterator ie=EC.begin();ie!=EC.end();++ie) {

                        if( member(processed, ie->first) ) continue;

                        // connected component with a single external edge
			if( ie->second.size()==1 ) {

			    pair<string,string> p = *(ie->second.begin());
			    pair<string,string> q;

                            // look for irregular edges
			    bool found = false;
    
			    for(set< pair<string,string> >::const_iterator ii=EI[ ie->first ].begin();ii!=EI[ ie->first ].end();++ii) {
				const pair<string,string>& t = *ii; // edge

                                // let check what would happen with (be-)edge e=(p.first,t.first)

				mcolor_t T = Q;
				for(int i=0;i<MBG.NumGen();++i) if( MBG.LG[i].defined(p.first) && MBG.LG[i][p.first]==t.first ) {
				    T.insert(i);
				}
                                // if e is enriched to T-consistent color, great!
				if( T.size()>Q.size() && member(MBG.Color,T) ) {
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
    
    
			    if( ! TwoBreak(p,q,Q).apply(MBG,true) ) continue;
			    nr++;
			    
			    outlog << "Stage 222.1: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second;

			    outlog << endl;

			    processed.insert( CC[p.first] );
			    if( q.first!=Infty ) processed.insert( CC[q.first] );
			    processed.insert( CC[p.second] );
			    if( q.second!=Infty ) processed.insert( CC[q.second] );

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

	if( Stage == 2222 ) {

            // search for 4-cycles

	    for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
		const string& x = *is;
		mularcs_t Mx = MBG.mularcs(x);

                bool next = false;

		for(mularcs_t::const_iterator im = Mx.begin();im!=Mx.end();++im) {

		    const string& y = im->first;
		    const mcolor_t& Q = im->second;

		    if( !member(MBG.DiColor,Q) || y==Infty ) continue;

		    mularcs_t My = MBG.mularcs(y);
		    My.erase(x);

		    for(mularcs_t::const_iterator jm = My.begin();jm!=My.end();++jm) {
			const string& z = jm->first;
			if( z==Infty ) continue;
			mulcols_t Cz = MBG.mulcols(z);

			if( member(Cz,Q) && member(Mx,Cz[Q]) ) {

			    pair< string, string > p = make_pair(x,y);
			    pair< string, string > q = make_pair(Cz[Q],z);

			    outlog << "Stage 2222: " << p.first << " - " << p.second << "\tX\t" << q.first << " - " << q.second << endl;

			    if( TwoBreak(p,q,Q).apply(MBG,true) ) nr++;

                            next = true;

			    break;

			}
		    }
                    if(next) break;
		}
	    }
	}


//        if( Stage != 222 )
	// Stage 1: loop over vertices
	for(set<string>::const_iterator is=MBG.VertexSet.begin();is!=MBG.VertexSet.end();++is) {
    
	    const string& x = *is;
	    mularcs_t M = MBG.mularcs(x);
	    
	    if( (Stage==1 || Stage==3) && M.size()==2 ) {

		list<string> path;
		path.push_back(x);
		
		set<string> processed;
		processed.insert(x);
		processed.insert(Infty);   // we count oo as already processed
		
	        for(map< string, set<int> >::const_iterator im = M.begin();im!=M.end();++im) {
		    const string& y = im->first;

    		    if( Stage==1 && !member(MBG.Color,im->second) ) continue; // not T-consistent

		    bool pathff = ( im==M.begin() );

		    string z0 = x;
		    string z = im->first;
		    
		    while( 1 ) {

			if( pathff ) path.push_front(z);
			else path.push_back(z);

			if( member(processed,z) ) break;
			processed.insert(z);
		
			mularcs_t Mz = MBG.mularcs(z);
			Mz.erase(z0);
		    
			if( Mz.size()!=1 || (Stage==1 && !member(MBG.Color,Mz.begin()->second)) ) break;

			z0 = z;
			z = Mz.begin()->first;
		    }
		    if( z==x ) break; // got a cycle from x to x, cannot extend it 
		    
		}

		nr += ProcessSimplePath(path);
	    }





	    // generalized reliable simple path
	    if( Stage==22 && M.size()>=2 ) {
		for(mularcs_t::const_iterator im = M.begin();im!=M.end();++im) {
		    const string& y = im->first;
		    const set<int>& Q = im->second;

                    if( y==Infty ) continue;
                    
		    mularcs_t My = MBG.mularcs(y);
		    My.erase(x);
		    map<string, set<int> > Mx = M;
		    Mx.erase(y);

		    if( member(Mx,Infty) && member(My,Infty) && Mx[Infty]==My[Infty] && member(MBG.Color,Mx[Infty]) ) {
			set<int> C;
			set_union(Q.begin(),Q.end(),Mx[Infty].begin(),Mx[Infty].end(),inserter(C,C.begin()));
			if( !member(MBG.Color,C) ) continue;
			for(set<int>::const_iterator iq=Mx[Infty].begin();iq!=Mx[Infty].end();++iq) {
			    MBG.LG[*iq].insert(x,y);
			}
			outlog << "Stage 22: fusion " << x << " + " << y << endl;
			nf++;
			break;
		    }
		    
		    
		    

		    mularcs_t::const_iterator Cx = Mx.end();
		    for(mularcs_t::const_iterator jm = Mx.begin();jm!=Mx.end();++jm) {
			if( !member(MBG.Color,jm->second) ) continue;
			set<int> C;
			set_union(Q.begin(),Q.end(),jm->second.begin(),jm->second.end(),inserter(C,C.begin()));
			if( member(MBG.Color,C) ) {
			    if( Cx!=Mx.end() ) { Cx=Mx.end(); break; }
			    Cx = jm;
			}
	    	    }
	    	    if( Cx == Mx.end() ) continue;
		    
		    mularcs_t::const_iterator Cy = My.end();
		    for(mularcs_t::const_iterator jm = My.begin();jm!=My.end();++jm) {
			if( !member(MBG.Color,jm->second) ) continue;
			set<int> C;
			set_union(Q.begin(),Q.end(),jm->second.begin(),jm->second.end(),inserter(C,C.begin()));
			if( member(MBG.Color,C) ) {
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

		    


	    if( Stage==2 ) {
		mulcols_t Cx = MBG.mulcols(x);

/*
		bool xgood = true;
		for(map< string, set<int> >::const_iterator im = M.begin();im!=M.end();++im) {
		    if( !member(Color, im->second) ) xgood = false;
		}
		if( !xgood ) continue;
*/

//if( M.size()==3 )
		for(map< string, set<int> >::const_iterator im = M.begin();im!=M.end();++im) {
		    const string& y = im->first;

		    const set<int>& Q = im->second; // color of central edge

                    if( y==Infty || Q.size()==MBG.NumGen() ) continue;

		    mulcols_t Cy = MBG.mulcols(y);

//if( MBG.mularcs(y).size()!=3 ) continue;

#ifdef VERBOSE
                    outlog << "Testing mobility of edge " << x << "-" << y << " " << MBG.set2name(Q) << " ";
#endif

		    // test "mobility" of central edge
		    // can it be ever find neighboring edge of the same multicolor
		    bool mobilQ = false;
		    for(mulcols_t::const_iterator jc = Cx.begin();jc!=Cx.end();++jc) {
			if( jc->second != y ) continue; // not a cental sub-edge
			const set<int>& QQ = jc->first; // color of central sub-edge (QQ is sub-multicolor of Q)
                        if( !member(MBG.DiColor,QQ) ) continue;


			for(mulcols_t::const_iterator ix=Cx.begin();ix!=Cx.end();++ix) {
			    if( ix->second==y /* || ix->second==Infty */ ) continue;
			    if( canformQ(ix->second,QQ) ) {
#ifdef VERBOSE
				outlog << "MOBIL: " << x << "-" << ix->second << " canForm: " << MBG.set2name(QQ) << endl;
#endif
				mobilQ = true;
				break;
			    }
			}
			if( mobilQ ) break;
    
			for(mulcols_t::const_iterator iy=Cy.begin();iy!=Cy.end();++iy) {
			    if( iy->second==x /* || iy->second==Infty */) continue;
			    if( canformQ(iy->second,QQ) ) {
#ifdef VERBOSE
				outlog << "MOBIL: " << y << "-" << iy->second << " canForm: " << MBG.set2name(QQ) << endl;
#endif
				mobilQ = true;
				break;
			    }
			}
			if( mobilQ ) break;
		    }
		    if( mobilQ ) continue;

#ifdef VERBOSE
		    outlog << "NOT MOBIL" << endl;
#endif

		    bool found = false;

		    for(mulcols_t::const_iterator ix=Cx.begin();ix!=Cx.end();++ix) {

			if( ix->second==y ) continue;
			const set<int>& QQ = ix->first;

#ifdef VERBOSE
                        outlog << " Sub-multiedge " << set2name(QQ) << endl;
#endif

			if( !member(MBG.DiColor,QQ) || !member(Cy,QQ) ) continue;

                        /*
			mcolor_t C;
                        set_union(Q.begin(),Q.end(),QQ.begin(),QQ.end(),inserter(C,C.begin()));

			// TODO: graph on central edges

			//if( !member(Color,C) ) continue; // do not create T-consistent color
                        */

			if( TwoBreak(x,ix->second,y,Cy[QQ],QQ).apply(MBG,true) ) {
                            found = true;
			    nf++;
			}
		    }

                    if( found ) break;
		}
	    }

            // cut hanging free ends
	    if( Stage==26 && M.size()==2 && member(M,Infty) ) {
		const mcolor_t& Q1 = M[Infty];
                vertex_t y;
		mcolor_t Q2;
		if( M.begin()->first==Infty ) {
                    y = M.rbegin()->first;
		    Q2 = M.rbegin()->second;
		}
		else {
                    y = M.begin()->first;
		    Q2 = M.begin()->second;
		}
		if( !member(MBG.DiColor,Q1) && member(MBG.DiColor,Q2) /* && !member(mularcs(y),Infty) */ ) {
		//if( member(DiColor,Q2) && !member(mularcs(y),Infty) ) {
		    outlog << "Unhanging fission:" << endl;
		    if( TwoBreak(x,y,Infty,Infty,Q2).apply(MBG,true) ) nf++;
		}
	    }


#if 0
	    if( Stage==2 && M.size()==3 ) {

/*
		bool xgood = true;
		for(map< string, set<int> >::const_iterator im = M.begin();im!=M.end();++im) {
		    if( !member(Color, im->second) ) xgood = false;
		}
		if( !xgood ) continue;
*/

		for(map< string, set<int> >::const_iterator im = M.begin();im!=M.end();++im) {
		    const string& y = im->first;

		    const set<int>& Q = im->second; // color of central edge

                    if( y==Infty ) continue;

		    mularcs_t My = MBG.mularcs(y);
		    My.erase(x);

		    map<string, set<int> > Mx = M;
		    Mx.erase(y);

		    const set<int>& Q1 = Mx.begin()->second;
		    const string& v1 = Mx.begin()->first;

		    const set<int>& Q2 = Mx.rbegin()->second;
		    const string& v2 = Mx.rbegin()->first;

		    // new processing code
		    if( member(MBG.DiColor,Q) || member(MBG.DiColorUnsafe,Q) ) {
			// test "mobility" of central edge
                        // can it be ever find close edge of the same multicolor
			bool mobilQ = false;
			for(mularcs_t::const_iterator ix=Mx.begin();ix!=Mx.end();++ix) {
			    if( ix->first!=Infty && canformQ(ix->first,Q) ) mobilQ = true;
			}
			for(mularcs_t::const_iterator ix=My.begin();ix!=My.end();++ix) {
			    if( ix->first!=Infty && canformQ(ix->first,Q) ) mobilQ = true;
			}

			if( mobilQ ) continue;
		    }

		    mulcols_t Mc = mulcols(y);

		    if( member(MBG.DiColor,Q1) && member(Mc,Q1) ) {
			if( TwoBreak(x,v1,y,Mc[Q1],Q1).apply(MBG,true) )	nf++;
		    }

		    if( member(MBG.DiColor,Q2) && member(Mc,Q2) ) {
			if( TwoBreak(x,v2,y,Mc[Q2],Q2).apply(MBG,true) ) nf++;
		    }
		    
		    
		    /*		    
		    if( !member(Color,Q1) || !member(Color,Q2) ) continue;

                    // if middle edge is T-consistent, do nothing
		    //if( member(Color,im->second) && (!member(Mx,Infty) || !member(My,Infty)) ) continue;

		    if( member(Color,Q) && (!member(Mx,Infty) || !member(My,Infty) ) && // || Mx[Infty]!=My[Infty]
			( canformQ(v1,Q) || canformQ(v2,Q) || canformQ(My.begin()->first,Q) || canformQ(My.rbegin()->first,Q) ) &&
			 (v1!=My.begin()->first || v2!=My.rbegin()->first) &&
			 (v2!=My.begin()->first || v1!=My.rbegin()->first)
			) {
			continue;
		    }

		    switch( My.size() ) {
		    case 1: 
		    {
			const set<int>& Q = My.begin()->second;
		        const string& z = My.begin()->first;
			
			// if z is good, don't untangle
			{
			    mularcs_t Mz = MBG.mularcs(z);
			    if( Mz.size() == 2 && member(Color,Q) ) {
				outlog << "possibility of 2-break detected" << endl;
				continue; // save x,y,z for a 2-break
			    }
			}
			
			TwoBreak(x,v1,y,z,Q1).apply(MBG);
			TwoBreak(x,v2,y,z,Q2).apply(MBG);

			break;
		    }
		    case 2:
		    {
                        string u1, u2;
			if( My.begin()->second == Q1 ) {
                            u1 = My.begin()->first;
			    u2 = My.rbegin()->first;
			}
                        else if( My.begin()->second == Q2 ) {
                            u1 = My.rbegin()->first;
			    u2 = My.begin()->first;
			}
			else continue;
			TwoBreak(x,v1,y,u1,Q1).apply(MBG);
			TwoBreak(x,v2,y,u2,Q2).apply(MBG);
                        break;
		    }
		    default: continue;
		    }
		    nf += 2;
		    */
		    
		    
		    break;
		}
	    }
#endif

	}
	outlog << nr << " 2-breaks performed" << endl;
	outlog << nf << " untangling 2-breaks performed" << endl;
	if(nr>0||nf>0) simplified = true;
    } while(nr>0||nf>0);
    outlog << endl;

    return(simplified);
}

