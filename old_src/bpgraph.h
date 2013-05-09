/* 
** Module: Breakpoint Graphs support 
** Version: 2.2
**
** This file is part of the 
** Multiple Genome Rearrangements and Ancestors (MGRA) 
** reconstruction software. 
** 
** Copyright (C) 2008,09 by Max Alekseyev <maxal@cse.sc.edu> 
**. 
** This program is free software; you can redistribute it and/or 
** modify it under the terms of the GNU General Public License 
** as published by the Free Software Foundation; either version 2 
** of the License, or (at your option) any later version. 
**. 
** You should have received a copy of the GNU General Public License 
** along with this program; if not, see http://www.gnu.org/licenses/gpl.html 
*/

#ifndef BPGRAPH_H
#define BPGRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <map>
using namespace std;

#include "gen_cls.h"

#include "symmap.h"
#include "equiv.h"

const string TH = "th";

typedef symmap<string> partgraph_t;

typedef set<string> edge_t;

map< edge_t, double > BEC, GEC;  // breakpoint reuse on black edges, gray edges, 
map< orf_t, double > BPR;	 // and endpoints
set< orf_t > USE, USE2;		 // used verticed and used at least 2 times
symmap< orf_t > ChrEnd;

class BPGraph {

public:

    symmap<string> BE, GE, OE;  // black and gray edges
    set<orf_t> CG; // common genes as vertices

    void add_edges(symmap<string>& EE, const Genome& GG, const bool closed = false) const {

	ChrEnd.clear();

    	string f,g;
        string prev_chr;
	for(genome_t::const_iterator ig=GG.begin();ig!=GG.end();++ig) {
	    if(CG.find(ig->second)==CG.end()) continue;
	    if(ig->first.first==prev_chr) {
		if(GG.sgn(ig->second)>0) {
		    EE.insert(g,ig->second + "t");
		}
		else {
		    EE.insert(g,ig->second + "h");
		}
	    }
	    else { // new chromosome detected
		if(!f.empty()) {
		    if( closed || member(GG.ChrCircular,prev_chr) ) EE.insert(f,g);
		    ChrEnd.insert(f,g);
		}
		// f is a starting vertex of a chromosome
		if(GG.sgn(ig->second)>0) f = ig->second + "t"; else f = ig->second + "h";
		prev_chr = ig->first.first;
	    }
	    // g is rightmost vertex of the previous block
	    if(GG.sgn(ig->second)>0) g = ig->second + "h"; else g = ig->second + "t";
	}
	if(!f.empty()) {
	    if( closed || member(GG.ChrCircular,prev_chr) ) EE.insert(f,g);
	    ChrEnd.insert(f,g);
	}
    }

    map<int,int> L, Lbb, Lgg, Lbg;
    // count chains
    int bb, gg, bg, ncyc, nblk;

    BPGraph(const Genome&,const Genome&,const bool);

    void count_segments(const bool);
    
    void get_blocks(map< string, set<string> >&, const bool) const;
    
    size_t num_blocks() const { return nblk; }
//	map<int,int>::const_iterator il = L.find(2);
//	if( il==L.end() ) return CG.size();
//	else return CG.size() - il->second;
//    }

    string node2gene(const string& g) const {
	return g.substr(0,g.length()-1);
    }
};



BPGraph::BPGraph(const Genome& G1,const Genome& G2,const bool closed = false) {
    { // construct set of common genes
	set<orf_t> S1;
	for(genome_t::const_iterator ig=G1.begin();ig!=G1.end();++ig) {
	    S1.insert(ig->second);
	}
	for(genome_t::const_iterator ig=G2.begin();ig!=G2.end();++ig) {
	    if(S1.find(ig->second)!=S1.end()) {
		CG.insert(ig->second);
		OE.insert(ig->second + "t",ig->second + "h");
	    }
	}
	//cout << "Common blocks: " << CG.size() << endl;
    }

    // adding black and gray edges
    add_edges(BE,G1,closed);
    add_edges(GE,G2,closed);

    count_segments(closed);
}

void BPGraph::count_segments(const bool closed = false) {
BEC.clear(); 
GEC.clear();
BPR.clear();
USE.clear();
USE2.clear();

    L.clear();
    Lbb.clear();
    Lgg.clear();
    Lbg.clear();
    {   // count chains
        ncyc = bb = gg = bg = 0;
	set<string> processed;
	for(set<orf_t>::const_iterator ig=CG.begin();ig!=CG.end();++ig) for(int j=0;j<2;++j) {
	    string gene = *ig + TH[j];
	    if(processed.find(gene)!=processed.end()) continue;
	    if(!BE.defined(gene) && !GE.defined(gene)) {
		processed.insert(gene); 
		
		BPR[gene] += 2; USE.insert( gene );  // FIXME: what to do w/ reuses (USE2) here?

		//cout << "Isolated " << gene << endl;
		Lbg[1]++;
		bg++;
		continue;
	    }

set<orf_t> CurrG;
	    int len = 0;
	    if(!BE.defined(gene) && GE.defined(gene)) {
		string y = gene;
		while(processed.find(y)==processed.end()) {
		    processed.insert(y); CurrG.insert(y);
		    len++;
		    if(GE.defined(y)) y = GE[y];
		    else {
			bg++;
			// cout << "bg:" << len << " ";
			Lbg[len]++;
			break;
		    }
		    processed.insert(y); CurrG.insert(y);
		    len++;
		    if(BE.defined(y)) y = BE[y];
		    else {
			gg++;
			// cout << "gg:" << len << " ";
			Lgg[len]++;
			break;
		    }
		}
	    }
	    if(BE.defined(gene) && !GE.defined(gene)) {
		string y = gene;
		while(processed.find(y)==processed.end()) {
		    processed.insert(y); CurrG.insert(y);
		    len++;
		    if(BE.defined(y)) y = BE[y];
		    else {
			bg++;
			// cout << "bg:" << len << " ";
			Lbg[len]++;
			break;
		    }
		    processed.insert(y); CurrG.insert(y);
		    len++;
		    if(GE.defined(y)) y = GE[y];
		    else {
			bb++;
			// cout << "bb:" << len << " ";
			Lbb[len]++;
			break;
		    }
		}
	    }
	    for(set<orf_t>::const_iterator ic=CurrG.begin();ic!=CurrG.end();++ic) {
		BPR[*ic] += (len-4)/(2.*len) + 2;
		USE.insert( *ic );
		if(len>4) USE2.insert(*ic);
	    }	    
	}		    
#if 1

	for(set<orf_t>::const_iterator ig=CG.begin();ig!=CG.end();++ig) for(int j=0;j<2;++j) {
	    string y = *ig + TH[j];
	    if(processed.find(y)!=processed.end()) continue;
	    ncyc++;
	    int len = 0;

	    set< edge_t > BES, GES;

    	    while(processed.find(y)==processed.end()) {
		
    edge_t e1; e1.insert(y); e1.insert(BE[y]);
    BES.insert(e1);
		
	        processed.insert(y);
		y = BE[y];
	
    edge_t e2; e2.insert(y); e2.insert(GE[y]);
    GES.insert(e2);

		processed.insert(y);
		y = GE[y];
		len+=2;
	    }
	    //cout << "cyc:" << len << " ";
	    L[len]++;

	      for(set<edge_t>::const_iterator is=BES.begin();is!=BES.end();++is) {
		if( len>4 ) {
	    	  BEC[*is] = (len-4)/(double)len;
		  BPR[*(is->begin())] += (len-4)/(2.*len);
		  BPR[*(is->rbegin())] += (len-4)/(2.*len);
		}
		if( len>=4 ) {
		    USE.insert( *(is->begin()) );
		    USE.insert( *(is->rbegin()) );
		}    
		if( len>4 ) {
		    USE2.insert( *(is->begin()) );
		    USE2.insert( *(is->rbegin()) );
		}    
	      }
	      for(set<edge_t>::const_iterator is=GES.begin();is!=GES.end();++is) {
		if( len>4 ) {
	          GEC[*is] = (len-4)/(double)len;
		  BPR[*(is->begin())] += (len-4)/(2.*len);
		  BPR[*(is->rbegin())] += (len-4)/(2.*len);
	        }		
		if( len>=4 ) {
		    USE.insert( *(is->begin()) );
		    USE.insert( *(is->rbegin()) );
		}    
		if( len>4 ) {
		    USE2.insert( *(is->begin()) );
		    USE2.insert( *(is->rbegin()) );
		}    
	      }


	}
	//cout << endl;
	
	//cout << "# bg/bb/gg: " << bg << " " << bb << " " << gg << endl;
	//cout << "# cycles: " << ncyc << "\ttrivial: " << L[2] << endl;
	
	nblk = CG.size() - L[2];
	ncyc -= L[2];
	L.erase(2);
	
	//cout << "# blocks: " << nblk << endl;
	//cout << endl;
	
	int tot = 0, tL = 0, ncyc_odd = 0;

	//cout << "cycles: ";
	for(map<int,int>::const_iterator il=L.begin();il!=L.end();++il) {
	    //cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    tL += il->second;
	    if((il->first/2)%2) ncyc_odd += il->second;
	}
	//cout << "(total vertices: " << tot << " " << "cycles: " << tL << ")" << endl;
	//cout << "bodd/godd: " << ncyc_odd << "\tbeven/geven: " << tL - ncyc_odd << endl;
	//cout << endl;

if( !closed ) {
	cout << "bb-chains: ";
	tot = 0; 
	int nlbb = 0, nlbb_bodd = 0;
	for(map<int,int>::const_iterator il=Lbb.begin();il!=Lbb.end();++il) {
	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlbb += il->second;
	    if( (il->first/2)%2 ) nlbb_bodd += il->second;
	}
	cout << "(total vertices: " << tot << " " << "chains: " << nlbb << ")" << endl;
	cout << "bodd/geven: " << nlbb_bodd << "\tbeven/godd: " << nlbb - nlbb_bodd << endl;
	cout << endl;
	
	cout << "gg-chains: ";
	tot = 0;
	int nlgg = 0, nlgg_bodd = 0;
	for(map<int,int>::const_iterator il=Lgg.begin();il!=Lgg.end();++il) {
	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlgg += il->second;
	    if( !((il->first/2)%2) ) nlgg_bodd += il->second;
	}
	cout << "(total vertices: " << tot << " " << "chains: " << nlgg << ")" << endl;
	cout << "bodd/geven: " << nlgg_bodd << "\tbeven/godd: " << nlgg - nlgg_bodd << endl;
	cout << endl;

	cout << "bg-chains: ";
	tot = 0;
	int nlbg = 0, nlbg_bodd = 0;
	for(map<int,int>::const_iterator il=Lbg.begin();il!=Lbg.end();++il) {
	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlbg += il->second;
	    if( ((il->first-1)/2)%2 ) nlbg_bodd += il->second;
	}
	cout << "(total vertices: " << tot << " " << "chains: " << nlbg << ")" << endl;
	cout << "bodd/godd: " << nlbg_bodd << "\tbeven/geven: " << nlbg - nlbg_bodd << endl;
	cout << endl;
	
	int d1 = (nlgg - nlgg_bodd) - (2*nlbg_bodd - nlbg)/2;
	int d2 = nlbb_bodd - (2*nlbg_bodd - nlbg)/2;
	if( d1>0 ) d1 %= 2; else d1 = 0;
	if( d2>0 ) d2 %= 2; else d2 = 0;

	cout << "B_2: " << nblk - ncyc - max(1,nlbg/2) - nlbb << "\t" << nblk - ncyc - max(1,nlbg/2) - nlgg << endl;

	int L3PQ = ncyc_odd + nlbb_bodd + d1 + max(0,abs(2*nlbg_bodd - nlbg)/2 - (nlgg - nlgg_bodd));
	int L3QP = ncyc_odd + (nlgg - nlgg_bodd) + d2 + max(0,abs(2*nlbg_bodd - nlbg)/2 - nlbb_bodd);
	
	cout << "B_3: " << (nblk - L3PQ)/2 << "\t" << (nblk - L3QP)/2 << endl;
	cout << "B_23: " << (3*nblk - 2*nlbb - 2*ncyc - nlbg - L3PQ)/2 << "\t" 
	<< (3*nblk - 2*nlgg - 2*ncyc - nlbg - L3QP)/2 << endl;
	cout << endl;
}	
	cout << endl;
#endif
    }		    
}

// compute conserved segments, reversing single genes if needed
void BPGraph::get_blocks(map< string, set<string> >& C, const bool GeneRev = false) const {
    equivalence<string> BLK;

    for(symmap<string>::const_iterator io=OE.begin();io!=OE.end();++io) {
	BLK.insert(node2gene(io->first));
    }
    for(symmap<string>::const_iterator ib=BE.begin();ib!=BE.end();++ib) {
	if(GE.defined(ib->first) && GE[ib->first]==ib->second) {
	    BLK.addrel(node2gene(ib->first),node2gene(ib->second));
	}
	if( GeneRev ) {
	    const string& t = ib->first;
	    const string& h = OE[ib->first];
	    if(GE.defined(ib->second) && GE[ib->second]==h
		&& GE.defined(t) && BE.defined(GE[t]) && BE[GE[t]]==h)
	    {
		BLK.addrel(node2gene(t),node2gene(ib->second));
		BLK.addrel(node2gene(t),node2gene(GE[t]));
	    }
	}
    }
    cout << "Classes: " << BLK.classes() << endl;
    BLK.get_eclasses(C);
}

    














// returns: # trivial cycles, # (true) synteny blocks, d2, d3
vector<size_t> genome_dist(const partgraph_t& BE, const partgraph_t& GE, const partgraph_t& OE, const bool closed = false) {

vector<size_t> res(4);

BEC.clear(); 
GEC.clear();
BPR.clear();
USE.clear();
USE2.clear();

/*
    L.clear();
    Lbb.clear();
    Lgg.clear();
    Lbg.clear();
*/
    map<int,int> L, Lbb, Lgg, Lbg;
    int  ncyc, bb, gg, bg;
    {   // count chains
        ncyc = bb = gg = bg = 0;
	set<string> processed;
	for(partgraph_t::const_iterator oi=OE.begin();oi!=OE.end();++oi) for(int j=0;j<2;++j) {
	    const string& gene = (j==0)?(oi->first):(oi->second);

	    if(processed.find(gene)!=processed.end()) continue;
	    if(!BE.defined(gene) && !GE.defined(gene)) {
		processed.insert(gene); 
/*		
		BPR[gene] += 2; USE.insert( gene );
*/
		//cout << "Isolated " << gene << endl;
		Lbg[1]++;
		bg++;
		continue;
	    }

//set<orf_t> CurrG;
	    int len = 0;
	    if(!BE.defined(gene) && GE.defined(gene)) {
		string y = gene;
		while(processed.find(y)==processed.end()) {
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if(GE.defined(y)) y = GE[y];
		    else {
			bg++;
			// cout << "bg:" << len << " ";
			Lbg[len]++;
			break;
		    }
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if(BE.defined(y)) y = BE[y];
		    else {
			gg++;
			// cout << "gg:" << len << " ";
			Lgg[len]++;
			break;
		    }
		}
	    }
	    if(BE.defined(gene) && !GE.defined(gene)) {
		string y = gene;
		while(processed.find(y)==processed.end()) {
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if(BE.defined(y)) y = BE[y];
		    else {
			bg++;
			// cout << "bg:" << len << " ";
			Lbg[len]++;
			break;
		    }
		    processed.insert(y); //CurrG.insert(y);
		    len++;
		    if(GE.defined(y)) y = GE[y];
		    else {
			bb++;
			// cout << "bb:" << len << " ";
			Lbb[len]++;
			break;
		    }
		}
	    }
/*
	    for(set<orf_t>::const_iterator ic=CurrG.begin();ic!=CurrG.end();++ic) {
		BPR[*ic] += (len-4)/(2.*len) + 2;
		USE.insert( *ic );
	    }	    
*/
	}		    
#if 1

	for(partgraph_t::const_iterator oi=OE.begin();oi!=OE.end();++oi) for(int j=0;j<2;++j) {
	    string y = (j==0)?(oi->first):(oi->second);

	    if(processed.find(y)!=processed.end()) continue;
	    ncyc++;
	    int len = 0;

    set< edge_t > BES, GES;

    	    while(processed.find(y)==processed.end()) {
		
    edge_t e1; e1.insert(y); e1.insert(BE[y]);
    BES.insert(e1);
		
	        processed.insert(y);
		y = BE[y];
	
    edge_t e2; e2.insert(y); e2.insert(GE[y]);
    GES.insert(e2);

		processed.insert(y);
		y = GE[y];
		len+=2;
	    }
	    //cout << "cyc:" << len << " ";
	    L[len]++;

	      for(set<edge_t>::const_iterator is=BES.begin();is!=BES.end();++is) {
		if( len>4 ) {
	    	  BEC[*is] = (len-4)/(double)len;
		  BPR[*(is->begin())] += (len-4)/(2.*len);
		  BPR[*(is->rbegin())] += (len-4)/(2.*len);
		}
		if( len>=4 ) {
		    USE.insert( *(is->begin()) );
		    USE.insert( *(is->rbegin()) );
		}    
		if( len>4 ) {
		    USE2.insert( *(is->begin()) );
		    USE2.insert( *(is->rbegin()) );
		}    
	      }
	      for(set<edge_t>::const_iterator is=GES.begin();is!=GES.end();++is) {
		if( len>4 ) {
	          GEC[*is] = (len-4)/(double)len;
		  BPR[*(is->begin())] += (len-4)/(2.*len);
		  BPR[*(is->rbegin())] += (len-4)/(2.*len);
	        }		
		if( len>=4 ) {
		    USE.insert( *(is->begin()) );
		    USE.insert( *(is->rbegin()) );
		}    
		if( len>4 ) {
		    USE2.insert( *(is->begin()) );
		    USE2.insert( *(is->rbegin()) );
		}    
	      }


	}
	//cout << endl;
//	cout << "# bg/bb/gg: " << bg << " " << bb << " " << gg << endl;
//	cout << "# cycles: " << ncyc << "\ttrivial: " << L[2] << endl;

res[0] = L[2];
    	
	int nblk = OE.size() - L[2];
	ncyc -= L[2];
	L.erase(2);

res[1] = nblk;
	
//	cout << "# blocks: " << nblk << endl;
//	cout << endl;
	
	int tot = 0, tL = 0, ncyc_odd = 0;

//	cout << "cycles: ";
	for(map<int,int>::const_iterator il=L.begin();il!=L.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    tL += il->second;
	    if((il->first/2)%2) ncyc_odd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "cycles: " << tL << ")" << endl;
//	cout << "bodd/godd: " << ncyc_odd << "\tbeven/geven: " << tL - ncyc_odd << endl;
//	cout << endl;
if( closed ) {
    res[2] = nblk - tL;
    res[3] = (nblk - ncyc_odd) / 2;
}

if( !closed ) {
//	cout << "bb-chains: ";
	tot = 0; 
	int nlbb = 0, nlbb_bodd = 0;
	for(map<int,int>::const_iterator il=Lbb.begin();il!=Lbb.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlbb += il->second;
	    if( (il->first/2)%2 ) nlbb_bodd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "chains: " << nlbb << ")" << endl;
//	cout << "bodd/geven: " << nlbb_bodd << "\tbeven/godd: " << nlbb - nlbb_bodd << endl;
//	cout << endl;
	
//	cout << "gg-chains: ";
	tot = 0;
	int nlgg = 0, nlgg_bodd = 0;
	for(map<int,int>::const_iterator il=Lgg.begin();il!=Lgg.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlgg += il->second;
	    if( !((il->first/2)%2) ) nlgg_bodd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "chains: " << nlgg << ")" << endl;
//	cout << "bodd/geven: " << nlgg_bodd << "\tbeven/godd: " << nlgg - nlgg_bodd << endl;
//	cout << endl;

//	cout << "bg-chains: ";
	tot = 0;
	int nlbg = 0, nlbg_bodd = 0;
	for(map<int,int>::const_iterator il=Lbg.begin();il!=Lbg.end();++il) {
//	    cout << il->first << "^" << il->second << " ";
	    tot += il->first * il->second;
	    nlbg += il->second;
	    if( ((il->first-1)/2)%2 ) nlbg_bodd += il->second;
	}
//	cout << "(total vertices: " << tot << " " << "chains: " << nlbg << ")" << endl;
//	cout << "bodd/godd: " << nlbg_bodd << "\tbeven/geven: " << nlbg - nlbg_bodd << endl;
//	cout << endl;
	
	int d1 = (nlgg - nlgg_bodd) - (2*nlbg_bodd - nlbg)/2;
	int d2 = nlbb_bodd - (2*nlbg_bodd - nlbg)/2;
	if( d1>0 ) d1 %= 2; else d1 = 0;
	if( d2>0 ) d2 %= 2; else d2 = 0;

//	cout << "B_2: " << nblk - ncyc - max(1,nlbg/2) - nlbb << "\t" << nblk - ncyc - max(1,nlbg/2) - nlgg << endl;

	res[2] = max( nblk - ncyc - max(1,nlbg/2) - nlbb, nblk - ncyc - max(1,nlbg/2) - nlgg );

	int L3PQ = ncyc_odd + nlbb_bodd + d1 + max(0,abs(2*nlbg_bodd - nlbg)/2 - (nlgg - nlgg_bodd));
	int L3QP = ncyc_odd + (nlgg - nlgg_bodd) + d2 + max(0,abs(2*nlbg_bodd - nlbg)/2 - nlbb_bodd);
	
//	cout << "B_3: " << (nblk - L3PQ)/2 << "\t" << (nblk - L3QP)/2 << endl;

	res[3] = max( (nblk - L3PQ)/2, (nblk - L3QP)/2 );

//	cout << "B_23: " << (3*nblk - 2*nlbb - 2*ncyc - nlbg - L3PQ)/2 << "\t" 
//	<< (3*nblk - 2*nlgg - 2*ncyc - nlbg - L3QP)/2 << endl;
//	cout << endl;
}	
//	cout << endl;
#endif
    }		    
    return res;
}

#endif
