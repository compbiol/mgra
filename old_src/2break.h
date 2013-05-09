/* 
** Module: 2-Breaks and Transformations support
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

#ifndef TWOBREAK_H
#define TWOBREAK_H

#include "mcolor.h"

class TwoBreak {
public:
    arc_t OldArc[2];     // (x1,y1) x (x2,y2) = (x1,x2) + (y1,y2)
    mcolor_t MultiColor;

    static list< TwoBreak > History;

    TwoBreak() {};

    TwoBreak(const arc_t& a1, const arc_t& a2, const mcolor_t& Q) {
	OldArc[0] = a1;
	OldArc[1] = a2;
	MultiColor = Q;
    }

    TwoBreak(const vertex_t& x1, const vertex_t& y1, const vertex_t& x2, const vertex_t& y2, const mcolor_t& Q) {
	OldArc[0] = make_pair(x1,y1);
	OldArc[1] = make_pair(x2,y2);
	MultiColor = Q;
    }

    bool islinear(MBGraph& M) const;

    bool apply(MBGraph& M, bool record = false) const {
	//clog << " (" << OldArc[0].first << "," << OldArc[0].second << ")x(" << OldArc[1].first << "," << OldArc[1].second << "):{" << set2name(MultiColor) << "} ";

        /*
	set<vertex_t> SS;
	SS.insert(OldArc[0].first);
	SS.insert(OldArc[0].second);
	SS.insert(OldArc[1].first);
	SS.insert(OldArc[1].second);
	if( member(SS,"774t") || member(SS,"774h") ) {
	    clog << endl << "DEBUG: (" << OldArc[0].first << "," << OldArc[0].second << ")x(" << OldArc[1].first << "," << OldArc[1].second << "):{" << set2name(MultiColor) << "} " << endl;
	}
        */

	if( record ) {
	    outlog << " " << *this << " ";
#ifdef PREVENT_CIRCULAR
	    if( !islinear(M) ) {
                outlog << "REJECTED ";
		return false;
	    }
#endif
	    History.push_back(*this);
	}

	for(mcolor_t::const_iterator ic=MultiColor.begin();ic!=MultiColor.end();++ic) {
	    for(size_t i=0;i<2;++i) {
		if(OldArc[i].first!=Infty && OldArc[i].second!=Infty) {
		    M.LG[*ic].erase(OldArc[i].first,OldArc[i].second);
		}
	    }
      	    if(OldArc[0].first!=Infty && OldArc[1].first!=Infty) {
	        M.LG[*ic].insert(OldArc[0].first,OldArc[1].first);
	    }
      	    if(OldArc[0].second!=Infty && OldArc[1].second!=Infty) {
	        M.LG[*ic].insert(OldArc[0].second,OldArc[1].second);
	    }
	}
        return true;
    }

    void applySingle(partgraph_t& SG) const {
	for(size_t i=0;i<2;++i) {
	    if(OldArc[i].first!=Infty && OldArc[i].second!=Infty) {
		SG.erase(OldArc[i].first,OldArc[i].second);
	    }
	}
	if(OldArc[0].first!=Infty && OldArc[1].first!=Infty) {
	    SG.insert(OldArc[0].first,OldArc[1].first);
	}
	if(OldArc[0].second!=Infty && OldArc[1].second!=Infty) {
	    SG.insert(OldArc[0].second,OldArc[1].second);
	}
    }

    inline void revert(MBGraph& M) const {
	TwoBreak(OldArc[0].first,OldArc[1].first,OldArc[0].second,OldArc[1].second,MultiColor).apply(M);
    }

    inline void revertSingle(partgraph_t& SG) const {
	TwoBreak(OldArc[0].first,OldArc[1].first,OldArc[0].second,OldArc[1].second,MultiColor).applySingle(SG);
    }

    inline TwoBreak inverse() const {
	return TwoBreak(OldArc[0].first,OldArc[1].first,OldArc[0].second,OldArc[1].second,MultiColor);
    }

    // make OldArc[0].first the smallest possible
    void normalize() {
	while( true ) {
	    if( OldArc[0].first > OldArc[0].second ) {
		OldArc[0] = make_pair( OldArc[0].second, OldArc[0].first );
		OldArc[1] = make_pair( OldArc[1].second, OldArc[1].first );
                continue;
	    }
	    if( OldArc[0].first > OldArc[1].first || OldArc[0].first > OldArc[1].second ) {
		arc_t temp = OldArc[0];
		OldArc[0] = OldArc[1];
		OldArc[1] = temp;
                continue;
	    }
	    break;
	}
    }

    friend ostream& operator << (ostream& os,const TwoBreak& t) {
        os << "(" << t.OldArc[0].first << "," << t.OldArc[0].second << ")x(" << t.OldArc[1].first << "," << t.OldArc[1].second << "):{" << MBG.set2name(t.MultiColor) << "}";
	return os;
    }


};
list< TwoBreak > TwoBreak::History;
//list< TwoBreak > TwoBreak::HistoryX;


// check if a 2-break creates a circular chromosome
bool TwoBreak::islinear(MBGraph& M) const {

    apply(M);

    for(int i=0;i<2;++i) {
	vertex_t x;
	if( i==0 ) x = OldArc[0].first;
        else x = OldArc[0].second;

        if( x==Infty ) continue;

	for(mcolor_t::const_iterator ic=MultiColor.begin();ic!=MultiColor.end();++ic) {

	    const partgraph_t& PG = M.LG[*ic];
	    bool circular = false;

	    for(string y = M.OB[x];PG.defined(y);) {
		y = PG[y];

		if( y!=x ) y = M.OB[y];

		if( y==x ) {
		    circular = true;
                    break;
		}
	    }
		
	    if( !circular && PG.defined(x) ) {
		for(string y = x;PG.defined(y);) {
		    y = PG[y];

		    if(y!=x) y = M.OB[y];

		    if(y==x) {
			circular = true;
			break;
		    }
		}
	    }

	    if( circular ) {
		revert(M);
		return false;
	    }
	}
    }
    revert(M);
    return true;
}

typedef list<TwoBreak> transform_t;

void ApplyAll(MBGraph& M, const transform_t& T, bool record = false) {
    for(transform_t::const_iterator it=T.begin();it!=T.end();++it) {
	it->apply(M,record);
    }
}

void RevertAll(MBGraph& M, const transform_t& T) {
    for(transform_t::const_reverse_iterator it=T.rbegin();it!=T.rend();++it) {
	it->revert(M);
    }
}


transform_t extractX(const transform_t& T) {
    transform_t X;
    for(transform_t::const_iterator it=T.begin();it!=T.end();++it) {
	if( member(chrX,it->OldArc[0].first) || member(chrX,it->OldArc[0].second) || member(chrX,it->OldArc[1].first) || member(chrX,it->OldArc[1].second) ) {
	    X.push_back(*it);
	}
    }
    return X;
}

/*
void readtrans(const string& fn,  transform_t& T) {


    T.clear();

    ifstream fh(fn.c_str());
    
    if(!fh) {
	cerr << "Unable to open " << fn << endl;
	exit(1);
    }

    cout << "Reading " << fn << endl;

    while(true) {

	string line;
	getline(fh,line);
	if(fh.eof()) break;

	
	if(line.empty() || line[0]=='#') continue;

	for(size_t j=0;j<line.size();++j) {
	    switch(line[j]) {
	    case '(':
	    case ',':
	    case ')':
	    case 'x':
	    case ':':
	    case '{':
	    case '}':
		line[j] = ' ';
	    }
	}

	istringstream istr(line);

	TwoBreak t;
        string MC;

	istr >> t.OldArc[0].first >> t.OldArc[0].second >> t.OldArc[1].first >> t.OldArc[1].second >> MC;
	t.MultiColor = MGB.name2set(MC);

	T.push_back(t);
    }
    
    cout << "2-Breaks read: " << T.size() << endl;
}
*/


#endif

