/* 
** Module: Multicolors and Mutiple Breakpoint Graphs support
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

#ifndef MCOLOR_H
#define MCOLOR_H

#include "bpgraph.h"

typedef string vertex_t;
typedef pair<vertex_t,vertex_t> arc_t;
typedef list<vertex_t> path_t;

//typedef vector<partgraph_t> MBGraph_t; 

typedef set<int> mcolor_t;

typedef map<vertex_t, mcolor_t > mularcs_t;
typedef map<mcolor_t, vertex_t > mulcols_t;

class MBGraph {

private:
    symmap<mcolor_t> CColorM;      // complement multicolor

    mcolor_t addtree(const string&, ofstream&);

public:
    // local graphs of each color
    vector<partgraph_t> LG;

    partgraph_t OB;         // OBverse relation
    set<string> VertexSet;  // set of vertices

    // ordered genomes
    static vector<string> num2gen;
    map<string,size_t> gen2num;    

    set< mcolor_t > Color;
    set< mcolor_t > DiColor;   // directed colors
    map< mcolor_t, int > OColor;   // dot output color
    vector< mcolor_t > TColor; // colors corresponding to ancestral genomes
    set< mcolor_t > DiColorUnsafe;

    set< pair<string,string> > LeafAdj;

    static bool SplitBadColors;

    void init(const vector<Genome>&, const list<string>&, const string&);
    
    size_t NumGen() const {
	return LG.size();
    }

    const mcolor_t& CColor(const mcolor_t& S) {
	if( !CColorM.defined(S) ) {
	    mcolor_t T;
	    for(size_t j=0; j<NumGen(); ++j) {
		if( !member(S,j) ) T.insert(j);
	    }
	    CColorM.insert(S,T);
	}
	return CColorM[S];
    }

    static string set2name(const mcolor_t&);
    mcolor_t name2set(const string&) const;

    mularcs_t mularcs(const vertex_t&) const;

    set< mcolor_t > SplitColor(const mcolor_t&) const;

    mulcols_t mulcols(const vertex_t&) const;

    friend ostream& operator << (ostream&,const mcolor_t&);

    // complementary colors representative
    mcolor_t CColorRep(const mcolor_t& c) {
	mcolor_t Q = CColor(c);
	if( Q.size()>(c).size() || (Q.size()==(c).size() && Q>(c)) ) return c;
	return Q;
    }

    bool AreAdjacentBranches(const mcolor_t&, const mcolor_t&) const;
};

bool MBGraph::SplitBadColors = false;
vector<string> MBGraph::num2gen;




///////////////////////////////////////////////////////////////////////////////


mcolor_t MBGraph::addtree(const string& t, ofstream& fleg) {
    if( t[0]=='(' ) {
	if( t[t.size()-1] != ')' ) {
	    cerr << "ERROR: Malformed input (sub)tree 1" << endl;
	    exit(3);
	}

	// non-trivial tree
	
	int p=0;
	for(size_t j=1;j<t.size() - 1; ++j) {
	    switch( t[j] ) {
	    case '(': p++; break;
	    case ')': p--; break;
	    case ',': if( p==0 )
		{

		    mcolor_t Q1 = addtree( t.substr(1,j-1), fleg );
		    mcolor_t Q2 = addtree( t.substr(j+1,t.size()-j-2), fleg );

		    DiColor.insert( Q1 );
		    DiColor.insert( Q2 );

		    mcolor_t Q = Q1;
		    copy(Q2.begin(),Q2.end(),inserter(Q,Q.end()));
		    
		    fleg << "\t\"" << set2name(Q) << "\"\t->\t\"" << set2name(Q1) << "\";" << endl;
		    fleg << "\t\"" << set2name(Q) << "\"\t->\t\"" << set2name(Q2) << "\";" << endl;

		    return Q;
		}
	    }

	    if( p<0 ) {
		cerr << "ERROR: Malformed input (sub)tree 2" << endl;
		exit(3);
	    }
	}
	if( p!=0 ) {
	    cerr << "ERROR: Malformed input (sub)tree 3" << endl;
	    exit(3);
	}
    }
    else {
	// single node
	mcolor_t Q;
	for(size_t j=0;j<t.size();++j) {
	    string c = t.substr(j,1);
	    if( !member(gen2num,c) ) {
		cerr << "ERROR: Unknown genome in (sub)tree: " << t << endl;
		exit(3);
	    }
	    Q.insert( gen2num[c] );
	}

	return Q;
    }
}



void MBGraph::init(const vector<Genome>& G, const list<string>& trees, const string& target) {

    num2gen.resize(G.size());
    LG.resize(G.size());

    for(size_t i=0;i<G.size();++i) {
	num2gen[i] = G[i].name;
        gen2num[G[i].name] = i;

	clog << "Genome " << G[i].name << " blocks: " << G[i].size() << endl;
    }

    {
	BPGraph B(G[0],G[1],true);
    
	// obverse edges
	OB = B.OE;
    
	for(partgraph_t::const_iterator ig=OB.begin();ig!=OB.end();++ig) {
	    VertexSet.insert(ig->first);
	    VertexSet.insert(ig->second);
	}
    
	size_t tcls = 0;
	for(size_t i=0;i<G.size();++i) {
	    B.add_edges(LG[i],G[i],false);
	    tcls += ChrEnd.size();
	    for(partgraph_t::const_iterator ia=LG[i].begin();ia!=LG[i].end();++ia) {
		LeafAdj.insert(*ia);
	    }
	}
	//clog << "# irregular edges: " << tcls << endl;
    }

#ifdef HAVE_BLOCK_LENGTH
    // initialize block length
    for(map<orf_t,cpan_t>::const_iterator ic=G[REFGENOME].orf2cpan.begin();ic!=G[REFGENOME].orf2cpan.end();++ic) {
	const orf_t& x = ic->first;
	BL[x+"h"] = BL[x+"t"] = ic->second.second.second - ic->second.second.first;
    }
#endif

    // init colors
    Color.clear();
    CColorM.clear();
    DiColor.clear();  DiColorUnsafe.clear(); 
    OColor.clear();
    TColor.clear();


    // parse given trees and save legend.dot
    {
	ofstream flegend("legend.dot");
	flegend << "digraph Legend {" << endl;

	flegend << "\tnode [style=filled];" << endl;

	// add terminal branches
	for( size_t j = 0; j<G.size(); ++j ) {
	    mcolor_t C;
	    C.insert(j);
	    DiColor.insert(C);
	
	    flegend << "\t\"" << set2name(C) << "\"\t[fillcolor=" <<  RGBcols[RGBcoeff * j]  << "];" << endl;
	}

	
	for(list<string>::const_iterator it=trees.begin(); it!=trees.end(); ++it) {
	    mcolor_t C = addtree(*it,flegend);
	    if( C.size() < G.size() ) DiColor.insert(C); // complete multicolor is excluded
	}
	
	flegend << "}" << endl;
	flegend.close();
    }

    if( !target.empty() ) DiColor.erase( name2set(target) );

    // check consistency
    for( set<mcolor_t>::const_iterator id=DiColor.begin(); id!=DiColor.end(); ++id ) {
	for( set<mcolor_t>::const_iterator jd=id; ++jd!=DiColor.end(); ) {
	    mcolor_t C;
	    set_intersection(id->begin(),id->end(),jd->begin(),jd->end(),inserter(C,C.begin()));
	    if( !C.empty() && C.size()!=id->size() && C.size()!=jd->size() ) {
		clog << "Multicolors " << set2name(*id) << " " << set2name(*jd) << " have nontrivial intersection, removing the latter" << endl;
		DiColor.erase(jd++);
		--jd;
	    }
	}
    }
	    
    


    TColor.resize( DiColor.size() );

    clog << "vecT-consistent colors: " << DiColor.size() << endl;

    size_t col = 1;
    for( set<mcolor_t>::const_iterator id=DiColor.begin(); id!=DiColor.end(); ++id ) {

	clog << "\t" << set2name(*id);
	
	Color.insert(*id);

	// compute complement to *id
	mcolor_t C;
	for(size_t j=0;j<G.size();++j) {
	    if( !member(*id,j) ) C.insert(j);
	}
	Color.insert(C);

	TColor[col-1] = *id;

	OColor[*id] = OColor[C] = col++;
    }
    clog << endl;

    // tell where the root resides
    clog << "the root resides in between:";
    set< mcolor_t > T = DiColor;
    for( set<mcolor_t>::const_iterator it=T.begin(); it!=T.end(); ++it ) {
	for( set<mcolor_t>::const_iterator jt=it; ++jt!=T.end(); ) {
	    mcolor_t C;
	    set_intersection(it->begin(),it->end(),jt->begin(),jt->end(),inserter(C,C.begin()));
	    if( C.size() == it->size() ) {
		T.erase(it++);
		jt = it;
		continue;
	    }
	    if( C.size() == jt->size() ) {
		T.erase(jt++);
		--jt;
	    }
	}
	clog << " " << set2name(*it);
    }
    clog << endl;
}



string MBGraph::set2name(const mcolor_t& S) {
    if( S.empty() ) return "\\ensuremath{\\emptyset}";
    ostringstream os;
    for(mcolor_t::const_iterator is=S.begin();is!=S.end();++is) {
	os << num2gen[*is];
    }
    return os.str();
}

mcolor_t MBGraph::name2set(const string& s) const {
    mcolor_t C;
    for(size_t j=0;j<s.size();++j) {
	string t = s.substr(j,1);
	if( gen2num.find(t) == gen2num.end() ) {
	    cerr << "ERROR: Malformed multicolor " << s << endl;
	    exit(1);
	}
	C.insert( gen2num.find(t)->second );
    }
    return C;
}



mularcs_t MBGraph::mularcs(const vertex_t& y) const {
    if( y==Infty ) {
	cerr << "mularcs ERROR: Infinite input" << endl;
	exit(1);
    }
    mularcs_t My;
    for(int i=0;i<NumGen();++i) {
	if(LG[i].defined(y)) My[LG[i][y]].insert(i);
	else My[Infty].insert(i);
    }
    return My;
}


set< mcolor_t > MBGraph::SplitColor(const mcolor_t& Q) const {
    set< mcolor_t > S;

    if( member(DiColor,Q) ) {
	S.insert(Q);
	return S;
    }

    if( !SplitBadColors ) {
        return S;
    }

    equivalence< int > EQ;
    for(mcolor_t::const_iterator iq=Q.begin();iq!=Q.end();++iq) EQ.addrel(*iq,*iq);

    for(set<mcolor_t>::const_iterator ic=DiColor.begin();ic!=DiColor.end();++ic) {
	mcolor_t C;
	set_intersection(ic->begin(),ic->end(),Q.begin(),Q.end(),inserter(C,C.begin()));
	if( C.size()>=2 && C.size()==ic->size() ) {
	    for(mcolor_t::const_iterator iq=C.begin();iq!=C.end();++iq) {
		EQ.addrel(*iq,*C.begin());
	    }
	}
    }

    EQ.update();
    map< int,set<int> > cls;
    EQ.get_eclasses(cls);
    for(map< int,set<int> >::const_iterator ic=cls.begin();ic!=cls.end();++ic) {
	S.insert(ic->second);
    }
    return S;
}



mulcols_t MBGraph::mulcols(const vertex_t& y) const {
    mularcs_t My = mularcs(y);
    mulcols_t R;
    for(mularcs_t::const_iterator im=My.begin();im!=My.end();++im) {
	const mcolor_t& Q = im->second;
	if( SplitBadColors && !member(DiColor,Q) && Q.size()<NumGen() ) {
            set< mcolor_t > C = SplitColor(Q);
	    for(set<mcolor_t>::const_iterator ic=C.begin();ic!=C.end();++ic) {
		R[*ic] = im->first;
	    }
	}
        else R[Q] = im->first;
    }
    return R;
}


ostream& operator << (ostream& os, const mcolor_t& C) {
    os << MBGraph::set2name(C);
    //copy( S.begin(), S.end(), ostream_iterator<int>( os, " " ) );
    return os;
}

bool MBGraph::AreAdjacentBranches(const mcolor_t& A, const mcolor_t& B) const {

    if( !member(Color,A) || !member(Color,B) ) return false;

    mcolor_t Q1, Q2;

    if( A.size() >= B.size() ) {
	Q1 = A;
	Q2 = B;
    }
    else {
	Q1 = B;
	Q2 = A;
    }

    mcolor_t C;
    set_difference(Q1.begin(),Q1.end(),Q2.begin(),Q2.end(),inserter(C,C.begin()));
    if( C.size()==Q1.size()-Q2.size() && member(Color,C) ) return true;

    C.clear();
    set_union(Q1.begin(),Q1.end(),Q2.begin(),Q2.end(),inserter(C,C.begin()));
    if( C.size()==Q1.size()+Q2.size() && member(Color,C) ) return true;

    return false;
}

#endif
