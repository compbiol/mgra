/* 
** Module: Genome generic class
** Version: 1.1
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


#ifndef GEN_CLS_H
#define GEN_CLS_H

#include <map>
#include <string>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;


#define CHR_INTERSPACE 100000

typedef string orf_t;
//typedef couple<orf_t> match_t;
typedef pair<int,int> span_t; 		// interval in absolute coordinates
typedef pair<string,int> coord_t;   	// (contig,offset)
typedef pair<string,span_t> cpan_t;  	// (chr,start,end)
typedef map<coord_t,orf_t> genome_t;
typedef multimap<orf_t,orf_t> genemap_t;

class Genome : public genome_t {
public:

    using genome_t::find;
    using genome_t::begin;
    using genome_t::end;
    

    string name;

    /* Absolute nucleotide scale coordinates */
/*
    //map<string,int> ChrLen;   // chromosome -> length
    map<string,int> ChrAbsOfs;   // chromosome -> absolute coordinate
    int AbsLen;
    map<span_t,orf_t> as2orf;    // was: Kabs, Sabs
    map<orf_t,span_t> orf2as;    // was: rKSabs
*/
    map<orf_t,cpan_t> orf2cpan;

    /* Absolute gene scale coordinates */
    map<string,int> ChrGeneOfs;   // chromosome -> absolute coordinate
    int GeneLen;
    map<int,orf_t> gs2orf;         // was: Kgene, Sgene
    map<orf_t,int> orf2gs;         // was: rKSgene

    inline orf_t gs2orf0(int i) const {
        map<int,orf_t>::const_iterator it = gs2orf.find(i);
	if(it==gs2orf.end()) cerr << "gs2orf0: unable to find " << i << endl;
	return it->second;
    }
    inline int orf2gs0(const orf_t& orf) const {
        map<orf_t,int>::const_iterator it = orf2gs.find(orf);
	if(it==orf2gs.end()) cerr << "orf2gs0: unable to find " << orf << endl;
	return it->second;
    }


//    genome_t coo2orf;
    map<orf_t,coord_t> orf2coo;   // was: rKS;  // ORF -> coordinates

    inline coord_t orf2coo0(const orf_t& orf) const {
        map<orf_t,coord_t>::const_iterator it = orf2coo.find(orf);
	if(it==orf2coo.end()) cerr << "orf2coo0: unable to find " << orf << endl;
	return it->second;
    }


    map<orf_t,int> sign;     // ORF orientation
    set<string> ChrCircular;  // set of curcular chromosomes

    int sgn(const orf_t& g) const {
	map<orf_t,int>::const_iterator is = sign.find(g);
	if( is == sign.end() ) return 0;
	else return is->second;
    }


//    map<orf_t,int> orf_mt;    // was: SO multiplicities of matched ORFs in S.

    /* Change coordinates to ORF order */
    void as2gs() {
	orf_t prev;
	int k=0;
	
	genome_t G = *this;
	clear();

	orf2coo.clear();
	orf2gs.clear();
	gs2orf.clear();
	
	GeneLen = 0;
	for(genome_t::const_iterator ik=G.begin();ik!=G.end();++ik) {
	    pair<coord_t,orf_t> p = *ik;
    
	    if(p.first.first!=prev) {
		k = 0;
		prev = p.first.first;
	    }
	    //erase(ik++);
	    p.first.second = k++;
	    insert(p);
	    orf2coo[p.second] = p.first;

	    gs2orf[GeneLen] = p.second;
	    orf2gs[p.second] = GeneLen++;
	}

/*	
	for(map<span_t,orf_t>::const_iterator ik=as2orf.begin();ik!=as2orf.end();++ik) {
	    const orf_t& orf = ik->second;
	    gs2orf[GeneLen] = orf;
	    orf2gs[orf] = GeneLen++;
	}
*/
    }
    
    void del_gene(const orf_t& g) {
	map<orf_t,coord_t>::iterator it = orf2coo.find(g);
	if(it != orf2coo.end()) {
	    erase(it->second);
	    orf2coo.erase(it);
	}
    }
    
    

};

//map<orf_t,span_t> Genome::orf2as;
//map<orf_t,int> Genome::orf2gs;
//map<orf_t,coord_t> Genome::orf2coo;

#endif
