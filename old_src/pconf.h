/* 
** Module: MGRA Configurations support
** Version: 1.1
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

#ifndef PCONF_H
#define PCONF_H


#include <vector>
#include <set>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#include "bpgraph.h"
#include "gen_cls.h"


string trim(string s, const string& drop = " \t\r\n") {
 s = s.erase(s.find_last_not_of(drop)+1);
 return s.erase(0,s.find_first_not_of(drop));
}
  

class ProblemInstance {
public:
    size_t NumGen;

    map<string,size_t> gen2num;

    vector<string> GenomeName;
    map<string,string> GenomeAlias;

    string blkformat, blkfile;

    list<string> trees;

    string graphfname;

    string colorscheme;

    size_t stages;

    string target;

    list< vector<string> > completion;
};


bool ReadProblem( const char* pfile, ProblemInstance& PI ) {

    map< string, list< string > > ProblemConfig;
    {
	ifstream PC(pfile);
	if( !PC ) {
	    cerr << "Cannot open " << pfile << endl;
	    return false;
	}

        string section;

	while(true) {
            string line;
	    getline( PC, line );
	    if( PC.eof() ) break;

            line = trim( line );
	    if( line.empty() ) continue;
	    if( line[0] == '#' ) continue;

	    if( line[0] == '[' && line[line.size()-1] == ']' ) {
		section = line;
	    }
	    else {
		ProblemConfig[section].push_back( line );
	    }
	}
	PC.close();
    }

//clog << "Sections: " << ProblemConfig.size() << endl;

    for( map< string, list< string > >::const_iterator ip=ProblemConfig.begin(); ip!=ProblemConfig.end(); ++ip ) {

//clog << ip->first << endl;

	if( ip->first == "[Genomes]" ) {
	    
	    PI.NumGen = ip->second.size();
	    PI.GenomeName.resize(PI.NumGen);
	    
	    size_t k = 0;
	    
	    for( list<string>::const_iterator js=ip->second.begin(); js!=ip->second.end(); ++js, ++k ) {
		istringstream is( *js );
		string name;
		is >> name;
		
		PI.gen2num[ name ] = k;
		PI.GenomeName[ k ] = name;

		if( PI.GenomeAlias.find( name ) != PI.GenomeAlias.end() ) {
		    cerr << "ERROR: Genome identificator " << name << " is not unique!" << endl;
		    return false;
		}
		PI.GenomeAlias[ name ] = name;
		while( !is.eof() ) {
		    string alias;
		    is >> alias;

		    if( PI.GenomeAlias.find( alias ) != PI.GenomeAlias.end() ) {
			cerr << "ERROR: Duplicate alias " << alias << endl;
			return false;
		    }
		    PI.GenomeAlias[ alias ] = name;
		    PI.gen2num[ alias ] = k;
		}
	    }

	    continue;
	}

	if( ip->first == "[Blocks]" ) {

	    for( list<string>::const_iterator js=ip->second.begin(); js!=ip->second.end(); ++js ) {
		istringstream is( *js );
		string name;
		is >> name;

		if( name == "format" ) is >> PI.blkformat;
		else if( name == "file" ) is >> PI.blkfile;
		else {
		    cerr << "Unknown option " << name << endl;
		    return false;
		}
	    }

            continue;
	}

	if( ip->first == "[Trees]" ) {

	    PI.trees = ip->second;

            continue;
	}

	if( ip->first == "[Graphs]" ) {

	    for( list<string>::const_iterator js=ip->second.begin(); js!=ip->second.end(); ++js ) {
		istringstream is( *js );
		string name;
		is >> name;

		if( name == "filename" ) is >> PI.graphfname;
		else if( name == "colorscheme" ) is >> PI.colorscheme;
		else {
		    cerr << "Unknown option " << name << endl;
		    return false;
		}
	    }


            continue;
	}

	if( ip->first == "[Algorithm]" ) {

	    for( list<string>::const_iterator js=ip->second.begin(); js!=ip->second.end(); ++js ) {
		istringstream is( *js );
		string name;
		is >> name;

		if( name == "stages" ) is >> PI.stages;
		else if( name == "target" ) is >> PI.target;
		else {
		    cerr << "Unknown option " << name << endl;
		    return false;
		}
	    }

            continue;
	}

	if( ip->first == "[Completion]" ) {

	    for( list<string>::const_iterator js=ip->second.begin(); js!=ip->second.end(); ++js ) {

                vector<string> mc(5);

		istringstream is( *js );

		is >> mc[0] >> mc[1] >> mc[2] >> mc[3] >> mc[4];

                PI.completion.push_back(mc);
	    }

            continue;
	}

	cerr << "Unknown section " << ip->first << endl;
        return false;

    }
    return true;
}


void read_infercars(const ProblemInstance& PI, vector<Genome>& G) {

    const size_t NN = PI.NumGen;
    G.resize(NN);

    map<size_t,int> SegLen;
    
    for(int i=0;i<NN;++i) {
	G[i].name = PI.GenomeName[i];
    }

    // chromosome coverage
    map<string,size_t> chrlen;
    size_t chrcov = 0;

    {
	ifstream fh(PI.blkfile.c_str());
	if(!fh) {
	    cerr << "Unable to open " << PI.blkfile << endl;
	    exit(1);
	}
    
	int ngen = 0;
	string gene;
	vector<string> chr(NN), sign(NN);
	vector<int> st(NN), en(NN), c(NN);
	while(true) {
	    string line;
	    getline(fh,line);
	    if(fh.eof()) break;
	    
	    line = trim(line);
	    
	    if(line.empty()) {
	    
		if( gene.empty() ) continue;
	    
		int n = 0;
		size_t len = 0;
		
		bool ambig = false;
		for(int i=0;i<NN;++i) if(c[i]!=1) ambig = true;
		if( ambig ) {
		    clog << "Ambiguous block: " << gene << endl;
		    continue;
		}

/*
#ifdef MRDQHC
		size_t Xcount = 0;
		for(int i=0;i<NN;++i) if(chr[i]=="chrX") Xcount++;
		if( Xcount!=0 && Xcount!=NN ) {
		    clog << "Inconsistent X chromosome block: " << gene << endl;
		    continue;
		}
#endif
*/

		for(int i=0;i<NN;++i) if(c[i] == 1) {
		    ++n;
		    len += abs(en[i] - st[i]) + 1;
#ifdef BLKLEN_THRESHOLD
if(i==0 && len < BLKLEN_THRESHOLD) break;
#endif
		    coord_t p = make_pair(chr[i],(st[i]+en[i])/2);                                                       
        	    G[i][p] = gene;                                                                                                 
        	    if( sign[i] == "+" ) G[i].sign[gene] = +1;
		    else G[i].sign[gene] = -1;
		    // saving block span info
		    G[i].orf2cpan[gene] = make_pair(chr[i],make_pair(min(st[i],en[i]),max(st[i],en[i])));
		}
		SegLen[len/n]++;

                /*
		if( chr[0]=="chrX" ) {
                    chrX.insert(gene+"t");
		    chrX.insert(gene+"h");
		}
                */

		chrcov += abs(en[0] - st[0]) + 1;
		chrlen[chr[0]] = max(chrlen[chr[0]],(size_t)max(abs(st[0]),abs(en[0])));
		
		gene.clear();
		
		continue;
	    }
	    
	    if(line[0] == '#') continue;
	    
	    if(line[0] == '>') {
		gene = line.substr(1);
		for(int i=0;i<NN;++i) c[i] = 0;
		continue;
	    }
	    
	    line[line.find(".")] = ' ';
	    line[line.find(":")] = ' ';
	    line[line.find("-")] = ' ';
	    
	    istringstream istr(line);
	    
	    string gename;
	    istr >> gename;
	    if( PI.gen2num.find(gename)==PI.gen2num.end() ) {
		cerr << "Unknown genome name: " << gename << endl;
		continue;
	    }
	    size_t k = PI.gen2num.find(gename)->second;

	    /*
	    {
		string chr1;
		int st1, en1, sign1;
		istr >> chr1 >> st1 >> en1 >> sign1;
		if( c[k]>0 ) {
		    if(chr1 != chr[k]) {
		      clog << "Multiple gene: " << gene << endl;
		      c[k]++;
		      
		st[k] = min(st[k],st1);
		en[k] = max(en[k],en1);
	    }
	    */

	    istr >> chr[k] >> st[k] >> en[k] >> sign[k];
	    ++c[k];
	}
	fh.close();
    }

    /*
    size_t totlen = 0;
    for(map<string,size_t>::const_iterator ic=chrlen.begin();ic!=chrlen.end();++ic) {
	totlen += ic->second;
    }
    cout << "Coverage: " << (double)chrcov/(double)totlen << endl;
    */
}






void read_grimm(const ProblemInstance& PI, vector<Genome>& G) {

    G.resize( PI.NumGen );
    for(int i=0;i<PI.NumGen;++i) {
	G[i].name = PI.GenomeName[i];
    }
                    
    int n;

    ifstream fh(PI.blkfile.c_str());
    
    //ofstream joins;

    if(!fh) {
	cerr << "Unable to open " << PI.blkfile << endl;
	exit(1);
    }

    //cout << "Reading " << PI.blkfile << endl;

    int nchr = 0;
    int ngen = 0;

    bool genomefound = false;

    while(true) {


	string line;
	getline(fh,line);
	if(fh.eof()) break;

	line = trim(line);

	static string chr;
	
	if(line.empty()) continue;

	if( line[0]=='#' ) {
	    //chr = trim(line.substr(1));
	    continue;
	}

	if( line[0]=='>' ) {
	    line = trim(line.substr(1));

	    genomefound = false;

	    if( PI.gen2num.find(line)!=PI.gen2num.end() ) {
		n = PI.gen2num.find(line)->second;
		genomefound = true;
		continue;
	    }

	    clog << "Unknown genome: " << line << endl;
	}

	if( !genomefound ) continue;

	if( chr.empty() ) {
	    ostringstream ostr;
	    ostr << "chr" << ++nchr;
	    chr = ostr.str();
	    //cerr << "Empty chromosome name" << endl;
	    //exit(1);
	}

	istringstream istr(line);
	
	//string prev = "0";
	int genord = 0;
	while( true ) {
	    string gen;
	    istr >> gen;
	    if( istr.eof() && gen.empty() ) break;
	    if( gen[0]=='$' ) break;
	    if( gen[0]=='@' ) {
		G[n].ChrCircular.insert(chr);
		break;
	    }

	    //joins << "\t" << prev << "\t" << gen << endl;
	    //prev = gen;
	    
	    bool minus = gen[0]=='-';
	    if( minus || gen[0]=='+' ) {
		gen = gen.substr(1);
	    }

	    coord_t p = make_pair(chr, ++genord );
	    G[n][p] = gen;
	    if( minus ) G[n].sign[gen] = -1; else G[n].sign[gen] = +1;
	    G[n].orf2cpan[gen] = make_pair(chr,make_pair(genord,genord));

	    /*
	    if( chr == "Chromosome X" ) {
		chrX.insert(gen+"t");
		chrX.insert(gen+"h");
	    }
	    */

	    ngen++;
	}
	
	//joins << "\t" << prev << "\t0" << endl;
	//prev = "0";
	
	chr.clear();
    }
    
    //if( genomefound ) joins.close();
    
    fh.close();

    //cout << "blocks read: " << ngen << endl;

    return;
}


#endif
