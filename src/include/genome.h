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


#ifndef GENOME_H_
#define GENOME_H_

#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <map>
#include <string>

typedef std::string orf_t;
typedef std::pair<int, int> span_t; 		// interval in absolute coordinates
typedef std::pair<std::string, int> coord_t;   	// (contig, offset)
typedef std::pair<std::string, span_t> cpan_t;  // (chr, start, end)

struct Genome {
	void insert(std::string gene, std::string chromosome, size_t offset, int sign, size_t start, size_t end);  
	
	/*inline void set_name(const std::string& name_) { 
		name = name_;
	} */

	inline void registrate_circular_chr(std::string name) { 
		circular_chromosome.insert(name);
	}  

	bool isCircular(std::string name) const;

	inline int get_sign(const orf_t& name) const { 
		auto is = sign.find(name);
		if (is == sign.end()) { 
			return 0;
		} else { 
			return is->second;
		} 
	}

	inline size_t size() const { 
		return main_genome.size();
	} 

	inline std::map<coord_t, orf_t>::const_iterator cbegin() const { 
		return main_genome.cbegin();	
	} 

	inline std::map<coord_t, orf_t>::const_iterator cend() const { 
		return main_genome.cend();
	} 

	/*inline std::string get_name() const { 
		return name;
	} */

	std::unordered_map<orf_t, cpan_t> orf2cpan;
private: 
	//std::string name;
	std::map<coord_t, orf_t> main_genome; //why not vector  
	std::unordered_map<orf_t, int> sign; 
	std::unordered_set<std::string> circular_chromosome;  //set of circular chromosomes
};
#endif
