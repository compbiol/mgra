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

#include "defined.h"
 
struct Genome {
  typedef std::string orf_t;
  typedef std::pair<size_t, size_t> span_t; 	// interval in absolute coordinates
  typedef std::pair<std::string, size_t> coord_t; // (contig, offset)
  typedef std::pair<orf_t, int> gene_t;
  typedef std::pair<std::string, span_t> cpan_t;  // (chr, start, end)

  inline void insert(const orf_t& gene, const std::string& chromosome, size_t offset, int sign, size_t start, size_t end) {
    const coord_t& p = std::make_pair(chromosome, offset);
    const gene_t& orf = std::make_pair(gene, sign);   
    main_genome.insert(std::make_pair(p, orf));
  }   
	
  inline void registrate_circular_chr(std::string name) { 
    circular_chromosome.insert(name);
  }  

  inline bool isCircular(std::string name) const {
    return (circular_chromosome.count(name) != 0);
  } 

  inline size_t size() const { 
    return main_genome.size();
  } 

  inline std::map<coord_t, gene_t>::const_iterator cbegin() const { 
    return main_genome.cbegin();	
  } 

  inline std::map<coord_t, gene_t>::const_iterator cend() const { 
    return main_genome.cend();
  } 
private: 
  std::map<coord_t, gene_t> main_genome; 
  std::unordered_set<std::string> circular_chromosome;  //set of circular chromosomes
};

#endif
