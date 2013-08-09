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

struct Chromosome { 
  typedef std::string orf_t;
  typedef std::pair<orf_t, int> gene_t;

  Chromosome () 
  : isCircular(false) 
  { 
  }

  inline void insert(const orf_t& gene, size_t offset, int sign) { 
	auto gen = gene_t(gene, sign);
	main_chromosome.insert(std::make_pair(offset, gen));
  } 

  inline void do_circular() { 
	isCircular = true; 
  } 

  inline std::map<size_t, gene_t>::const_iterator begin() const { 
    return main_chromosome.cbegin();	
  } 

  inline std::map<size_t, gene_t>::const_iterator end() const { 
    return main_chromosome.cend();
  } 

  inline size_t size() const {
    return main_chromosome.size();      
  } 

  inline bool is_circular() const { 
    return isCircular;
  }
private: 
  bool isCircular; 
  std::map<size_t, gene_t> main_chromosome;
}; 

struct Genome {
  typedef std::string orf_t;
  typedef std::pair<std::string, size_t> coord_t; // (contig, offset)
  typedef std::pair<orf_t, int> gene_t;
  
  inline void insert(const orf_t& gene, const std::string& chromosome, size_t offset, int sign, size_t start, size_t end) {
    const coord_t& p = std::make_pair(chromosome, offset);
    const gene_t& orf = std::make_pair(gene, sign);   
    main_genome[chromosome].insert(gene, offset, sign);
  }   
	
  inline void registrate_circular_chr(std::string name) { 
    main_genome[name].do_circular();
  }  

  inline size_t size() const { 
    size_t sz = 0; 
    for (const auto &chromosome: main_genome) { 
      sz += chromosome.second.size(); 
    } 
    return sz; 
  }

  inline std::map<std::string, Chromosome>::const_iterator begin() const { 
    return main_genome.cbegin();	
  } 

  inline std::map<std::string, Chromosome>::const_iterator end() const { 
    return main_genome.cend();
  } 
private: 
  std::map<std::string, Chromosome> main_genome;
};

#endif
