#ifndef DEFINED_H_
#define DEFINED_H_

#include <list>
#include <vector>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <set> 
#include <map>

#include <algorithm>

#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include <memory>
#include <limits>
#include <utility>
#include <cassert>

#include <functional>

#include "utility/sym_multi_hashmap.h"
#include "utility/equivalence.h"

typedef std::string vertex_t;
typedef utility::sym_multi_hashmap<vertex_t> edges_t;
typedef utility::sym_multi_hashmap<vertex_t> partgraph_t;
typedef std::list<vertex_t> path_t;
typedef std::pair<vertex_t, vertex_t> arc_t;

const vertex_t Infty = "oo"; 

template <typename T>
__attribute__((always_inline)) inline std::string toString(T val) { 
  std::ostringstream oss; 
  oss << val; 
  return oss.str();
}  

__attribute__((always_inline)) inline bool check_symbol(char c) { 
  if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9')) { 
    return true;
  } else { 
    return false;
  } 
} 

__attribute__((always_inline)) inline std::string trim(std::string s, const std::string& drop = " \t\r\n"){
  s = s.erase(s.find_last_not_of(drop) + 1);
  return s.erase(0, s.find_first_not_of(drop));
}

#include "genome.h"
#include "mcolor.h"
typedef std::pair<vertex_t, structure::Mcolor> edge_t; //FIXME CHECK on all
#include "mularcs.h"
//#include pconf and tree

#endif


