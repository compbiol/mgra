#ifndef DEFINED_H_
#define DEFINED_H_

#include <list>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set> 
#include <map>

#include <algorithm>

#include <iostream>
#include <fstream>
#include <sstream>

#include <cassert>

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


#endif


