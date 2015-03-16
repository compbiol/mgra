//
// Created by Nikita Kartashov on 16/03/2015.
//

#ifndef _MGRA_MCOLOR_HASH_HPP_
#define _MGRA_MCOLOR_HASH_HPP_

#include "mcolor.hpp"

namespace std {
  template <> struct hash<structure::Mcolor>
  {
    size_t operator()(structure::Mcolor const& mcolor) const
    {
      return mcolor.make_hash();
    }
  };
}


#endif //_MGRA_MCOLOR_HASH_HPP_
