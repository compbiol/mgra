//
// Created by Nikita Kartashov on 16/03/2015.
//

#ifndef HASH_HPP__
#define HASH_HPP__
namespace util {
  template <class T>
  inline void hash_combine(std::size_t& seed, const T& v) {
    std::hash <T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
}

#endif //HASH_HPP__
