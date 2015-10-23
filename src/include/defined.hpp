//
// Created by pavel on 10/21/15.
//

#ifndef MGRA_DEFINED_HPP
#define MGRA_DEFINED_HPP

#include <list>
#include <vector>
#include <array>
#include <queue>
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
#include <tuple>

#include <functional>

#include <cassert>

#include <boost/algorithm/string.hpp>

using vertex_t = std::string;
vertex_t const Infty = "oo";

#include "utility/config_singl.hpp"
#include "utility/property.hpp"
#include "utility/sym_multihashmap.hpp"
#include "utility/equivalence.hpp"

#include "logger/logger.hpp"

#include "structures/tree.hpp"
#include "structures/mcolor.hpp"
#include "structures/mularcs.hpp"
#include "structures/genome.hpp"

enum build_type { default_algo, wgd_algo, assembly_algo, target_algo };

enum block_file_type_t { infercars, grimm };

namespace algo {

enum kind_stage {balance_k,
    simple_path_k,
    four_cycles_k,
    fair_edge_k,
    clone_k,
    fair_clone_edge_k,
    components_k,
    change_canform_k,
    bruteforce_k,
    blossomv_k,
    completion_k
};

}


#endif //MGRA_DEFINED_HPP
