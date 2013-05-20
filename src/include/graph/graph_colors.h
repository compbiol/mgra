#ifndef GRAPH_COLOR_H_
#define GRAPH_COLOR_H_

#include "genome_match.h"
#include "utility/sym_map.h"

template<class mcolor_t>
struct ColorsGraph { 
  ColorsGraph(size_t size, const ProblemInstance& cfg); 

  void update_complement_color(const std::vector<mcolor_t>& colors);

  inline mcolor_t get_complement_color(const mcolor_t& color) const { 
    assert(CColorM.find(color) != CColorM.end());
    return CColorM.find(color)->second;
  } 
  
  inline mcolor_t get_min_complement_color(const mcolor_t& color) const {
    mcolor_t temp = get_complement_color(color);
    if (temp.size() > color.size() || (temp.size() == color.size() && temp > color)) {	
      return color;
    } 
    return temp;
  }	
  
  inline bool is_T_consistent_color(const mcolor_t& col) const { 
    return (all_T_color.find(col) != all_T_color.end());
  } 

  inline bool is_good_color(const mcolor_t& Q, const mcolor_t& Q1) const { 
	if (!is_T_consistent_color(Q)) { 
		return false; 
	} 

	if (!is_T_consistent_color(Q1)) { 
		return false; 
	} 

	if (get_complement_color(Q) == Q1) { 
		return true; 
	} 

	return false;
  } 


  bool are_adjacent_branches(const mcolor_t& A, const mcolor_t & B) const;

  mcolor_t CColorRep(const mcolor_t& c);
  const mcolor_t& CColor(const mcolor_t& S);

private: 
  void parsing_tree(size_t size, const ProblemInstance& cfg);
  mcolor_t add_tree(const std::string& tree, std::vector<std::string>& output);

private:
  size_t count_local_color;
  sym_map<mcolor_t> CColorM; //complementary multicolor
  std::set<mcolor_t> all_T_color;	//all T-consistent colors
public:
  std::set<mcolor_t> DiColor;   // directed colors
  std::vector<mcolor_t> TColor; // colors corresponding to ancestral genomes. nodes in tree
}; 

template<class mcolor_t>
ColorsGraph<mcolor_t>::ColorsGraph(size_t size, const ProblemInstance& cfg) 
:count_local_color(size) {
  parsing_tree(count_local_color, cfg);

  if (!cfg.get_target().empty()) { 
    DiColor.erase(genome_match::name_to_mcolor(cfg.get_target()));
  } 
  
  //check consistency
  for (auto id = DiColor.cbegin(); id != DiColor.cend(); ++id) {
    for(auto jd = id; jd != DiColor.end(); ++jd) {
      mcolor_t C(*id, *jd, mcolor_t::Intersection);
      if (!C.empty() && C.size() != id->size() && C.size() != jd->size()) {
	std::clog << "Multicolors " << genome_match::mcolor_to_name(*id) << " " << genome_match::mcolor_to_name(*jd) << " have nontrivial intersection, removing the latter" << std::endl;
	DiColor.erase(jd++);
	--jd;
      }
    }
  }
  
  TColor.resize(DiColor.size());
  
  std::clog << "vecT-consistent colors: " << DiColor.size() << std::endl;

  size_t col = 1;
  for (auto id = DiColor.begin(); id != DiColor.end(); ++id) {
    std::clog << "\t" << genome_match::mcolor_to_name(*id);
    all_T_color.insert(*id);
    
    // compute complement to *id
    mcolor_t C;
    for(size_t j = 0; j < size; ++j) {
      if (!(id->mymember(j))) { //IS HERE MEMBER
	C.insert(j);
      } 
    }
    all_T_color.insert(C);
    
    TColor[col-1] = *id;
    
    col++;
  }
  std::clog << std::endl;
  
  // tell where the root resides
  std::clog << "the root resides in between:";
  auto T = DiColor;
  for (auto it = T.begin(); it != T.end(); ++it) {
    for (auto jt = it; ++jt != T.end(); ) {
      mcolor_t C(*it, *jt, mcolor_t::Intersection);
      if (C.size() == it->size()) {
	T.erase(it++);
	jt = it;
	continue;
      }
      if (C.size() == jt->size()) {
	T.erase(jt++);
	--jt;
      }
    }
    std::clog << " " << genome_match::mcolor_to_name(*it);
  }
  std::clog << std::endl;
}

template<class mcolor_t>
mcolor_t ColorsGraph<mcolor_t>::add_tree(const std::string& tree, std::vector<std::string>& output) {
  if (tree[0] == '(') {
    //non-trivial tree
    if (tree[tree.size() - 1] != ')') {
      std::cerr << "ERROR: Malformed input (sub)tree 1" << std::endl;
      exit(3);
    }

    int p = 0;
    for(size_t j = 1; j < tree.size() - 1; ++j) {
      if (tree[j] == '(') { 
	++p; 
      } else if (tree[j] == ')') {
	--p;
      } else if (tree[j] == ',') { 
	if (p == 0) { 
	  mcolor_t Q1 = add_tree(tree.substr(1, j - 1), output);
	  mcolor_t Q2 = add_tree(tree.substr(j + 1, tree.size() - j - 2), output);

	  DiColor.insert(Q1);
	  DiColor.insert(Q2);

	  mcolor_t Q(Q1, Q2, mcolor_t::Union);
		    
	  output.push_back("\t\"" + genome_match::mcolor_to_name(Q) + "\"\t->\t\"" + genome_match::mcolor_to_name(Q1) + "\";");
	  output.push_back("\t\"" + genome_match::mcolor_to_name(Q) + "\"\t->\t\"" + genome_match::mcolor_to_name(Q2) + "\";");

	  return Q;
	} 
      } 
      if (p < 0) {
	std::cerr << "ERROR: Malformed input (sub)tree 2" << std::endl;
	exit(3);
      }
    }
    if (p != 0) {
      std::cerr << "ERROR: Malformed input (sub)tree 3" << std::endl;
      exit(3);
    }
  } else {
    //single node
    mcolor_t Q;
    for(size_t j = 0; j < tree.size(); ++j) {
      std::string c = tree.substr(j, 1);
      if (!genome_match::member_name(c)) {
	std::cerr << "ERROR: Unknown genome in (sub)tree: " << tree << std::endl;
	exit(3);
      }
      Q.insert(genome_match::get_number(c));
    }
    return Q;
  }
}

template<class mcolor_t>
const mcolor_t& ColorsGraph<mcolor_t>::CColor(const mcolor_t& S) {
  if (!CColorM.defined(S)) {
    mcolor_t T;
    for (size_t j = 0; j < count_local_color; ++j) {
      if (!S.mymember(j)) {  
	T.insert(j);
      } 
    }
    CColorM.insert(S, T);
  }
  return CColorM[S];
}

template<class mcolor_t>
mcolor_t ColorsGraph<mcolor_t>::CColorRep(const mcolor_t& c) {
  mcolor_t Q = CColor(c);
  if (Q.size() > c.size() || (Q.size() == c.size() && Q > c)) 	
    return c;
  return Q;
}

template<class mcolor_t>
void ColorsGraph<mcolor_t>::parsing_tree(size_t size, const ProblemInstance& cfg) { 
  std::vector<std::string> trees = cfg.get_trees();

  // add terminal branches	
  for (size_t j = 0; j < size; ++j) {
    DiColor.insert(mcolor_t(j));
  }

  std::vector<std::string> output;
  for(auto it = trees.cbegin(); it != trees.cend(); ++it) {
    mcolor_t C = add_tree(*it, output);
    if (C.size() < size) { 
      DiColor.insert(C); // complete multicolor is excluded
    } 
  }

  std::ofstream flegend("legend.dot");
  flegend << "digraph Legend {" << std::endl;

  flegend << "\tnode [style=filled];" << std::endl;

  for (size_t j = 0; j < size; ++j) {
    flegend << "\t\"" << cfg.get_name(j) << "\"\t[fillcolor=" <<  cfg.get_RGBcolor(cfg.get_RGBcoeff() * j)  << "];" << std::endl;
  } 


  for(auto it = output.cbegin(); it != output.cend(); ++it) {
    flegend << *it << std::endl;
  } 
  flegend << "}" << std::endl;
  flegend.close();
} 

template<class mcolor_t>
void ColorsGraph<mcolor_t>::update_complement_color(const std::vector<mcolor_t>& colors) {
  for(auto it = colors.begin(); it != colors.end(); ++it) { 
    if (!CColorM.defined(*it)) { 
      Mcolor temp; 
      for(size_t j = 0; j < count_local_color; ++j) { 
	if ((!it->mymember(j))) { 
	  temp.insert(j);
	} 
      } 
      CColorM.insert(*it, temp);
    } 
  } 
} 

template<class mcolor_t>
bool ColorsGraph<mcolor_t>::are_adjacent_branches(const mcolor_t& A, const mcolor_t & B) const {
	if (!is_T_consistent_color(A) || !is_T_consistent_color(B)) { 
		return false;
	} 

	mcolor_t Q1; 
	mcolor_t Q2;

	if (A.size() >= B.size()) {
		Q1 = A;
		Q2 = B;
	} else {
		Q1 = B;
		Q2 = A;
	}

	mcolor_t C(Q1, Q2, Mcolor::Difference);
	
	if (C.size() == Q1.size() - Q2.size() && is_T_consistent_color(C)) { 		
		return true;
	} 

	mcolor_t M(Q1, Q2, Mcolor::Union); 
	if (M.size() == Q1.size() + Q2.size() && is_T_consistent_color(M)) { 	
		return true;
	} 

	return false;
}


#endif

