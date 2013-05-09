#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>

typedef std::map<std::string, size_t> myhash; 

myhash counts; 

void read (std::vector<std::string>& name_all, std::vector<std::vector<std::string> >& genomes) { 
  std::ifstream myfile("ecoli.txt"); 
  
  if (myfile.is_open()) { 
    std::string temp; 
    std::getline(myfile, temp);

    while(myfile.good()) { 
      if (temp.find('>') != std::string::npos) { 
	name_all.push_back(temp); 
	std::getline(myfile, temp);

	std::string current = ""; 
	while(temp.find('>') == std::string::npos && myfile.good()) { 
	  current += temp; 
	  std::getline(myfile, temp); 
	} 

	std::vector<std::string> vcurr; 
	std::string cur = "";
	std::stringstream str(current);

	while(str.good()) { 
	  str >> cur; 
	  vcurr.push_back(cur);
	}
 
	genomes.push_back(vcurr);
      } 
    }
  }
  
  myfile.close();
}

void solve(const std::vector<std::vector<std::string> >& genomes) { 
  for (size_t i = 0; i < genomes.size(); ++i) { 
    for (auto it = genomes[i].cbegin(); it != genomes[i].cend(); ++it) { 
      if (counts.find(*it) == counts.end()) { 
	counts.insert(std::make_pair(*it, 1));
      } else { 
	counts[*it] += 1;
      }  
    } 
  } 
} 

void write(const std::vector<std::string>& name_all, const std::vector<std::vector<std::string> >& genomes) { 
  std::ofstream out("ecoli1.txt");
  
  for (size_t i = 0; i < genomes.size(); ++i) { 
    out << name_all[i] << std::endl;
    for (auto it = genomes[i].cbegin(); it != genomes[i].cend(); ++it) { 
      if (counts.find(*it)->second != 1) { 
	out << *it << " "; 
      }  
    } 
    
    out << std::endl; 
  } 

    
} 

int main() { 
  std::vector<std::string> name; 
  std::vector<std::vector<std::string> > genomes; 

  read(name, genomes);
  solve(genomes);
  write(name, genomes);
}
