#include "UHH2/JERCStudies/include/Utils.hpp"

bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}

int FindInVector(const std::vector<std::string>& vec, const std::string& el) {
  int index = -1;
  // Find given element in vector
  auto it = std::find(vec.begin(), vec.end(), el);
  if (it != vec.end()) index = distance(vec.begin(), it);
  return index;
}
