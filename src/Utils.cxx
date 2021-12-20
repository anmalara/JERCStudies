#include "UHH2/JERCStudies/include/Utils.hpp"

bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}

int FindInVector(const std::vector<std::string>& vec, const std::string& el) {
  int index = -1;
  // Find given element in vector
  auto it = std::find(vec.begin(), vec.end(), el);
  if (it != vec.end()) index = distance(vec.begin(), it);
  return index;
}

void FilterVector(std::vector<double>& vec, double val, std::string invert) {
  static double val_ = 0;
  static bool invert_ = false;
  val_ = val;
  invert_ = (invert=="invert")? true : false;
  // std::cout << "Check " << (invert_? "x >= "+std::to_string(val_) : "x < "+std::to_string(val_)) << std::endl;
  vec.erase(std::remove_if( vec.begin(), vec.end(),[](const double& x) { return invert_? x > val_ : x < val_; }), vec.end());
}
