#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"

bool FindInString(const std::string& search, const std::string& str);
int FindInVector(const std::vector<std::string>& vec, const std::string& el);

inline const char* BoolToString(bool b) { return b ? "true" : "false";}
inline std::string GetStringFromFloat(float x) {return std::to_string(int(x))+"p"+std::to_string(int((x-int(x))*10));}
