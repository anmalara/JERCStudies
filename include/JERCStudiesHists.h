#pragma once

#include "UHH2/core/include/Hists.h"

#include "UHH2/JERCStudies/include/HistsBase.hpp"
#include "UHH2/JERCStudies/include/Utils.hpp"

class JERCStudiesHists: public HistsBase {
public:
  JERCStudiesHists(uhh2::Context & ctx, const std::string & dirname, const std::string & collection_, const bool isPositive_=false, const bool isCentral_=false, const bool isHF_=false);
  virtual void fill(const uhh2::Event&) override;
  virtual ~JERCStudiesHists();

private:
  std::string collection;
  uhh2::Event::Handle<std::vector<Jet> > h_jets;
  bool isPositive, isCentral, isHF;
  std::vector<std::string> flav_bins = {"all", "light","heavy", "uds", "c", "b", "g", "pu", "else"};
  // std::vector<double> eta_bins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0};
  // std::vector<double> pt_bins = {10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0};

  std::vector<double> eta_bins = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, +0.000, +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, +4.363, +4.538, +4.716, +4.889, +5.191};
  std::vector<double> pt_bins = {15.0, 17.0, 20.0, 23.0, 27.0, 30.0, 35.0, 40.0, 45.0, 57.0, 72.0, 90.0, 120.0, 150.0, 200.0, 300.0, 400.0, 550.0, 750.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0};

};
