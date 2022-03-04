#pragma once

#include "UHH2/core/include/Hists.h"

#include "UHH2/JERCStudies/include/HistsBase.hpp"
#include "UHH2/JERCStudies/include/Utils.hpp"

class JERCStudiesHists: public HistsBase {
public:
  JERCStudiesHists(uhh2::Context & ctx, const std::string & dname, const std::string & collection_, const float eta_min_, const float eta_max_, const float pu_min_, const float pu_max_, const std::string flav_);
  virtual void fill(const uhh2::Event&) override;
  virtual ~JERCStudiesHists();

private:
  std::string collection;
  uhh2::Event::Handle<std::vector<Jet> > h_jets;
  std::vector<double> pt_bins;
  float eta_min, eta_max, pu_min, pu_max;
  std::string flav;
  float dr_min = 0.2;

};
