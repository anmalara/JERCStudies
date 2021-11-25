#pragma once

#include "UHH2/core/include/Hists.h"

#include "UHH2/JERCStudies/include/HistsBase.hpp"
#include "UHH2/JERCStudies/include/Utils.hpp"

class JERCStudiesHists: public HistsBase {
public:
  JERCStudiesHists(uhh2::Context & ctx, const std::string & dirname, const std::string & collection_);
  virtual void fill(const uhh2::Event&) override;
  virtual ~JERCStudiesHists();

private:
  std::string collection;
  uhh2::Event::Handle<std::vector<Jet> > h_jets;
  std::vector<std::string> flav_bins = {"uds", "c", "b", "g", "pu", "else"};
  std::vector<double> eta_bins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0};
  std::vector<double> pt_bins = {10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0};

};
