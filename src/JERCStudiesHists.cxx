#include "TH1F.h"
#include "TH2F.h"

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/JERCStudies/include/JERCStudiesHists.h"

using namespace std;
using namespace uhh2;


JERCStudiesHists::JERCStudiesHists(Context & ctx, const string & dname, const string & collection_, const bool isPositive_, const bool isCentral_, const bool isHF_):
HistsBase(ctx, dname), collection(collection_), isPositive(isPositive_), isCentral(isCentral_), isHF(isHF_) {

  FilterVector(eta_bins, 0,                             isPositive? "":"invert");
  FilterVector(eta_bins, (isPositive?  1:-1)*etaBarrel, isPositive==isCentral? "invert":"");
  FilterVector(eta_bins, (isPositive? 1:-1)*etaHF,      isPositive==isHF? "":"invert");

  h_jets = ctx.get_handle<vector<Jet>>(collection);

  for(unsigned int f = 0; f <flav_bins.size(); f++){
    for(unsigned int e = 0; e <eta_bins.size(); e++){
      for(unsigned int p = 0; p <pt_bins.size(); p++){
        std::string name = "Resp_";
        name += "flav_"+flav_bins[f];
        name += "eta_"+GetStringFromFloat(eta_bins[e])+"to"+GetStringFromFloat(eta_bins[e+1]);
        name += "pt_"+GetStringFromFloat(pt_bins[p])+"to"+GetStringFromFloat(pt_bins[p+1]);
        book_TH1F(name, "R", 100, -1, 3);
      }
    }
  }
}


void JERCStudiesHists::fill(const Event & event){
  double weight = event.weight;

  for (const auto & jet : event.get(h_jets)) {

    auto next_genjet = closestParticle(jet, *event.genjets);
    auto drmin = next_genjet ? deltaR(jet, *next_genjet) : numeric_limits<float>::infinity();
    if (drmin>0.2) continue;

    double gen_flav = abs(next_genjet->partonFlavour());
    double eta_genjet = next_genjet->eta();

    if (isPositive && eta_genjet<0) continue;
    if (isCentral && eta_genjet<0) continue;

    double pt_jet = jet.pt();
    double pt_genjet = next_genjet->pt();
    double resp = pt_jet/pt_genjet;

    string flav = "else";
    if (gen_flav==0) flav = "pu";
    else if (gen_flav<=3) flav = "uds";
    else if (gen_flav==4) flav = "c";
    else if (gen_flav==5) flav = "b";
    else if (gen_flav==21) flav = "g";

    string flav2 = "";
    if (gen_flav==4 || gen_flav==5) flav2 = "heavy";
    else if (gen_flav>=1 || gen_flav<=3 || gen_flav==21) flav2 = "light";

    string flav3 = "";
    if (gen_flav!=0) flav3 = "all";

    for(unsigned int e = 0; e <eta_bins.size(); e++){
      if (eta_bins[e]<eta_genjet && eta_genjet<eta_bins[e+1]){
        for(unsigned int p = 0; p <pt_bins.size(); p++){
          if (pt_bins[p]<pt_genjet && pt_genjet<pt_bins[p+1]){
            std::string name = "Resp_";
            name += "flav_"+flav;
            name += "eta_"+GetStringFromFloat(eta_bins[e])+"to"+GetStringFromFloat(eta_bins[e+1]);
            name += "pt_"+GetStringFromFloat(pt_bins[p])+"to"+GetStringFromFloat(pt_bins[p+1]);
            fill_H1(name, resp, weight);
            if (flav2!="") {
              name = "Resp_";
              name += "flav_"+flav2;
              name += "eta_"+GetStringFromFloat(eta_bins[e])+"to"+GetStringFromFloat(eta_bins[e+1]);
              name += "pt_"+GetStringFromFloat(pt_bins[p])+"to"+GetStringFromFloat(pt_bins[p+1]);
              fill_H1(name, resp, weight);
            }
            if (flav3!="") {
              name = "Resp_";
              name += "flav_"+flav3;
              name += "eta_"+GetStringFromFloat(eta_bins[e])+"to"+GetStringFromFloat(eta_bins[e+1]);
              name += "pt_"+GetStringFromFloat(pt_bins[p])+"to"+GetStringFromFloat(pt_bins[p+1]);
              fill_H1(name, resp, weight);
            }
          }
        }
      }
    }
  }
}

JERCStudiesHists::~JERCStudiesHists(){}
