#include "TH1F.h"
#include "TH2F.h"

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/JERCStudies/include/JERCStudiesHists.h"

using namespace std;
using namespace uhh2;


JERCStudiesHists::JERCStudiesHists(Context & ctx, const string & dname, const string & collection_): HistsBase(ctx, dname), collection(collection_){

  h_jets = ctx.get_handle<vector<Jet>>(collection);

  for(unsigned int f = 0; f <flav_bins.size(); f++){
    for(unsigned int i = 0; i <eta_bins.size(); i++){
      for(unsigned int j = 0; j <pt_bins.size(); j++){
        std::string name = "Resp_";
        name += "flav_"+flav_bins[f];
        name += "eta_"+GetStringFromFloat(eta_bins[i])+"to"+GetStringFromFloat(eta_bins[i+1]);
        name += "pt_"+GetStringFromFloat(pt_bins[j])+"to"+GetStringFromFloat(pt_bins[j+1]);
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
    double eta_jet = jet.eta();
    double eta_genjet = next_genjet->eta();
    double pt_jet = jet.pt();
    double pt_genjet = next_genjet->pt();
    double resp = pt_jet/pt_genjet;

    string flav = "else";
    if (gen_flav==0) flav = "pu";
    else if (gen_flav<=3) flav = "uds";
    else if (gen_flav==4) flav = "c";
    else if (gen_flav==5) flav = "b";
    else if (gen_flav==21) flav = "g";

    for(unsigned int i = 0; i <eta_bins.size(); i++){
      if (eta_bins[i]<eta_genjet && eta_genjet<eta_bins[i+1]){
        for(unsigned int j = 0; j <pt_bins.size(); j++){
          if (pt_bins[j]<pt_genjet && pt_genjet<pt_bins[j+1]){
            std::string name = "Resp_";
            name += "flav_"+flav;
            name += "eta_"+GetStringFromFloat(eta_bins[i])+"to"+GetStringFromFloat(eta_bins[i+1]);
            name += "pt_"+GetStringFromFloat(pt_bins[j])+"to"+GetStringFromFloat(pt_bins[j+1]);
            fill_H1(name, resp, weight);
          }
        }
      }
    }
  }
}

JERCStudiesHists::~JERCStudiesHists(){}
