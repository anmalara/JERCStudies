#include "TH1F.h"
#include "TH2F.h"

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/JERCStudies/include/JERCStudiesHists.h"
#include "UHH2/JERCProtoLab/macros/common_info/common_binning.hpp"

using namespace std;
using namespace uhh2;


JERCStudiesHists::JERCStudiesHists(Context & ctx, const string & dname, const string & collection_,  const float eta_min_, const float eta_max_, const float pu_min_, const float pu_max_, const string flav_):
HistsBase(ctx, dname), collection(collection_), eta_min(eta_min_), eta_max(eta_max_), pu_min(pu_min_), pu_max(pu_max_), flav(flav_) {

  pt_bins.reserve(JERC_Constants::pt_bins_MCTruth.size());
  pt_bins.insert(pt_bins.end(), JERC_Constants::pt_bins_MCTruth.begin(), JERC_Constants::pt_bins_MCTruth.end());

  h_jets = ctx.get_handle<vector<Jet>>(collection);

  for(unsigned int p = 0; p <pt_bins.size(); p++){
    std::string name = "pt_"+GetStringFromFloat(pt_bins[p])+"to"+GetStringFromFloat(pt_bins[p+1]);
    book_TH1F("Resp_"+name, "R", 100, -1, 3);
    book_TH1F("pt_"+name, "pt", 250,10,510);
  }

}


void JERCStudiesHists::fill(const Event & event){
  double weight = event.weight;

  for (const auto & jet : event.get(h_jets)) {

    auto next_genjet = closestParticle(jet, *event.genjets);
    auto dr_min_ = next_genjet ? deltaR(jet, *next_genjet) : numeric_limits<float>::infinity();
    if (dr_min_>dr_min) continue;

    double eta_genjet = next_genjet->eta();
    double gen_flav = abs(next_genjet->partonFlavour());
    double pu_gen = event.genInfo->pileup_TrueNumInteractions();

    if (eta_min>eta_genjet || eta_genjet>eta_max) continue;
    if (pu_min>pu_gen || pu_gen>pu_max) continue;
    if (gen_flav==0 && flav != "pu" && flav !="all") continue;
    if ((gen_flav==4 || gen_flav==5) && flav != "heavy" && flav !="all") continue;
    if ((gen_flav>=1 && gen_flav<=3)  && flav != "light" && flav !="all") continue;
    if (gen_flav==21 && flav != "g" && flav !="all") continue;

    double pt_jet = jet.pt();
    double pt_genjet = next_genjet->pt();
    double resp = pt_jet/pt_genjet;

    for(unsigned int p = 0; p <pt_bins.size(); p++){
      if (pt_bins[p]<pt_genjet && pt_genjet<pt_bins[p+1]){
        std::string name = "pt_"+GetStringFromFloat(pt_bins[p])+"to"+GetStringFromFloat(pt_bins[p+1]);
        fill_H1("Resp_"+name, resp,    weight);
        fill_H1("pt_"+name,   pt_jet,  weight);
      }
    }

  }
}

JERCStudiesHists::~JERCStudiesHists(){}
