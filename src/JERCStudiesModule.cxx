#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/JERCStudies/include/ModuleBase.h"
#include "UHH2/JERCStudies/include/Utils.hpp"
#include "UHH2/JERCStudies/include/JERCStudiesSelections.h"
#include "UHH2/JERCStudies/include/JERCStudiesHists.h"
#include "UHH2/JERCStudies/include/GenericJetCleaner.h"

using namespace std;

class JERCStudiesModule: public ModuleBASE {
public:

  explicit JERCStudiesModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  // void book_handles(uhh2::Context&);
  void PrintInputs();

protected:

  // Define variables
  std::string NameModule = "JERCStudiesModule";
  std::vector<std::string> histogram_tags = { "nocuts", "weights", "cleaned"};
  std::vector<string> myTags;
  std::unordered_map<std::string, std::string> myStrings;
  std::unordered_map<std::string, bool> myBools;

  Event::Handle<std::vector<Jet> > h_jets_CHS;
  Event::Handle<std::vector<Jet> > h_jets_Puppi;
  // Define common modules

  std::vector<std::unique_ptr<AnalysisModule>> weightsmodules, modules;
  std::unique_ptr<uhh2::AndSelection> metfilters_selection;

  std::unique_ptr<GenericJetCleaner> jet_cleaner_CHS, jet_cleaner_Puppi;

};


void JERCStudiesModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "             "+NameModule+"             " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : myStrings) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : myBools) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << BoolToString(x.second) << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void JERCStudiesModule::book_histograms(uhh2::Context& ctx) {
  for(const auto & tag : histogram_tags){
    string tag_, mytag;
    for(const string & algo : {"CHS","Puppi"}){
      tag_ = "Response_"+algo+"_Neg_HF_"; myTags.push_back(tag_); mytag = tag_ + tag; book_HFolder(mytag, new JERCStudiesHists(ctx,mytag, myStrings["jetLabel_"+algo], false, false, true));
      tag_ = "Response_"+algo+"_Neg_EC_"; myTags.push_back(tag_); mytag = tag_ + tag; book_HFolder(mytag, new JERCStudiesHists(ctx,mytag, myStrings["jetLabel_"+algo], false, false, false));
      tag_ = "Response_"+algo+"_Neg_BA_"; myTags.push_back(tag_); mytag = tag_ + tag; book_HFolder(mytag, new JERCStudiesHists(ctx,mytag, myStrings["jetLabel_"+algo], false, true,  false));
      tag_ = "Response_"+algo+"_Pos_BA_"; myTags.push_back(tag_); mytag = tag_ + tag; book_HFolder(mytag, new JERCStudiesHists(ctx,mytag, myStrings["jetLabel_"+algo], true,  true,  false));
      tag_ = "Response_"+algo+"_Pos_EC_"; myTags.push_back(tag_); mytag = tag_ + tag; book_HFolder(mytag, new JERCStudiesHists(ctx,mytag, myStrings["jetLabel_"+algo], true,  false, false));
      tag_ = "Response_"+algo+"_Pos_HF_"; myTags.push_back(tag_); mytag = tag_ + tag; book_HFolder(mytag, new JERCStudiesHists(ctx,mytag, myStrings["jetLabel_"+algo], true,  false, true));
    }
  }
}

void JERCStudiesModule::fill_histograms(uhh2::Event& event, string tag){
  for (auto& mytag : myTags) HFolder(mytag+ tag)->fill(event);
}

JERCStudiesModule::JERCStudiesModule(Context & ctx){

  // Set up variables
  myStrings["year"]            = ctx.get("year");
  myStrings["dataset_version"] = ctx.get("dataset_version");
  myBools["is_mc"]             = ctx.get("dataset_type") == "MC";
  myStrings["SysType_PU"]      = ctx.get("SysType_PU");
  myBools["lumisel"]           = string2bool(ctx.get("lumisel"));
  myBools["mclumiweight"]      = string2bool(ctx.get("mclumiweight"));
  myBools["mcpileupreweight"]  = string2bool(ctx.get("mcpileupreweight"));
  myBools["eleid"]             = string2bool(ctx.get("eleid"));
  myBools["muid"]              = string2bool(ctx.get("muid"));
  myBools["tauid"]             = string2bool(ctx.get("tauid"));
  myBools["metfilters"]        = string2bool(ctx.get("metfilters"));

  myStrings["jetLabel_CHS"] = "jets";
  myStrings["jetLabel_Puppi"] = "jetsAk4Puppi";
  // myStrings["topjetLabel"] = myBools["isPuppi"]? "toppuppijets": "topjets";

  // h_jets = ctx.get_handle<std::vector<Jet>>(myStrings["jetLabel"]);
  // h_topjets = ctx.get_handle<std::vector<TopJet>>(myStrings["topjetLabel"]);

  // Set up histograms:

  book_histograms(ctx);
  // book_handles(ctx);
  PrintInputs();

  const MuonId muoId = AndId<Muon>(PtEtaCut(5, 2.4), MuonID(Muon::PFIsoTight));
  const ElectronId eleId = AndId<Electron>(ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_loose), PtEtaSCCut(5, 2.4));
  const JetId jetId_CHS = AndId<Jet> (JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(10, 5.2));
  const JetId jetId_Puppi = AndId<Jet> (JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(10, 5.2));

  // Set up selections
  PrimaryVertexId pvid = StandardPrimaryVertexId();
  modules.emplace_back(new PrimaryVertexCleaner(pvid));
  if(myBools["is_mc"]) {
    if(myBools["mclumiweight"])  weightsmodules.emplace_back(new MCLumiWeight(ctx));
    if(myBools["mcpileupreweight"]) weightsmodules.emplace_back(new MCPileupReweight(ctx,myStrings["SysType_PU"]));
  }
  if(myBools["eleid"]) modules.emplace_back(new ElectronCleaner(eleId));
  if(myBools["muid"])  modules.emplace_back(new MuonCleaner(muoId));


  jet_cleaner_CHS.reset( new GenericJetCleaner(ctx, myStrings["jetLabel_CHS"], false, jetId_CHS, jetId_CHS, muoId, eleId));
  jet_cleaner_Puppi.reset( new GenericJetCleaner(ctx, myStrings["jetLabel_Puppi"], false, jetId_Puppi, jetId_Puppi, muoId, eleId));

  h_jets_CHS = ctx.get_handle<vector<Jet>>(myStrings["jetLabel_CHS"]);
  h_jets_Puppi = ctx.get_handle<vector<Jet>>(myStrings["jetLabel_Puppi"]);

}


bool JERCStudiesModule::process(Event & event) {

  bool debug = false;

  if (debug) cout << "New event " << event.event << " " << event.rho << endl;

  fill_histograms(event, "nocuts");

  for(auto & m : weightsmodules) m->process(event);
  for(auto & m : modules) m->process(event);

  for (auto & jet : event.get(h_jets_CHS)) {
    if (debug) cout << "new Jet CHS -- eta = " << jet.eta() << " -- phi = " << jet.phi() << " pT_GT = " << jet.pt() << " pT_Original = ";
    jet.set_v4(jet.v4() * jet.JEC_factor_raw() );
    jet.set_JEC_factor_raw(1);
    if (debug) cout << jet.pt() << endl;
  }

  for (auto & jet : event.get(h_jets_Puppi)) {
    if (debug) cout << "new Jet Puppi -- eta = " << jet.eta() << " -- phi = " << jet.phi() << " pT_GT = " << jet.pt() << " pT_Original = ";
    jet.set_v4(jet.v4() * jet.JEC_factor_raw() );
    jet.set_JEC_factor_raw(1);
    if (debug) cout << jet.pt() << endl;
  }

  fill_histograms(event, "weights");

  jet_cleaner_CHS->process(event);
  jet_cleaner_Puppi->process(event);

  if (debug) {

    for (auto & jet : event.get(h_jets_CHS)) {
      cout << "new Jet CHS -- eta = " << jet.eta() << " -- phi = " << jet.phi() << " pT_Final = " << jet.pt() << endl;
    }

    for (auto & jet : event.get(h_jets_Puppi)) {
      cout << "new Jet Puppi -- eta = " << jet.eta() << " -- phi = " << jet.phi() << " pT_Final = " << jet.pt() << endl;
    }
  }

  if (debug) cout << "End event" << endl;

  fill_histograms(event, "cleaned");

  return false;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the JERCStudiesModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(JERCStudiesModule)
