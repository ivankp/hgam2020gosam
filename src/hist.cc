#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <filesystem>
#include <optional>
#include <memory>

#include <TFile.h>
#include <TKey.h>
#include <TChain.h>
#include <TEnv.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#define STR1(x) #x
#define STR(x) STR1(x)

#define TEST(var) std::cout << \
  "\033[33m" STR(__LINE__) ": " \
  "\033[36m" #var ":\033[0m " << (var) << std::endl;

#include "ivanp/error.hh"
#include "ivanp/branch_reader.hh"
#include "ivanp/enumerate.hh"
#include "ivanp/timed_counter.hh"
#include "ivanp/hist/histograms.hh"
#include "ivanp/hist/bins.hh"
#include "ivanp/hist/json.hh"

#include "reweighter.hh"
#include "Higgs2diphoton.hh"

using namespace ivanp::map::operators;

template <ivanp::map::Container C>
decltype(auto) operator+=(std::vector<auto,auto>& v, C&& r) {
  v.reserve(v.size()+ivanp::cont::size(r));
  std::forward<C>(r) | [&]<typename T>(T&& x){
    v.emplace_back(std::forward<T>(x));
  };
  return v;
}

using std::cout;
using std::cerr;
using std::endl;
using nlohmann::json;
namespace fj = fastjet;
using namespace ivanp;

using hist_t = ivanp::hist::histogram<
  ivanp::hist::nlo_mc_multibin,
  std::vector<std::vector<ivanp::hist::list_axis<std::vector<double>>>>,
  ivanp::hist::perbin_axes_spec<true>
>;
using bin_t = typename hist_t::bin_type;

bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cerr << "usage: " << argv[0] <<
      " config.json hists.json ntuple1.root ... \n";
    return 1;
  }

  const auto conf = json::parse(std::ifstream(argv[1]));
  cout << conf.dump(2) <<'\n'<< endl;

  // Chain input files
  std::unique_ptr<TChain> chain;
  { auto file = std::make_unique<TFile>(argv[3]);
    TTree* tree = nullptr;
    for (auto* _key : *file->GetListOfKeys()) {
      auto* key = static_cast<TKey*>(_key);
      const auto* key_class = TClass::GetClass(key->GetClassName(),true);
      if (!key_class) continue;
      if (key_class->InheritsFrom(TTree::Class())) {
        if (!tree)
          tree = dynamic_cast<TTree*>(key->ReadObj());
        else THROW("multiple trees in file \"",file->GetName(),"\"");
      }
    }
    chain = std::make_unique<TChain>(tree->GetName());
    cout << "Tree name: " << chain->GetName() << '\n';
  }
  for (int i=3; i<argc; ++i) {
    cout << argv[i] << endl;
    if (!chain->Add(argv[i],0)) THROW("failed to add file to chain");
  }
  cout << endl;

  // Read branches
  TTreeReader reader(&*chain);
  branch_reader<int>
    _id(reader,"id"),
    _nparticle(reader,"nparticle");
  branch_reader<double[],float[]>
    _px(reader,"px"),
    _py(reader,"py"),
    _pz(reader,"pz"),
    _E (reader,"E" );
  branch_reader<int[]> _kf(reader,"kf");
  branch_reader<double> _weight(reader,"weight");
  branch_reader<double> _weight2(reader,"weight2");

  std::optional<branch_reader<int>> _ncount;
  for (auto* b : *reader.GetTree()->GetListOfBranches()) {
    if (!strcmp(b->GetName(),"ncount")) {
      _ncount.emplace(reader,"ncount");
      break;
    }
  }

  std::vector<std::string> weights_names {
    "weight", "weight2"
  };

  // Make reweighters
  std::vector<reweighter> reweighters;
  // auto reweighters = conf.at("reweighting") | [&](const auto& def){
  //   return reweighter(reader,def);
  // };

  for (auto& rew : reweighters) // get weights names
    weights_names += rew.weights_names();

  // weight vector needs to be resized before histograms are created
  bin_t::weight.resize( weights_names.size() );

  // create histograms ----------------------------------------------
  std::map<const char*,hist_t*,chars_less> hists;
#define h_(NAME) \
  hist_t h_##NAME; \
  hists[#NAME] = &h_##NAME;

#include ".build/punch.hh" // defines cards_names

  // read histograms' binning
  for (const char* card_name : cards_names) {
    cout << card_name << '\n';
    const auto file_name = cat("punchcards/",card_name,".punch");
    auto get = [card=TEnv(file_name.c_str())]
      (const auto& key, const auto& x0) -> decltype(auto) {
        return card.GetValue(ivanp::cstr(key),x0);
      };

    hist_t::axes_type axes;

    std::string name, prefix;
    for (int j=0; ; ++j) {
      if (j>0) prefix = cat("Var",j,".");

      name = get(prefix+"VarName","");
      if (name.empty()) { if (j==0) continue; else break; }

      auto& dim = axes.emplace_back();

      auto binning = [&](const auto& suffix) {
        const char* val = get(cat(prefix,"Binning",suffix),"");
        if (*val=='\0') return false;
        auto& edges = dim.emplace_back().edges();
        std::stringstream ss(val);
        for (double x; ss >> x;) edges.emplace_back(x);
        if (edges.size()<2)
          THROW("must provide Binning with at least 2 edges");
        return true;
      };
      if (!binning(""))
        for (int i=1; binning(i); ++i);

      if (j==0) break;
    }
    if (axes.empty()) THROW("no binning in file \"",file_name,"\"");
    auto& h = *hists[card_name] = hist_t(std::move(axes));
    cout << json(h.axes()) << '\n';
  }
  cout << endl;

  std::vector<fj::PseudoJet> partons, jets;
  Higgs2diphoton higgs_decay(
    conf.value("higgs_decay_seed",Higgs2diphoton::seed_type(0)));
  Higgs2diphoton::photons_type photons;
  TLorentzVector higgs;

  const fj::JetDefinition jet_def(
    fj::antikt_algorithm, conf.value("jet_dR",0.4) );
  fj::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;

  const double
    jet_pt_cut = conf.value("jet_pt_cut",30.),
    jet_eta_cut = conf.value("jet_eta_cut",4.4);
  const bool apply_photon_cuts = conf.value("apply_photon_cuts",true);

  cout << endl;

  // EVENT LOOP =====================================================
  for (timed_counter cnt(reader.GetEntries()); reader.Next(); ++cnt) {
    const bool new_id = [id=*_id]{
      return (bin_t::id != id) ? ((bin_t::id = id),true) : false;
    }();

    // read 4-momenta -----------------------------------------------
    partons.clear();
    const unsigned np = *_nparticle;
    bool got_higgs = false;
    for (unsigned i=0; i<np; ++i) {
      if (_kf[i] == 25) {
        higgs.SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
        got_higgs = true;
      } else {
        partons.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
      }
    }
    if (!got_higgs) THROW("event without Higgs boson (kf==25)");

    // H -> γγ ----------------------------------------------------
    photons = higgs_decay(higgs,new_id);
    auto A_pT = photons | [](const auto& p){ return p.Pt(); };
    if (A_pT[0] < A_pT[1]) {
      std::swap(A_pT[0],A_pT[1]);
      std::swap(photons[0],photons[1]);
    }
    const auto A_eta = photons | [](const auto& p){ return p.Eta(); };

    if (apply_photon_cuts && (
      (A_pT[0] < 0.35*125.) or
      (A_pT[1] < 0.25*125.) or
      photon_eta_cut(std::abs(A_eta[0])) or
      photon_eta_cut(std::abs(A_eta[1]))
    )) continue;

    // Jets ---------------------------------------------------------
    jets = fj::ClusterSequence(partons,jet_def)
          .inclusive_jets(); // get clustered jets

    jets.erase( std::remove_if( jets.begin(), jets.end(), // apply jet cuts
      [=](const auto& jet){
        return (jet.pt() < jet_pt_cut)
        or (std::abs(jet.eta()) > jet_eta_cut);
      }), jets.end() );
    std::sort( jets.begin(), jets.end(), // sort by pT
      [](const auto& a, const auto& b){ return ( a.pt() > b.pt() ); });
    const unsigned njets = jets.size(); // number of clustered jets

    // set weights --------------------------------------------------
    bin_t::weight[0] = *_weight;
    bin_t::weight[1] = *_weight2;

    // Observables **************************************************
    const auto pT_yy = higgs.Pt();
    h_pT_yy(pT_yy);

    const auto yAbs_yy = std::abs(higgs.Rapidity());
    h_yAbs_yy_vs_pT_yy(yAbs_yy,pT_yy);

    h_N_j_30(njets);

  } // end event loop
  // ================================================================
  cout << endl;

  for (auto& [name,h] : hists)
    for (auto& bin : *h)
      bin.finalize();

  std::ofstream out(argv[2]);
  if (ends_with(argv[2],".cbor")) {
    auto cbor = json::to_cbor(hists);
    out.write(
      reinterpret_cast<const char*>(cbor.data()),
      cbor.size() * sizeof(decltype(cbor[0]))
    );
  } else {
    out << json(hists) << '\n';
  }
}
