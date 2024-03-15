#include "TChain.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TF1.h"

#include <iostream>
#include <vector>

void MoveStat(TH1F * h,
	      const float x1, const float y1,
	      const float x2, const float y2,
	      const Color_t text_colour)
{
  h->GetListOfFunctions()->Print();
  TPaveStats * s = (TPaveStats*)h->GetListOfFunctions()->FindObject("stats");
  if(s != nullptr) {
    s->SetX1NDC(x1);
    s->SetY1NDC(y1);
    s->SetX2NDC(x2);
    s->SetY2NDC(y2);
    s->SetTextColor(text_colour);
  }

  h->GetFunction("gaus")->SetLineColor(text_colour);
}

TH1F * get_resolution(const char * filename,
		      const char * treename,
		      const char * var_opt,
		      const char * hname,
		      const char * binning,
		      const char * cut_opt,
		      const char * units,
		      const char * draw_opt,
		      const Color_t colour
		      )
{
  TChain * t = new TChain(treename);
  const int nfiles = t->Add(filename);
  if(nfiles < 1) {
    cerr << "Could not open tree: " << treename << " in files: " << filename << endl;
    exit(-1);
  }
  const int total_entries = t->GetEntries();
  const int cut_entries = t->GetEntries(cut_opt);
  cout << nfiles << " files from " << filename << " added to TChain " << treename << endl
       << "Chain has " << total_entries << " total entries" << endl
       << cut_entries << " after cut(s)" << endl
       << total_entries - cut_entries << " outside cuts = "
       << 100. * (total_entries - cut_entries) / total_entries << "%" << endl;
  
  TString full_var_exp = TString::Format("(%s)>>%s%s", var_opt, hname, binning);
  TString full_draw_opt = TString::Format("%s", draw_opt);
  cout << endl
       << "Drawing: " << full_var_exp << endl
       << "    Cut: " << cut_opt << endl
       << "   Draw: " << full_draw_opt << endl;
  const int nentries = t->Draw(full_var_exp, cut_opt, full_draw_opt);
  cout << "Plot has " << nentries << " entries" << endl;
  return nullptr;
  auto htemp = (TH1F*)gPad->GetPrimitive(hname); // 1D
  htemp->SetTitle(TString::Format(";%s %s;Number in bin", var_opt, units));

  htemp->SetLineWidth(2);
  htemp->SetLineColor(colour);
  htemp->SetMarkerColor(colour);
  cout << htemp->GetMarkerSize() << endl << endl;

  for(int ix = 1; ix < htemp->GetNbinsX(); ix++) {
    htemp->SetBinError(ix, TMath::Sqrt(htemp->GetBinContent(ix)));
  }

  return htemp;
}

void compare_resolution_one(float energy = 10,
			    const char * resolution = "true_energy-reco_energy",
			    const char * units = "(MeV)")
{
  gStyle->SetOptFit(0111);
  
  TCanvas * c = new TCanvas();
  c->SetTopMargin(0.2);
  c->SetLeftMargin(0.15);
  TH1F * h1 = get_resolution(TString::Format("%.1f/wcsim_output_%.1f_*_bonsai_flower_oldB_oldF.root", energy, energy), "lowEreco",
			     resolution, "Old Tune", "", "reco_energy > -98", units, "", kBlack);
  if(h1 == nullptr)
    return;
  TH1F * h2 = get_resolution(TString::Format("%.1f/wcsim_output_%.1f_*_bonsai_flower_oldB_oldF.root", energy, energy), "lowEreco",
			     resolution, "New Tune", "", "reco_energy > -98", units, "SAME", kRed);
  if(h2 == nullptr)
    return;
  c->BuildLegend(0.2, 0.8, 0.5, 1.0);

  h1->Fit("gaus");
  MoveStat(h1, 0.50, 0.8, 0.75, 1, kBlack);
  h2->Fit("gaus");
  MoveStat(h2, 0.75, 0.8, 1.00, 1,  kRed);
}

void compare_resolution()
{
  vector<float> energies = {10};
  
  std::vector<std::pair<string, string> > resolutions;
  resolutions.push_back(std::make_pair("true_energy-reco_energy", "MeV"));
  resolutions.push_back(std::make_pair("true_direction.Angle(reco_direction)", "degrees??"));
  resolutions.push_back(std::make_pair("(true_position-reco_position).Mag()", "cm??"));
  resolutions.push_back(std::make_pair("true_position.X()-reco_position.X()", "cm??"));
  resolutions.push_back(std::make_pair("true_position.Y()-reco_position.Y()", "cm??"));
  resolutions.push_back(std::make_pair("true_position.Z()-reco_position.Z()", "cm??"));
  resolutions.push_back(std::make_pair("true_time-reco_time", "s???"));
  
  for(const auto resolution : resolutions) {
    for(const auto energy : energies) {
      compare_resolution_one(energy, resolution.first.c_str(), resolution.second.c_str());
    }//energies
  }//resolutions

  cerr << "TODOOOOOOOOO" << endl
       << " Get direction & position magnitude plots working" << endl
       << " Get fit box working (I presume need to create a function)" << endl
       << " Fix units" << endl
       << " Make summary plots (as a function of energy)" << endl;
}
