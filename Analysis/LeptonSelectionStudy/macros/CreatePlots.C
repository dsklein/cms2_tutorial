#include <iostream>
#include <string>
#include <map>
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLegend.h"
#include "THStack.h"
#include "HistTools.h"

// simple type to hold selection info
struct Selection
{
    std::string name;
    std::string title;
    TCut cut;
};

// a list of selections that are applied cumulatively
const std::vector<Selection>& GetMuonSelections()
{
    static std::vector<Selection> selections
    {
        {"base", "Base selections", ""  },

        {"fid" , "Passes Fiducial", "mu_is_fiducial"},
		{"pt25", "p_{T} > 25 GeV" , "p4.pt() > 25"  },

		{"globmu", "Is global muon", "mu_is_global_mu"},
		{"pfmu", "Is PF muon",    "mu_is_pf_mu" },
		{"chi2", "#chi^{2} / ndof < 10", "mu_gfit_chi2ndof < 10."},
		{"vstahits", ">0 valid station hits", "mu_nvalid_stahits > 0"},
		{"mstations", ">1 matched station", "mu_nmatched_stations > 1"},
		{"trklayrs", ">5 tracker layers", "mu_trk_nlayers > 5"},
		{"pxhits", ">0 valid pixel hits", "mu_trk_nvalid_pixhits > 0"},

		{"d0", "|d_{0,PV}| < 0.02 cm", "abs(d0) < 0.02"},
		{"dz", "|d_{z,PV}| < 0.5 cm", "abs(dz) < 0.5"},

		{"pfiso", "PFIso / pt < 0.12", "pfiso < 0.12"}



    };
    return selections;
}

// convenenice typedef for a simple histogram container
typedef std::map<std::string, TH1*> TH1Map;

// book the histograms
void BookHists(TH1Map& hm)
{
    TH1::SetDefaultSumw2(true);

    // muon cut-flow
    AddHist(hm, new TH1D("h_mu_cutflow", "Muon Selection Cut-Flow", GetMuonSelections().size(), 0, GetMuonSelections().size()));
    TH1& h_mu_cutflow = *hm["h_mu_cutflow"];
    h_mu_cutflow.GetXaxis()->SetLabelSize(0.05);
    h_mu_cutflow.SetStats(false);

    // muon plots
    for (size_t i = 0; i < GetMuonSelections().size(); i++)
    {
        const Selection& s = GetMuonSelections().at(i);  

        // book hists per selection
        AddHist(hm, new TH1D(Form("h_mu_count_%s"    , s.name.c_str()), Form("Muon count (%s);count"      , s.title.c_str()),   3,  0,   3));
        AddHist(hm, new TH1D(Form("h_mu_pt_%s"       , s.name.c_str()), Form("Muon p_{T} (%s);p_{T} (GeV)", s.title.c_str()), 50,  0, 100));
        AddHist(hm, new TH1D(Form("h_mu_eta_%s"      , s.name.c_str()), Form("Muon #eta (%s);#eta"        , s.title.c_str()), 30, 0,   3));
		//No histo for "is global" and "is PF"
		AddHist(hm, new TH1D(Form("h_mu_chi2_%s"     , s.name.c_str()), Form("Muon #chi^{2} / ndof (%s);#chi^{2} / ndof"    , s.title.c_str()), 30, 0,  25));
		AddHist(hm, new TH1D(Form("h_mu_vstahits_%s" , s.name.c_str()), Form("# of valid station hits (%s);Number of hits"  , s.title.c_str()), 52, -1.5,  50.5));
		AddHist(hm, new TH1D(Form("h_mu_mstations_%s" , s.name.c_str()), Form("# of matched stations (%s);Number of stations", s.title.c_str()), 12, -1.5,  10.5));
		AddHist(hm, new TH1D(Form("h_mu_trklayers_%s", s.name.c_str()), Form("# of tracker layers (%s);Number of layers"    , s.title.c_str()), 22, -1.5,  20.5));
		AddHist(hm, new TH1D(Form("h_mu_pxhits_%s"   , s.name.c_str()), Form("# of pixel hits (%s);Number of hits"          , s.title.c_str()), 17, -1.5,    15.5));
		AddHist(hm, new TH1D(Form("h_mu_d0_%s"       , s.name.c_str()), Form("Muon d_{0,PV} (%s);d_{0} (cm)"                , s.title.c_str()), 30, 0, 0.1));
		AddHist(hm, new TH1D(Form("h_mu_dz_%s"       , s.name.c_str()), Form("Muon d_{z,PV} (%s);d_{z} (cm)"                , s.title.c_str()), 30, 0,   2));
		AddHist(hm, new TH1D(Form("h_mu_pfiso_%s"    , s.name.c_str()), Form("Muon PFIso / p_{T} (%s);PFIso / p_{T}"        , s.title.c_str()), 50, 0,  0.5));

        // add cutflow bin
        h_mu_cutflow.GetXaxis()->SetBinLabel(i+1, s.name.c_str());
    }
}

void UnderOverFlow(TH1* hist) {

  const int nbinsx = hist->GetNbinsX();

  const double lowcenter = hist->GetBinCenter(1);
  const double highcenter = hist->GetBinCenter(nbinsx);

  const double underflow = hist->GetBinContent(0);
  const double overflow = hist->GetBinContent(nbinsx+1);

  const double undererr = hist->GetBinError(0);
  const double overerr = hist->GetBinError(nbinsx+1);
  const double lowerr = hist->GetBinError(1);
  const double higherr = hist->GetBinError(nbinsx);

  hist->Fill(lowcenter, underflow);
  hist->Fill(highcenter, overflow);

  hist->SetBinError( 1, sqrt(lowerr*lowerr + undererr*undererr) );
  hist->SetBinError( nbinsx, sqrt(higherr*higherr + overerr*overerr) );
}

// use TTree::Draw to create all the hists for each selection
void CreateMuonHists
(
    const std::string& sample,
    const bool require_prompt,
    const std::string& label
)
{
    // ntuple
    // ----------------------------- // 

    TChain chain("tree");
    chain.Add(Form("babies/%s_baby.root", sample.c_str()));

    // hists
    // ----------------------------- // 

    TH1Map hm;
    BookHists(hm);

    // Fill hists
    // ----------------------------- // 

    TCut selection = (require_prompt ? "is_prompt && is_mu" : "!is_prompt && is_mu");
	const double num_entries_orig = chain.GetEntries(selection);

    SetDirectory(hm, gDirectory);
    for (const auto& s : GetMuonSelections())
    {
        // cumulative selection
        selection = selection && s.cut;
        // normalizing to 100%
        //const TCut sel_scaled = selection;
		const double num_entries = chain.GetEntries(selection);
        const TCut sel_scaled = Form("%1.5e*(%s)", 100.0/num_entries, selection.GetTitle());
        const TCut sel_scaled_orig = Form("%1.5e*(%s)", 100.0/num_entries_orig, selection.GetTitle());

        std::cout << sel_scaled << std::endl;

        // fill the hists
        chain.Draw(Form("1>>h_mu_count_%s"                      , s.name.c_str()), sel_scaled_orig, "goff");
        chain.Draw(Form("p4.pt()>>h_mu_pt_%s"                   , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("abs(p4.eta())>>h_mu_eta_%s"            , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("mu_gfit_chi2ndof>>h_mu_chi2_%s"        , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("mu_nvalid_stahits>>h_mu_vstahits_%s"   , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("mu_nmatched_stations>>h_mu_mstations_%s", s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("mu_trk_nlayers>>h_mu_trklayers_%s"     , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("mu_trk_nvalid_pixhits>>h_mu_pxhits_%s" , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("abs(d0)>>h_mu_d0_%s"                   , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("abs(dz)>>h_mu_dz_%s"                   , s.name.c_str()), sel_scaled, "goff");
        chain.Draw(Form("pfiso>>h_mu_pfiso_%s"                  , s.name.c_str()), sel_scaled, "goff");

		TH1* h1 =  hm.at( Form("h_mu_pt_%s"       ,s.name.c_str()) );
		TH1* h2 =  hm.at( Form("h_mu_eta_%s"      ,s.name.c_str()) );
		TH1* h3 =  hm.at( Form("h_mu_chi2_%s"     ,s.name.c_str()) );
		TH1* h4 =  hm.at( Form("h_mu_vstahits_%s" ,s.name.c_str()) );
		TH1* h5 =  hm.at( Form("h_mu_mstations_%s" ,s.name.c_str()) );
		TH1* h6 =  hm.at( Form("h_mu_trklayers_%s",s.name.c_str()) );
		TH1* h7 =  hm.at( Form("h_mu_pxhits_%s"   ,s.name.c_str()) );
		TH1* h8 =  hm.at( Form("h_mu_d0_%s"       ,s.name.c_str()) );
		TH1* h9 =  hm.at( Form("h_mu_dz_%s"       ,s.name.c_str()) );
		TH1* h10 = hm.at( Form("h_mu_pfiso_%s"    ,s.name.c_str()) );

		UnderOverFlow(h1);
		UnderOverFlow(h2);
		UnderOverFlow(h3);
		UnderOverFlow(h4);
		UnderOverFlow(h5);
		UnderOverFlow(h6);
		UnderOverFlow(h7);
		UnderOverFlow(h8);
		UnderOverFlow(h9);
		UnderOverFlow(h10);
    }
    SetDirectory(hm, NULL);

    // fill cut-flow
    for (size_t i = 0; i < GetMuonSelections().size(); i++)
    {
        const Selection& s = GetMuonSelections().at(i);  
        double error = 0.0;
        const double value = hm["h_mu_count_"+s.name]->IntegralAndError(0, -1, error);
        hm["h_mu_cutflow"]->SetBinContent(i+1, value);
        hm["h_mu_cutflow"]->SetBinError  (i+1, error);
    }

    // write plots
    // ----------------------------- // 

    gSystem->Exec(Form("mkdir -p plots/%s", label.c_str()));
    SaveHists(hm, Form("plots/%s/%s_plots.root", label.c_str(), sample.c_str()));

    // cleanup
    ClearHists(hm);
}
    
// simple overlaying routine
void PrintOverlays
(
    const std::string label,
    const std::vector<Selection>& selections,
    TH1Map& hm_sig, 
    TH1Map& hm_bkg, 
    const std::string hist_stem,
    const bool log_plot = false,
    const std::string suffix = "png",
	const std::string drawoption = ""
)
{
    // canvas for drawing
    TCanvas c1;
	gStyle->SetOptStat("ioumen");

    for (const auto& s : selections)
    {
        // formatting
        TH1& h_sig = *hm_sig[hist_stem + "_" + s.name];
        h_sig.SetLineWidth(2);

        TH1& h_bkg = *hm_bkg[hist_stem + "_" + s.name];
        h_bkg.SetLineWidth(2);

        const std::string plot_name = Form("p_%s_%s", hist_stem.substr(2, std::string::npos).c_str(), s.name.c_str());

        // statbox
        SetStatBoxPosition(h_sig, 0.7, 0.7, 0.9, 0.9);
        SetStatBoxPosition(h_bkg, 0.7, 0.5, 0.9, 0.7);

        // legend
        TLegend legend(0.5, 0.9, 0.7, 0.75);
        legend.AddEntry(&h_sig, "Drell-Yan", "L");
        legend.AddEntry(&h_bkg, "Background"      , "L");
        legend.SetFillColor(0);  // 0 makes it the background clear on the pad
        legend.SetFillStyle(0);
        legend.SetBorderSize(0);

        // create stack
        THStack hs(plot_name.c_str(), h_sig.GetTitle());
        hs.Add(&h_sig);
        hs.Add(&h_bkg);
		const char* axistitle = h_bkg.GetXaxis()->GetTitle();

        // write
        hs.Draw( Form("nostack %s", drawoption.c_str()) );
        legend.Draw();
        c1.SetLogy(log_plot);
		hs.GetXaxis()->SetTitle(axistitle);
		gPad->Update();
        c1.Print(Form("plots/%s/overlays/%s/%s.%s", label.c_str(), suffix.c_str(), plot_name.c_str(), suffix.c_str()));
    }    
}

// simple overlaying routine for cutflow
void PrintCutflowOverlays
(
    const std::string label,
    TH1Map& hm_sig, 
    TH1Map& hm_bkg, 
    const bool log_plot = false,
    const std::string suffix = "png"
)
{
    // canvas for drawing
    TCanvas c1;

    // overlay cutflow
    TH1& h_sig = *hm_sig["h_mu_cutflow"]; h_sig.SetLineWidth(2);
    TH1& h_bkg = *hm_bkg["h_mu_cutflow"]; h_bkg.SetLineWidth(2);
    const std::string plot_name = "p_mu_cutflow"; 

    // legend
    TLegend legend(0.5, 0.9, 0.7, 0.75);
    legend.AddEntry(&h_sig, "Drell-Yan", "L");
    legend.AddEntry(&h_bkg, "Background"      , "L");
    legend.SetFillColor(0);  // 0 makes it the background clear on the pad
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);

    // create stack
    THStack hs(plot_name.c_str(), h_sig.GetTitle());
    hs.Add(&h_sig);
    hs.Add(&h_bkg);
    hs.SetMaximum(h_sig.GetMaximum()*1.5);

    // write
    hs.Draw("nostack");
    legend.Draw();
    c1.SetLogy(log_plot);
    c1.Print(Form("plots/%s/overlays/%s/%s.%s", label.c_str(), suffix.c_str(), plot_name.c_str(), suffix.c_str()));
}


// create the overlays for dyll and QCD samples
void CreateOverlays(const std::string& label, const std::string& suffix = "png")
{
    std::cout << "[CreatePlots] creating overlays plots" << std::endl;

    // get prompt plots
    TH1Map hm_sig;
    LoadHists(hm_sig, Form("plots/%s/dyll_plots.root", label.c_str()));
    SetLineColor(hm_sig, kBlue);
/*     NormalizeHists(hm_sig); */

    // get non prompt plots
    TH1Map hm_bkg;
    LoadHists(hm_bkg, Form("plots/%s/tthad_plots.root", label.c_str()));
    SetLineColor(hm_bkg, kRed);
/*     NormalizeHists(hm_bkg); */

    // create overlay
    gSystem->Exec(Form("mkdir -p plots/%s/overlays/%s", label.c_str(), suffix.c_str()));
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_pt" ,       /*log_plot=*/false, suffix, "");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_eta",       /*log_plot=*/false, suffix, "");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_chi2",      /*log_plot=*/true,  suffix, "");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_vstahits",  /*log_plot=*/false, suffix, "hist");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_mstations", /*log_plot=*/false, suffix, "hist");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_trklayers", /*log_plot=*/false, suffix, "hist");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_pxhits",    /*log_plot=*/false, suffix, "hist");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_d0",        /*log_plot=*/true,  suffix, "");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_dz",        /*log_plot=*/true,  suffix, "");
    PrintOverlays(label, GetMuonSelections(), hm_sig, hm_bkg, "h_mu_pfiso",     /*log_plot=*/true,  suffix, "");
    
    PrintCutflowOverlays(label, hm_sig, hm_bkg, /*log_plot=*/true, suffix);

    // cleanup
    ClearHists(hm_sig);
    ClearHists(hm_bkg);
}

// Do all the plots
void CreatePlots(const std::string& label = "test", const std::string& suffix = "png")
{
    std::cout << "[CreatePlots] creating dyll plots" << std::endl;
    CreateMuonHists("dyll", /*prompt=*/true , label);

    std::cout << "[CreatePlots] creating QCD plots" << std::endl;
    CreateMuonHists("tthad" , /*prompt=*/false, label);

    CreateOverlays(label, suffix);
}
