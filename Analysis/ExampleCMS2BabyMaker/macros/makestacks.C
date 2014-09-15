//ROOT includes
#include "TH1.h"
#include "TFile.h"
#include "THStack.h"


using std::cout;
using std::endl;


void makestacks() {

  TH1::SetDefaultSumw2();

  // -------------------------------------------------------------- //
  // Load up a whole bunch of histograms                            //
  // -------------------------------------------------------------- //

  cout << "Loading histograms..." << endl;

  // Open all our output files ///////////////////////////////////////////////////
  TFile* f_data    = TFile::Open("output/data_plots.root");
  TFile* f_dyll    = TFile::Open("output/dyll_plots.root");
  TFile* f_dytt    = TFile::Open("output/dytt_plots.root");
  TFile* f_wjets   = TFile::Open("output/wjets_plots.root");
  TFile* f_tt2l2v  = TFile::Open("output/ttdil_plots.root");
  TFile* f_ttlvjj  = TFile::Open("output/ttslq_plots.root");
  //TFile* f_tthadr  = TFile::Open("output/tthad_plots.root");
  //TFile* f_qcdmu   = TFile::Open("output/qcdmu15_plots.root");
  TFile* f_ww2l2v  = TFile::Open("output/ww2l2nu_plots.root");
  TFile* f_wz2l2q  = TFile::Open("output/wz2l2q_plots.root");
  TFile* f_wz3lv   = TFile::Open("output/wz3lnu_plots.root");
  TFile* f_zz2l2v  = TFile::Open("output/zz2l2nu_plots.root");
  TFile* f_zz2l2q  = TFile::Open("output/zz2l2q_plots.root");
  TFile* f_zz4l    = TFile::Open("output/zz4l_plots.root");

  // Get all the dilepton mass histograms from the output files /////////////////////////////////////////////
  //  Naming convention: e for electron, m for muon, g for gen, r for reco
  const TH1D* h_eg_dyll    = (TH1D*)(f_dyll->Get(   "ee_gen")->Clone("eg_dyll") );
  const TH1D* h_eg_wjets   = (TH1D*)(f_wjets->Get(  "ee_gen")->Clone("eg_wjets") );
  const TH1D* h_eg_tt2l2v  = (TH1D*)(f_tt2l2v->Get( "ee_gen")->Clone("eg_tt2l2v") );
  const TH1D* h_eg_ttlvjj  = (TH1D*)(f_ttlvjj->Get( "ee_gen")->Clone("eg_ttlvjj") );
  //const TH1D* h_eg_tthadr  = (TH1D*)(f_tthadr->Get( "ee_gen")->Clone("eg_tthadr") );
  //const TH1D* h_eg_qcdmu   = (TH1D*)(f_qcdmu->Get(  "ee_gen")->Clone("eg_qcdmu") );
  const TH1D* h_eg_ww2l2v  = (TH1D*)(f_ww2l2v->Get( "ee_gen")->Clone("eg_ww2l2v") );
  const TH1D* h_eg_wz2l2q  = (TH1D*)(f_wz2l2q->Get( "ee_gen")->Clone("eg_wz2l2q") );
  const TH1D* h_eg_wz3lv   = (TH1D*)(f_wz3lv->Get(  "ee_gen")->Clone("eg_wz3lv") );
  const TH1D* h_eg_zz2l2v  = (TH1D*)(f_zz2l2v->Get( "ee_gen")->Clone("eg_zz2l2v") );
  const TH1D* h_eg_zz2l2q  = (TH1D*)(f_zz2l2q->Get( "ee_gen")->Clone("eg_zz2l2q") );
  const TH1D* h_eg_zz4l    = (TH1D*)(f_zz4l->Get(   "ee_gen")->Clone("eg_zz4l") );

  const TH1D* h_mg_dyll    = (TH1D*)(f_dyll->Get(   "mm_gen")->Clone("mg_dyll") );
  const TH1D* h_mg_wjets   = (TH1D*)(f_wjets->Get(  "mm_gen")->Clone("mg_wjets") );
  const TH1D* h_mg_tt2l2v  = (TH1D*)(f_tt2l2v->Get( "mm_gen")->Clone("mg_tt2l2v") );
  const TH1D* h_mg_ttlvjj  = (TH1D*)(f_ttlvjj->Get( "mm_gen")->Clone("mg_ttlvjj") );
  //const TH1D* h_mg_tthadr  = (TH1D*)(f_tthadr->Get( "mm_gen")->Clone("mg_tthadr") );
  //const TH1D* h_mg_qcdmu   = (TH1D*)(f_qcdmu->Get(  "mm_gen")->Clone("mg_qcdmu") );
  const TH1D* h_mg_ww2l2v  = (TH1D*)(f_ww2l2v->Get( "mm_gen")->Clone("mg_ww2l2v") );
  const TH1D* h_mg_wz2l2q  = (TH1D*)(f_wz2l2q->Get( "mm_gen")->Clone("mg_ww2l2q") );
  const TH1D* h_mg_wz3lv   = (TH1D*)(f_wz3lv->Get(  "mm_gen")->Clone("mg_wz3lv") );
  const TH1D* h_mg_zz2l2v  = (TH1D*)(f_zz2l2v->Get( "mm_gen")->Clone("mg_zz2l2v") );
  const TH1D* h_mg_zz2l2q  = (TH1D*)(f_zz2l2q->Get( "mm_gen")->Clone("mg_zz2l2q") );
  const TH1D* h_mg_zz4l    = (TH1D*)(f_zz4l->Get(   "mm_gen")->Clone("mg_zz4l") );

  // Reco histograms also include data
  const TH1D* h_er_data    = (TH1D*)(f_data->Get(   "ee_reco")->Clone("er_data") );
  const TH1D* h_er_dyll    = (TH1D*)(f_dyll->Get(   "ee_reco")->Clone("er_dyll") );
  const TH1D* h_er_wjets   = (TH1D*)(f_wjets->Get(  "ee_reco")->Clone("er_wjets") );
  const TH1D* h_er_tt2l2v  = (TH1D*)(f_tt2l2v->Get( "ee_reco")->Clone("er_tt2l2v") );
  const TH1D* h_er_ttlvjj  = (TH1D*)(f_ttlvjj->Get( "ee_reco")->Clone("er_ttlvjj") );
  //const TH1D* h_er_tthadr  = (TH1D*)(f_tthadr->Get( "ee_reco")->Clone("er_tthadr") );
  //const TH1D* h_er_qcdmu   = (TH1D*)(f_qcdmu->Get(  "ee_reco")->Clone("er_qcdmu") );
  const TH1D* h_er_ww2l2v  = (TH1D*)(f_ww2l2v->Get( "ee_reco")->Clone("er_ww2l2v") );
  const TH1D* h_er_wz2l2q  = (TH1D*)(f_wz2l2q->Get( "ee_reco")->Clone("er_ww2l2q") );
  const TH1D* h_er_wz3lv   = (TH1D*)(f_wz3lv->Get(  "ee_reco")->Clone("er_wz3lv") );
  const TH1D* h_er_zz2l2v  = (TH1D*)(f_zz2l2v->Get( "ee_reco")->Clone("er_zz2l2v") );
  const TH1D* h_er_zz2l2q  = (TH1D*)(f_zz2l2q->Get( "ee_reco")->Clone("er_zz2l2q") );
  const TH1D* h_er_zz4l    = (TH1D*)(f_zz4l->Get(   "ee_reco")->Clone("er_zz4l") );

  const TH1D* h_mr_data    = (TH1D*)(f_data->Get(   "mm_reco")->Clone("mr_data") );
  const TH1D* h_mr_dyll    = (TH1D*)(f_dyll->Get(   "mm_reco")->Clone("mr_dyll") );
  const TH1D* h_mr_wjets   = (TH1D*)(f_wjets->Get(  "mm_reco")->Clone("mr_wjets") );
  const TH1D* h_mr_tt2l2v  = (TH1D*)(f_tt2l2v->Get( "mm_reco")->Clone("mr_tt2l2v") );
  const TH1D* h_mr_ttlvjj  = (TH1D*)(f_ttlvjj->Get( "mm_reco")->Clone("mr_ttlvjj") );
  //const TH1D* h_mr_tthadr  = (TH1D*)(f_tthadr->Get( "mm_reco")->Clone("mr_tthadr") );
  //const TH1D* h_mr_qcdmu   = (TH1D*)(f_qcdmu->Get(  "mm_reco")->Clone("mr_qcdmu") );
  const TH1D* h_mr_ww2l2v  = (TH1D*)(f_ww2l2v->Get( "mm_reco")->Clone("mr_ww2l2v") );
  const TH1D* h_mr_wz2l2q  = (TH1D*)(f_wz2l2q->Get( "mm_reco")->Clone("mr_wz2l2q") );
  const TH1D* h_mr_wz3lv   = (TH1D*)(f_wz3lv->Get(  "mm_reco")->Clone("mr_wz3lv") );
  const TH1D* h_mr_zz2l2v  = (TH1D*)(f_zz2l2v->Get( "mm_reco")->Clone("mr_zz2l2v") );
  const TH1D* h_mr_zz2l2q  = (TH1D*)(f_zz2l2q->Get( "mm_reco")->Clone("mr_zz2l2q") );
  const TH1D* h_mr_zz4l    = (TH1D*)(f_zz4l->Get(   "mm_reco")->Clone("mr_zz4l") );

  // Get yield histograms
  const TH1D* h_yield_data    = (TH1D*)(f_data->Get(   "yields_reco")->Clone("yield_data") );
  const TH1D* h_yield_dyll    = (TH1D*)(f_dyll->Get(   "yields_reco")->Clone("yield_dyll") );
  const TH1D* h_yield_dytt    = (TH1D*)(f_dytt->Get(   "yields_reco")->Clone("yield_dytt") );
  const TH1D* h_yield_wjets   = (TH1D*)(f_wjets->Get(  "yields_reco")->Clone("yield_wjets") );
  const TH1D* h_yield_tt2l2v  = (TH1D*)(f_tt2l2v->Get( "yields_reco")->Clone("yield_tt2l2v") );
  const TH1D* h_yield_ttlvjj  = (TH1D*)(f_ttlvjj->Get( "yields_reco")->Clone("yield_ttlvjj") );
  //const TH1D* h_yield_tthadr  = (TH1D*)(f_tthadr->Get( "yields_reco")->Clone("yield_tthadr") );
  //const TH1D* h_yield_qcdmu   = (TH1D*)(f_qcdmu->Get(  "yields_reco")->Clone("yield_qcdmu") );
  const TH1D* h_yield_ww2l2v  = (TH1D*)(f_ww2l2v->Get( "yields_reco")->Clone("yield_ww2l2v") );
  const TH1D* h_yield_wz2l2q  = (TH1D*)(f_wz2l2q->Get( "yields_reco")->Clone("yield_wz2l2q") );
  const TH1D* h_yield_wz3lv   = (TH1D*)(f_wz3lv->Get(  "yields_reco")->Clone("yield_wz3lv") );
  const TH1D* h_yield_zz2l2v  = (TH1D*)(f_zz2l2v->Get( "yields_reco")->Clone("yield_zz2l2v") );
  const TH1D* h_yield_zz2l2q  = (TH1D*)(f_zz2l2q->Get( "yields_reco")->Clone("yield_zz2l2q") );
  const TH1D* h_yield_zz4l    = (TH1D*)(f_zz4l->Get(   "yields_reco")->Clone("yield_zz4l") );

  // Get acceptance histogram
  const TH1D* h_accept = (TH1D*)(f_dyll->Get("accept_reco"));


  // -------------------------------------------------------------- //
  // Do some calculations for the yields and cross-sections         //
  // -------------------------------------------------------------- //

  cout << "Calculating cross-section...\n\n" << endl;

  h_yield_dyll->Add(h_yield_dytt, -1 ); //Don't double-count the dytt events! They should go in background MC, but not signal MC

  // Create background yield histogram
  TH1D* h_yield_bkg = h_yield_wjets->Clone("yield_bkg");
  h_yield_bkg->Add(h_yield_tt2l2v);
  h_yield_bkg->Add(h_yield_ttlvjj);
  //h_yield_bkg->Add(h_yield_tthadr);
  //h_yield_bkg->Add(h_yield_qcdmu);
  h_yield_bkg->Add(h_yield_ww2l2v);
  h_yield_bkg->Add(h_yield_wz2l2q);
  h_yield_bkg->Add(h_yield_wz3lv);
  h_yield_bkg->Add(h_yield_zz2l2v);
  h_yield_bkg->Add(h_yield_zz2l2q);
  h_yield_bkg->Add(h_yield_zz4l);
  h_yield_bkg->Add(h_yield_dytt);

  // Create total MC yield histogram (bkg + dyll)
  TH1D* h_yield_mc = h_yield_bkg->Clone("yield_mc");
  h_yield_mc->Add(h_yield_dyll);

  // Create signal histogram (observed - background)
  TH1D* h_yield_sig = h_yield_data->Clone("yield_sig");
  h_yield_sig->Add(h_yield_bkg, -1);

  // Calculate and print cross sections! ///////////////////
  TH1D* h_crosssec = h_yield_sig->Clone("crosssection");
  h_crosssec->Divide(h_crosssec, h_accept, 1, 1);
  h_crosssec->Scale( 1e-6/0.082 ); //Integrated luminosity = 0.082 fb^-1.  Multiply by 1e-6 to convert final answer to nanobarns

  printf("          \t|              ee\t|              mm\t|\n");
  printf("N_mc      \t| %7.1f +/- %3.1f\t| %7.1f +/- %3.1f\t|\n", h_yield_mc->GetBinContent(2),   h_yield_mc->GetBinError(2),   h_yield_mc->GetBinContent(4),   h_yield_mc->GetBinError(4)   );
  printf("N_dyll    \t| %7.1f +/- %3.1f\t| %7.1f +/- %3.1f\t|\n", h_yield_dyll->GetBinContent(2), h_yield_dyll->GetBinError(2), h_yield_dyll->GetBinContent(4), h_yield_dyll->GetBinError(4) );
  printf("N_bkg     \t| %7.1f +/- %3.1f\t| %7.1f +/- %3.1f\t|\n", h_yield_bkg->GetBinContent(2),  h_yield_bkg->GetBinError(2),  h_yield_bkg->GetBinContent(4),  h_yield_bkg->GetBinError(4)  );
  printf("N_obs     \t| %7.0f +/- %3.1f\t| %7.0f +/- %3.1f\t|\n", h_yield_data->GetBinContent(2), h_yield_data->GetBinError(2), h_yield_data->GetBinContent(4), h_yield_data->GetBinError(4) );
  printf("N_sig     \t| %7.1f +/- %3.1f\t| %7.1f +/- %3.1f\t|\n", h_yield_sig->GetBinContent(2),  h_yield_sig->GetBinError(2),  h_yield_sig->GetBinContent(4),  h_yield_sig->GetBinError(4)  );
  printf("Acc.      \t| %5.3f +/- %5.3f\t| %5.3f +/- %5.3f\t|\n", h_accept->GetBinContent(2),     h_accept->GetBinError(2),     h_accept->GetBinContent(4),     h_accept->GetBinError(4)     );
  printf("Sigma (nb)\t| %5.3f +/- %5.3f\t| %5.3f +/- %5.3f\t|\n", h_crosssec->GetBinContent(2),   h_crosssec->GetBinError(2),   h_crosssec->GetBinContent(4),   h_crosssec->GetBinError(4)   );

  cout << "\nMaking stacked histograms..." << endl;

  // -------------------------------------------------------------- //
  // Make pretty-looking stacked histograms                         //
  // -------------------------------------------------------------- //

  // Set histogram draw options (colors and markers) 
  h_er_data->SetMarkerStyle(20);
  h_mr_data->SetMarkerStyle(20);
  //Don't forget to make the marker size BIG!
  //Also, maybe rebin for ease of viewing?
  //Definitely have to pick the range carefully

  int color_dyll = kYellow-7;
  int color_wjets = kOrange-3;
  int color_tt2l2v = kRed+3;
  int color_ttlvjj = kRed-3;
  //int color_tthadr = kRed-9;
  int color_ww2l2v = kBlue+1;
  int color_wz2l2q = kAzure-2;
  int color_wz3lv = kCyan+1;
  int color_zz2l2v = kGreen+3;
  int color_zz2l2q = kGreen-3;
  //int color_qcdmu = kMagenta+3;
  int color_zz4l = kSpring-9;

  h_eg_dyll->SetFillColor( color_dyll );
  h_eg_wjets->SetFillColor( color_wjets);
  h_eg_tt2l2v->SetFillColor( color_tt2l2v );
  h_eg_ttlvjj->SetFillColor( color_ttlvjj );
  //h_eg_tthadr->SetFillColor( color_tthadr );
  h_eg_ww2l2v->SetFillColor( color_ww2l2v );  
  h_eg_wz2l2q->SetFillColor( color_wz2l2q );  
  h_eg_wz3lv->SetFillColor( color_wz3lv );   
  h_eg_zz2l2v->SetFillColor( color_zz2l2v );  
  h_eg_zz2l2q->SetFillColor( color_zz2l2q );  
  //h_eg_qcdmu->SetFillColor( color_qcdmu );
  h_eg_zz4l->SetFillColor( color_zz4l );    

  h_mg_dyll->SetFillColor( color_dyll );
  h_mg_wjets->SetFillColor( color_wjets );
  h_mg_tt2l2v->SetFillColor( color_tt2l2v );  
  h_mg_ttlvjj->SetFillColor( color_ttlvjj );
  //h_mg_tthadr->SetFillColor( color_tthadr );  
  h_mg_ww2l2v->SetFillColor( color_ww2l2v );  
  h_mg_wz2l2q->SetFillColor( color_wz2l2q );  
  h_mg_wz3lv->SetFillColor( color_wz3lv );   
  h_mg_zz2l2v->SetFillColor( color_zz2l2v );  
  h_mg_zz2l2q->SetFillColor( color_zz2l2q );
  //h_mg_qcdmu->SetFillColor( color_qcdmu );     
  h_mg_zz4l->SetFillColor( color_zz4l );    

  h_er_dyll->SetFillColor( color_dyll );
  h_er_wjets->SetFillColor( color_wjets );  
  h_er_tt2l2v->SetFillColor( color_tt2l2v ); 
  h_er_ttlvjj->SetFillColor( color_ttlvjj);
  //h_er_tthadr->SetFillColor( color_tthadr ); 
  h_er_ww2l2v->SetFillColor( color_ww2l2v ); 
  h_er_wz2l2q->SetFillColor( color_wz2l2q ); 
  h_er_wz3lv->SetFillColor( color_wz3lv );  
  h_er_zz2l2v->SetFillColor( color_zz2l2v ); 
  h_er_zz2l2q->SetFillColor( color_zz2l2q ); 
  //h_er_qcdmu->SetFillColor( color_qcdmu );  
  h_er_zz4l->SetFillColor( color_zz4l );   

  h_mr_dyll->SetFillColor( color_dyll );
  h_mr_wjets->SetFillColor( color_wjets );   
  h_mr_tt2l2v->SetFillColor( color_tt2l2v );  
  h_mr_ttlvjj->SetFillColor( color_ttlvjj );
  //h_mr_tthadr->SetFillColor( color_tthadr );  
  h_mr_ww2l2v->SetFillColor( color_ww2l2v );  
  h_mr_wz2l2q->SetFillColor( color_wz2l2q );  
  h_mr_wz3lv->SetFillColor( color_wz3lv );   
  h_mr_zz2l2v->SetFillColor( color_zz2l2v );  
  h_mr_zz2l2q->SetFillColor( color_zz2l2q );  
  //h_mr_qcdmu->SetFillColor( color_qcdmu );   
  h_mr_zz4l->SetFillColor( color_zz4l );    

  // Create the 4 different THStacks /////////////////////////////////////////////////////////
  THStack* eg_bkg = new THStack("eg_bkg", "Dielectron gen events;m_{ee} (GeV/c^{2});Events/2 GeV");
  THStack* mg_bkg = new THStack("mg_bkg", "Dimuon gen events;m_{#mu#mu} (GeV/c^{2});Events/2 GeV");
  THStack* er_bkg = new THStack("er_bkg", "Dielectron reco events;m_{ee} (GeV/c^{2});Events/2 GeV");
  THStack* mr_bkg = new THStack("mr_bkg", "Dimuon reco events;m_{#mu#mu} (GeV/c^{2});Events/2 GeV");

  // Add histograms to the stacks
  eg_bkg->Add(h_eg_zz4l);
  eg_bkg->Add(h_eg_zz2l2v);
  eg_bkg->Add(h_eg_ttlvjj);
  eg_bkg->Add(h_eg_ww2l2v);
  eg_bkg->Add(h_eg_wz3lv);
  eg_bkg->Add(h_eg_wz2l2q);
  eg_bkg->Add(h_eg_zz2l2q);
  eg_bkg->Add(h_eg_tt2l2v);
  //eg_bkg->Add(h_eg_tthadr);
  eg_bkg->Add(h_eg_wjets);
  //eg_bkg->Add(h_eg_qcdmu);
  eg_bkg->Add(h_eg_dyll);

  mg_bkg->Add(h_mg_zz4l);
  mg_bkg->Add(h_mg_zz2l2v);
  mg_bkg->Add(h_mg_ttlvjj);
  mg_bkg->Add(h_mg_ww2l2v);
  mg_bkg->Add(h_mg_wz3lv);
  mg_bkg->Add(h_mg_wz2l2q);
  mg_bkg->Add(h_mg_zz2l2q);
  mg_bkg->Add(h_mg_tt2l2v);
  //mg_bkg->Add(h_mg_tthadr);
  mg_bkg->Add(h_mg_wjets);
  //mg_bkg->Add(h_mg_qcdmu);
  mg_bkg->Add(h_mg_dyll);

  er_bkg->Add(h_er_zz4l);
  er_bkg->Add(h_er_zz2l2v);
  er_bkg->Add(h_er_ttlvjj);
  er_bkg->Add(h_er_ww2l2v);
  er_bkg->Add(h_er_wz3lv);
  er_bkg->Add(h_er_wz2l2q);
  er_bkg->Add(h_er_zz2l2q);
  er_bkg->Add(h_er_tt2l2v);
  //er_bkg->Add(h_er_tthadr);
  er_bkg->Add(h_er_wjets);
  //er_bkg->Add(h_er_qcdmu);
  er_bkg->Add(h_er_dyll);

  mr_bkg->Add(h_mr_zz4l);
  mr_bkg->Add(h_mr_zz2l2v);
  mr_bkg->Add(h_mr_ttlvjj);
  mr_bkg->Add(h_mr_ww2l2v);
  mr_bkg->Add(h_mr_wz3lv);
  mr_bkg->Add(h_mr_wz2l2q);
  mr_bkg->Add(h_mr_zz2l2q);
  mr_bkg->Add(h_mr_tt2l2v);
  //mr_bkg->Add(h_mr_tthadr);
  mr_bkg->Add(h_mr_wjets);
  //mr_bkg->Add(h_mr_qcdmu);
  mr_bkg->Add(h_mr_dyll);

  er_bkg->SetMinimum( 0.01 );
  mr_bkg->SetMinimum( 0.01 );
  er_bkg->SetMaximum( 3000 );
  mr_bkg->SetMaximum( 6000 );

  // Make legends
  TLegend* leg_eg = new TLegend(0.75, 0.42, 0.95, 0.97, "Sample");
  TLegend* leg_mg = new TLegend(0.75, 0.42, 0.95, 0.97, "Sample");
  TLegend* leg_er = new TLegend(0.75, 0.42, 0.95, 0.97, "Sample");
  TLegend* leg_mr = new TLegend(0.75, 0.42, 0.95, 0.97, "Sample");

  // Add histograms to legend
  leg_eg->AddEntry(h_eg_dyll, "DY #rightarrow 2l", "f");
  //leg_eg->AddEntry(h_eg_qcdmu, "QCD + #mu", "f");
  leg_eg->AddEntry(h_eg_wjets,"W + jets","f" );
  //leg_eg->AddEntry(h_eg_tthadr, "t#bar{t} #rightarrow hadrons", "f");
  leg_eg->AddEntry(h_eg_tt2l2v, "t#bar{t} #rightarrow 2l2#nu", "f");
  leg_eg->AddEntry(h_eg_zz2l2q, "ZZ #rightarrow 2l2q", "f");
  leg_eg->AddEntry(h_eg_wz2l2q, "WZ #rightarrow 2l2q", "f");
  leg_eg->AddEntry(h_eg_wz3lv, "WZ #rightarrow 3l#nu", "f");
  leg_eg->AddEntry(h_eg_ww2l2v, "WW #rightarrow 2l2#nu", "f");
  leg_eg->AddEntry(h_eg_ttlvjj, "t#bar{t} #rightarrow l#nu2j", "f");
  leg_eg->AddEntry(h_eg_zz2l2v, "ZZ #rightarrow 2l2#nu", "f");
  leg_eg->AddEntry(h_eg_zz4l, "ZZ #rightarrow 4l", "f");

  leg_mg->AddEntry(h_mg_dyll, "DY #rightarrow 2l", "f");
  //leg_mg->AddEntry(h_mg_qcdmu, "QCD + #mu", "f");
  leg_mg->AddEntry(h_mg_wjets, "W + jets","f");
  //leg_mg->AddEntry(h_mg_tthadr, "t#bar{t} #rightarrow hadrons", "f");
  leg_mg->AddEntry(h_mg_tt2l2v, "t#bar{t} #rightarrow 2l2#nu", "f");
  leg_mg->AddEntry(h_mg_zz2l2q, "ZZ #rightarrow 2l2q", "f");
  leg_mg->AddEntry(h_mg_wz2l2q, "WZ #rightarrow 2l2q", "f");
  leg_mg->AddEntry(h_mg_wz3lv, "WZ #rightarrow 3l#nu", "f");
  leg_mg->AddEntry(h_mg_ww2l2v, "WW #rightarrow 2l2#nu", "f");
  leg_mg->AddEntry(h_mg_ttlvjj, "t#bar{t} #rightarrow l#nu2j", "f");
  leg_mg->AddEntry(h_mg_zz2l2v, "ZZ #rightarrow 2l2#nu", "f");
  leg_mg->AddEntry(h_mg_zz4l, "ZZ #rightarrow 4l", "f");

  leg_er->AddEntry(h_er_data, "Data", "p");
  leg_er->AddEntry(h_er_dyll, "DY #rightarrow 2l", "f");
  //leg_er->AddEntry(h_er_qcdmu, "QCD + #mu", "f");
  leg_er->AddEntry(h_er_wjets, "W + jets","f");
  //leg_er->AddEntry(h_er_tthadr, "t#bar{t} #rightarrow hadrons", "f");
  leg_er->AddEntry(h_er_tt2l2v, "t#bar{t} #rightarrow 2l2#nu", "f");
  leg_er->AddEntry(h_er_zz2l2q, "ZZ #rightarrow 2l2q", "f");
  leg_er->AddEntry(h_er_wz2l2q, "WZ #rightarrow 2l2q", "f");
  leg_er->AddEntry(h_er_wz3lv, "WZ #rightarrow 3l#nu", "f");
  leg_er->AddEntry(h_er_ww2l2v, "WW #rightarrow 2l2#nu", "f");
  leg_er->AddEntry(h_er_ttlvjj, "t#bar{t} #rightarrow l#nu2j", "f");
  leg_er->AddEntry(h_er_zz2l2v, "ZZ #rightarrow 2l2#nu", "f");
  leg_er->AddEntry(h_er_zz4l, "ZZ #rightarrow 4l", "f");

  leg_mr->AddEntry(h_mr_data, "Data", "p");
  leg_mr->AddEntry(h_mr_dyll, "DY #rightarrow 2l", "f");
  //leg_mr->AddEntry(h_mr_qcdmu, "QCD + #mu", "f");
  leg_mr->AddEntry(h_mr_wjets, "W + jets","f");
  //leg_mr->AddEntry(h_mr_tthadr, "t#bar{t} #rightarrow hadrons", "f");
  leg_mr->AddEntry(h_mr_tt2l2v, "t#bar{t} #rightarrow 2l2#nu", "f");
  leg_mr->AddEntry(h_mr_zz2l2q, "ZZ #rightarrow 2l2q", "f");
  leg_mr->AddEntry(h_mr_wz2l2q, "WZ #rightarrow 2l2q", "f");
  leg_mr->AddEntry(h_mr_wz3lv, "WZ #rightarrow 3l#nu", "f");
  leg_mr->AddEntry(h_mr_ww2l2v, "WW #rightarrow 2l2#nu", "f");
  leg_mr->AddEntry(h_mr_ttlvjj, "t#bar{t} #rightarrow l#nu2j", "f");
  leg_mr->AddEntry(h_mr_zz2l2v, "ZZ #rightarrow 2l2#nu", "f");
  leg_mr->AddEntry(h_mr_zz4l, "ZZ #rightarrow 4l", "f");

  // Make canvases
  TCanvas* c_eg = new TCanvas("c_eg", "Dielectrons (gen)",  800, 600);
  TCanvas* c_mg = new TCanvas("C_mg", "Dimuons (gen)",      800, 600);
  TCanvas* c_er = new TCanvas("c_er", "Dielectrons (reco)", 800, 600);
  TCanvas* c_mr = new TCanvas("C_mr", "Dimuons (reco)",     800, 600);

  // Draw the histograms and legends!
  c_eg->cd();
  c_eg->SetLogy();
  eg_bkg->Draw("hist");
  leg_eg->Draw();

  c_mg->cd();
  c_mg->SetLogy();
  mg_bkg->Draw("hist");
  leg_mg->Draw();

  c_er->cd();
  c_er->SetLogy();
  er_bkg->Draw("hist");
  h_er_data->Draw("p,e,same");
  leg_er->Draw();

  c_mr->cd();
  c_mr->SetLogy();
  mr_bkg->Draw("hist");
  h_mr_data->Draw("p,e,same");
  leg_mr->Draw();

  c_eg->SaveAs("plots/stack_electron_gen.svg");
  c_mg->SaveAs("plots/stack_muon_gen.svg");
  c_er->SaveAs("plots/stack_electron_reco.svg");
  c_mr->SaveAs("plots/stack_muon_reco.svg");

  cout << "Done!" << endl;

}
