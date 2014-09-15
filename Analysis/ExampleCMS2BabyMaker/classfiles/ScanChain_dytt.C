// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "Math/LorentzVector.h"

// DrellYan
#include "DrellYan.cc"

using namespace std;
using namespace dy;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1 /* , string skimFilePrefix = "test"*/) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  cout << "\n-------------------------------------------------------------------------------" << endl;

  TH1::SetDefaultSumw2();

  // Create histograms and take ownership of them
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  TH1D* h_ee_gen       = new TH1D("ee_gen",       "Dielectron mass (gen);m_{ee} (GeV/c^{2})", 60, 60, 120);
  TH1D* h_mm_gen       = new TH1D("mm_gen",       "Dimuon mass (gen);m_{#mu#mu} (GeV/c^{2})", 60, 60, 120);
  TH1D* h_yield_gen    = new TH1D("yields_gen",   "Yields (gen)", 4, 9.5, 13.5);
  TH1D* h_acc_num_gen  = new TH1D("acc_num_gen",  "Acceptance numerator (gen)", 4, 9.5, 13.5);
  TH1D* h_acc_den      = new TH1D("acc_den",      "Acceptance denominator", 4, 9.5, 13.5);
  TH1D* h_ee_reco      = new TH1D("ee_reco",      "Dielectron mass (reco);m_{ee} (GeV/c^{2})", 60, 60, 120);
  TH1D* h_mm_reco      = new TH1D("mm_reco",      "Dimuon mass (reco);m_{#mu#mu} (GeV/c^{2})", 60, 60, 120);
  TH1D* h_yield_reco   = new TH1D("yields_reco",  "Yields (reco)", 4, 9.5, 13.5);
  TH1D* h_acc_num_reco = new TH1D("acc_num_reco", "Acceptance numerator (reco)", 4, 9.5, 13.5);

  TH1D* h_accept_gen   = (TH1D*) h_acc_num_gen->Clone("accept_gen");
  TH1D* h_accept_reco  = (TH1D*) h_acc_num_reco->Clone("accept_reco");

  h_accept_gen->SetTitle(  "Acceptance (gen)"  );
  h_accept_reco->SetTitle( "Acceptance (reco)" );

  h_ee_gen->SetDirectory(rootdir);
  h_mm_gen->SetDirectory(rootdir);
  h_yield_gen->SetDirectory(rootdir);
  h_acc_num_gen->SetDirectory(rootdir);
  h_acc_den->SetDirectory(rootdir);
  h_ee_reco->SetDirectory(rootdir);
  h_mm_reco->SetDirectory(rootdir);
  h_yield_reco->SetDirectory(rootdir);
  h_acc_num_reco->SetDirectory(rootdir);
  h_accept_gen->SetDirectory(rootdir);
  h_accept_reco->SetDirectory(rootdir);

  // Holder for the sample name
  string samplename = "dytt";

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("tree");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    dy_obj.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      dy_obj.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      DrellYan::progress( nEventsTotal, nEventsChain );

	  //---------------------------------//
      // Analysis Code
	  //---------------------------------//

	  double int_lumi = 1.0;
	  if( !is_real_data() ) int_lumi = 0.082;

	  if( !is_gen_tt() ) continue;

	  //Fill appropriate dilepton mass histogram, and yield histogram
	  if( is_gen_ee() ) {
		h_ee_gen->Fill( gen_p4().M(), scale1fb()*int_lumi );
		h_yield_gen->Fill( 11, scale1fb()*int_lumi );
	  }
	  else if( is_gen_mm() ) {
		h_mm_gen->Fill( gen_p4().M(), scale1fb()*int_lumi );
		h_yield_gen->Fill( 13, scale1fb()*int_lumi );
	  }

	  if( is_reco_ee() ) {
		h_ee_reco->Fill( reco_p4().M(), scale1fb()*int_lumi );
		h_yield_reco->Fill( 11, scale1fb()*int_lumi );
	  }
	  else if( is_reco_mm() ) {
		h_mm_reco->Fill( reco_p4().M(), scale1fb()*int_lumi );
		h_yield_reco->Fill( 13, scale1fb()*int_lumi );
	  }

	  //Fill appropriate yield histograms
	  if( is_reco_acc_num() ) {
		h_acc_num_reco->Fill(abs(reco_lep1_id()), scale1fb()*int_lumi);
		h_acc_num_reco->Fill(10, scale1fb()*int_lumi);
	  }
	  if( is_gen_acc_num() ) {
		h_acc_num_gen->Fill(abs(gen_lep1_id()), scale1fb()*int_lumi);
		h_acc_num_gen->Fill(10, scale1fb()*int_lumi);
	  }
	  if( is_acc_den() ) {
		h_acc_den->Fill(abs(gen_lep1_id()), scale1fb()*int_lumi);
		h_acc_den->Fill(10, scale1fb()*int_lumi);
	  }

	  //---------------------------------//

    } //End of loop over file
  
    // Clean Up file and tree
    delete tree;
    file->Close();
    delete file;

  } //End of loop over files in chain

  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  // Calculate yields and acceptance
  double yielderr_ee_gen = 0;
  double yielderr_mm_gen = 0;
  double yielderr_total_gen = 0;
  double yield_ee_gen = 0;
  double yield_mm_gen = 0;
  double yield_total_gen = 0;
  double yielderr_ee_reco = 0;
  double yielderr_mm_reco = 0;
  double yielderr_total_reco = 0;
  double yield_ee_reco = 0;
  double yield_mm_reco = 0;
  double yield_total_reco = 0;

  double acc_ee_gen = 0;
  double acc_mm_gen = 0;
  double acc_total_gen = 0;
  double accerr_ee_gen = 0;
  double accerr_mm_gen = 0;
  double accerr_total_gen = 0;
  double acc_ee_reco = 0;
  double acc_mm_reco = 0;
  double acc_total_reco = 0;
  double accerr_ee_reco = 0;
  double accerr_mm_reco = 0;
  double accerr_total_reco = 0;

  yield_ee_gen     = h_ee_gen->IntegralAndError(     0, -1, yielderr_ee_gen     );
  yield_mm_gen     = h_mm_gen->IntegralAndError(     0, -1, yielderr_mm_gen     );
  yield_total_gen  = h_yield_gen->IntegralAndError(  0, -1, yielderr_total_gen  );
  yield_ee_reco    = h_ee_reco->IntegralAndError(    0, -1, yielderr_ee_reco    );
  yield_mm_reco    = h_mm_reco->IntegralAndError(    0, -1, yielderr_mm_reco    );
  yield_total_reco = h_yield_reco->IntegralAndError( 0, -1, yielderr_total_reco );

  h_accept_gen->Divide(  h_acc_num_gen,  h_acc_den, 1., 1., "B");
  h_accept_reco->Divide( h_acc_num_reco, h_acc_den, 1., 1., "B");

  acc_ee_gen         = h_accept_gen->GetBinContent(2);
  acc_mm_gen         = h_accept_gen->GetBinContent(4);
  acc_total_gen      = h_accept_gen->GetBinContent(1);
  accerr_ee_gen      = h_accept_gen->GetBinError(2);
  accerr_mm_gen      = h_accept_gen->GetBinError(4);
  accerr_total_gen   = h_accept_gen->GetBinError(1);
  acc_ee_reco        = h_accept_reco->GetBinContent(2);
  acc_mm_reco        = h_accept_reco->GetBinContent(4);
  acc_total_reco     = h_accept_reco->GetBinContent(1);
  accerr_ee_reco     = h_accept_reco->GetBinError(2);
  accerr_mm_reco     = h_accept_reco->GetBinError(4);
  accerr_total_reco  = h_accept_reco->GetBinError(1);

  //Print out yields
  std::cout << "Sample: " << samplename << std::endl;
  // std::cout << "Gen EE yield:\t" << yield_ee_gen << " +/- " << yielderr_ee_gen << std::endl;
  // std::cout << "Gen MM yield:\t" << yield_mm_gen << " +/- " << yielderr_mm_gen << std::endl;
  // std::cout << "Total gen yield:\t" << yield_total_gen << " +/- " << yielderr_total_gen << std::endl;
  // std::cout << "Reco EE yield:\t" << yield_ee_reco << " +/- " << yielderr_ee_reco << std::endl;
  // std::cout << "Reco MM yield:\t" << yield_mm_reco << " +/- " << yielderr_mm_reco << std::endl;
  // std::cout << "Total reco yield:\t" << yield_total_reco << " +/- " << yielderr_total_reco << std::endl;

  // if( samplename == "dyll" ) {
  // 	std::cout << "\nGen EE acceptance:\t" << acc_ee_gen << " +/- " << accerr_ee_gen << std::endl;
  // 	std::cout << "Gen MM acceptance:\t" << acc_mm_gen << " +/- " << accerr_mm_gen << std::endl;
  // 	std::cout << "Total gen acceptance:\t" << acc_total_gen << " +/- " << accerr_total_gen << std::endl;
  // 	std::cout << "Reco EE acceptance:\t" << acc_ee_reco << " +/- " << accerr_ee_reco << std::endl;
  // 	std::cout << "Reco MM acceptance:\t" << acc_mm_reco << " +/- " << accerr_mm_reco << std::endl;
  // 	std::cout << "Total reco acceptance:\t" << acc_total_reco << " +/- " << accerr_total_reco << std::endl;
  // }

  printf( "Gen EE yield:    \t%7.1f +/- %3.1f\n",  yield_ee_gen,     yielderr_ee_gen );
  printf( "Gen MM yield:    \t%7.1f +/- %3.1f\n",  yield_mm_gen,     yielderr_mm_gen );
  printf( "Total gen yield: \t%7.1f +/- %3.1f\n",  yield_total_gen,  yielderr_total_gen );
  printf( "Reco EE yield:   \t%7.1f +/- %3.1f\n",  yield_ee_reco,    yielderr_ee_reco );
  printf( "Reco MM yield:   \t%7.1f +/- %3.1f\n",  yield_mm_reco,    yielderr_mm_reco );
  printf( "Total reco yield:\t%7.1f +/- %3.1f\n",  yield_total_reco, yielderr_total_reco );

  if( samplename == "dyll" ) {
	printf( "\nGen EE acceptance:  \t%6.3f +/- %5.3f\n",  acc_ee_gen,     accerr_ee_gen );
	printf( "Gen MM acceptance:    \t%6.3f +/- %5.3f\n",  acc_mm_gen,     accerr_mm_gen );
	printf( "Total gen acceptance: \t%6.3f +/- %5.3f\n",  acc_total_gen,  accerr_total_gen );
	printf( "Reco EE acceptance:   \t%6.3f +/- %5.3f\n",  acc_ee_reco,    accerr_ee_reco );
	printf( "Reco MM acceptance:   \t%6.3f +/- %5.3f\n",  acc_mm_reco,    accerr_mm_reco );
	printf( "Total reco acceptance:\t%6.3f +/- %5.3f\n",  acc_total_reco, accerr_total_reco );
  }

  //Open output file
  TFile* outfile = new TFile(Form("../output/%s_plots.root", samplename.c_str() ), "RECREATE");

  //Save histograms to output file
  h_ee_gen->Write();
  h_mm_gen->Write();
  h_yield_gen->Write();
  h_acc_num_gen->Write();
  h_acc_den->Write();
  h_ee_reco->Write();
  h_mm_reco->Write();
  h_yield_reco->Write();
  h_acc_num_reco->Write();

  if( samplename == "dyll" ) {
	h_accept_gen->Write();
	h_accept_reco->Write();
  }

  outfile->Close();

  //Clean up histograms. Necessary because we own them, and ROOT won't do it for us
  delete h_ee_gen;
  delete h_mm_gen;
  delete h_yield_gen;
  delete h_acc_num_gen;
  delete h_acc_den;
  delete h_ee_reco;
  delete h_mm_reco;
  delete h_yield_reco;
  delete h_acc_num_reco;
  delete h_accept_gen;
  delete h_accept_reco;
  
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
