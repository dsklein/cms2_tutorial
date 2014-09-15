{

  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch_data     = new TChain("tree"); 
  TChain *ch_dyll     = new TChain("tree"); 
  TChain *ch_qcd      = new TChain("tree"); 
  TChain *ch_ttdil    = new TChain("tree"); 
  TChain *ch_tthad    = new TChain("tree"); 
  TChain *ch_ttslq    = new TChain("tree"); 
  TChain *ch_wjets    = new TChain("tree"); 
  TChain *ch_ww2l2nu  = new TChain("tree"); 
  TChain *ch_wz2l2q   = new TChain("tree"); 
  TChain *ch_wz3lnu   = new TChain("tree"); 
  TChain *ch_zz2l2nu  = new TChain("tree"); 
  TChain *ch_zz2l2q   = new TChain("tree"); 
  TChain *ch_zz4l     = new TChain("tree"); 


  ch_data->Add(    "babies/data.root");
  ch_dyll->Add(    "babies/dyll.root");
  ch_qcd->Add(     "babies/qcdmu15.root");
  ch_ttdil->Add(   "babies/ttdil.root");
  ch_tthad->Add(   "babies/tthad.root");
  ch_ttslq->Add(   "babies/ttslq.root");
  ch_wjets->Add(   "babies/wjets.root");
  ch_ww2l2nu->Add( "babies/ww2l2nu.root");
  ch_wz2l2q->Add(  "babies/wz2l2q.root");
  ch_wz3lnu->Add(  "babies/wz3lnu.root");
  ch_zz2l2nu->Add( "babies/zz2l2nu.root");
  ch_zz2l2q->Add(  "babies/zz2l2q.root");
  ch_zz4l->Add(    "babies/zz4l.root");


  ScanChain(ch_data   ); 
  ScanChain(ch_dyll   ); 
  ScanChain(ch_qcd    ); 
  ScanChain(ch_ttdil  ); 
  ScanChain(ch_tthad  ); 
  ScanChain(ch_ttslq  ); 
  ScanChain(ch_wjets  ); 
  ScanChain(ch_ww2l2nu); 
  ScanChain(ch_wz2l2q ); 
  ScanChain(ch_wz3lnu ); 
  ScanChain(ch_zz2l2nu); 
  ScanChain(ch_zz2l2q ); 
  ScanChain(ch_zz4l   ); 
}
