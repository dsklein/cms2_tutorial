{

  gROOT->ProcessLine(".L ScanChain_dytt.C+");

  TChain *ch = new TChain("tree"); 
  ch->Add("../babies/dyll.root");
  ScanChain(ch); 
}
