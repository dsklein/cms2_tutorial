// -*- C++ -*-
#ifndef DrellYan_H
#define DrellYan_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
#include <string> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class DrellYan {
private: 
protected: 
	unsigned int index;
	int	run_;
	TBranch *run_branch;
	bool run_isLoaded;
	int	ls_;
	TBranch *ls_branch;
	bool ls_isLoaded;
	int	evt_;
	TBranch *evt_branch;
	bool evt_isLoaded;
	string *sample_;
	TBranch *sample_branch;
	bool sample_isLoaded;
	string *dataset_;
	TBranch *dataset_branch;
	bool dataset_isLoaded;
	string *filename_;
	TBranch *filename_branch;
	bool filename_isLoaded;
	bool	is_real_data_;
	TBranch *is_real_data_branch;
	bool is_real_data_isLoaded;
	double	scale1fb_;
	TBranch *scale1fb_branch;
	bool scale1fb_isLoaded;
	double	scale1fb_cms2_;
	TBranch *scale1fb_cms2_branch;
	bool scale1fb_cms2_isLoaded;
	double	lumi_;
	TBranch *lumi_branch;
	bool lumi_isLoaded;
	double	xsec_;
	TBranch *xsec_branch;
	bool xsec_isLoaded;
	double	filt_eff_;
	TBranch *filt_eff_branch;
	bool filt_eff_isLoaded;
	bool	is_gen_z_;
	TBranch *is_gen_z_branch;
	bool is_gen_z_isLoaded;
	bool	is_gen_tt_;
	TBranch *is_gen_tt_branch;
	bool is_gen_tt_isLoaded;
	bool	is_gen_ee_;
	TBranch *is_gen_ee_branch;
	bool is_gen_ee_isLoaded;
	bool	is_gen_mm_;
	TBranch *is_gen_mm_branch;
	bool is_gen_mm_isLoaded;
	bool	is_gen_acc_num_;
	TBranch *is_gen_acc_num_branch;
	bool is_gen_acc_num_isLoaded;
	bool	is_acc_den_;
	TBranch *is_acc_den_branch;
	bool is_acc_den_isLoaded;
	int	gen_lep1_id_;
	TBranch *gen_lep1_id_branch;
	bool gen_lep1_id_isLoaded;
	int	gen_lep1_charge_;
	TBranch *gen_lep1_charge_branch;
	bool gen_lep1_charge_isLoaded;
	int	gen_lep2_id_;
	TBranch *gen_lep2_id_branch;
	bool gen_lep2_id_isLoaded;
	int	gen_lep2_charge_;
	TBranch *gen_lep2_charge_branch;
	bool gen_lep2_charge_isLoaded;
	bool	is_reco_ee_;
	TBranch *is_reco_ee_branch;
	bool is_reco_ee_isLoaded;
	bool	is_reco_mm_;
	TBranch *is_reco_mm_branch;
	bool is_reco_mm_isLoaded;
	bool	is_reco_acc_num_;
	TBranch *is_reco_acc_num_branch;
	bool is_reco_acc_num_isLoaded;
	int	reco_lep1_id_;
	TBranch *reco_lep1_id_branch;
	bool reco_lep1_id_isLoaded;
	int	reco_lep1_charge_;
	TBranch *reco_lep1_charge_branch;
	bool reco_lep1_charge_isLoaded;
	int	reco_lep2_id_;
	TBranch *reco_lep2_id_branch;
	bool reco_lep2_id_isLoaded;
	int	reco_lep2_charge_;
	TBranch *reco_lep2_charge_branch;
	bool reco_lep2_charge_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *gen_p4_;
	TBranch *gen_p4_branch;
	bool gen_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *gen_lep1_p4_;
	TBranch *gen_lep1_p4_branch;
	bool gen_lep1_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *gen_lep2_p4_;
	TBranch *gen_lep2_p4_branch;
	bool gen_lep2_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *reco_p4_;
	TBranch *reco_p4_branch;
	bool reco_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *reco_lep1_p4_;
	TBranch *reco_lep1_p4_branch;
	bool reco_lep1_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > *reco_lep2_p4_;
	TBranch *reco_lep2_p4_branch;
	bool reco_lep2_p4_isLoaded;
public: 
void Init(TTree *tree) {
	gen_p4_branch = 0;
	if (tree->GetBranch("gen_p4") != 0) {
		gen_p4_branch = tree->GetBranch("gen_p4");
		if (gen_p4_branch) {gen_p4_branch->SetAddress(&gen_p4_);}
	}
	gen_lep1_p4_branch = 0;
	if (tree->GetBranch("gen_lep1_p4") != 0) {
		gen_lep1_p4_branch = tree->GetBranch("gen_lep1_p4");
		if (gen_lep1_p4_branch) {gen_lep1_p4_branch->SetAddress(&gen_lep1_p4_);}
	}
	gen_lep2_p4_branch = 0;
	if (tree->GetBranch("gen_lep2_p4") != 0) {
		gen_lep2_p4_branch = tree->GetBranch("gen_lep2_p4");
		if (gen_lep2_p4_branch) {gen_lep2_p4_branch->SetAddress(&gen_lep2_p4_);}
	}
	reco_p4_branch = 0;
	if (tree->GetBranch("reco_p4") != 0) {
		reco_p4_branch = tree->GetBranch("reco_p4");
		if (reco_p4_branch) {reco_p4_branch->SetAddress(&reco_p4_);}
	}
	reco_lep1_p4_branch = 0;
	if (tree->GetBranch("reco_lep1_p4") != 0) {
		reco_lep1_p4_branch = tree->GetBranch("reco_lep1_p4");
		if (reco_lep1_p4_branch) {reco_lep1_p4_branch->SetAddress(&reco_lep1_p4_);}
	}
	reco_lep2_p4_branch = 0;
	if (tree->GetBranch("reco_lep2_p4") != 0) {
		reco_lep2_p4_branch = tree->GetBranch("reco_lep2_p4");
		if (reco_lep2_p4_branch) {reco_lep2_p4_branch->SetAddress(&reco_lep2_p4_);}
	}
  tree->SetMakeClass(1);
	run_branch = 0;
	if (tree->GetBranch("run") != 0) {
		run_branch = tree->GetBranch("run");
		if (run_branch) {run_branch->SetAddress(&run_);}
	}
	ls_branch = 0;
	if (tree->GetBranch("ls") != 0) {
		ls_branch = tree->GetBranch("ls");
		if (ls_branch) {ls_branch->SetAddress(&ls_);}
	}
	evt_branch = 0;
	if (tree->GetBranch("evt") != 0) {
		evt_branch = tree->GetBranch("evt");
		if (evt_branch) {evt_branch->SetAddress(&evt_);}
	}
	sample_branch = 0;
	if (tree->GetBranch("sample") != 0) {
		sample_branch = tree->GetBranch("sample");
		if (sample_branch) {sample_branch->SetAddress(&sample_);}
	}
	dataset_branch = 0;
	if (tree->GetBranch("dataset") != 0) {
		dataset_branch = tree->GetBranch("dataset");
		if (dataset_branch) {dataset_branch->SetAddress(&dataset_);}
	}
	filename_branch = 0;
	if (tree->GetBranch("filename") != 0) {
		filename_branch = tree->GetBranch("filename");
		if (filename_branch) {filename_branch->SetAddress(&filename_);}
	}
	is_real_data_branch = 0;
	if (tree->GetBranch("is_real_data") != 0) {
		is_real_data_branch = tree->GetBranch("is_real_data");
		if (is_real_data_branch) {is_real_data_branch->SetAddress(&is_real_data_);}
	}
	scale1fb_branch = 0;
	if (tree->GetBranch("scale1fb") != 0) {
		scale1fb_branch = tree->GetBranch("scale1fb");
		if (scale1fb_branch) {scale1fb_branch->SetAddress(&scale1fb_);}
	}
	scale1fb_cms2_branch = 0;
	if (tree->GetBranch("scale1fb_cms2") != 0) {
		scale1fb_cms2_branch = tree->GetBranch("scale1fb_cms2");
		if (scale1fb_cms2_branch) {scale1fb_cms2_branch->SetAddress(&scale1fb_cms2_);}
	}
	lumi_branch = 0;
	if (tree->GetBranch("lumi") != 0) {
		lumi_branch = tree->GetBranch("lumi");
		if (lumi_branch) {lumi_branch->SetAddress(&lumi_);}
	}
	xsec_branch = 0;
	if (tree->GetBranch("xsec") != 0) {
		xsec_branch = tree->GetBranch("xsec");
		if (xsec_branch) {xsec_branch->SetAddress(&xsec_);}
	}
	filt_eff_branch = 0;
	if (tree->GetBranch("filt_eff") != 0) {
		filt_eff_branch = tree->GetBranch("filt_eff");
		if (filt_eff_branch) {filt_eff_branch->SetAddress(&filt_eff_);}
	}
	is_gen_z_branch = 0;
	if (tree->GetBranch("is_gen_z") != 0) {
		is_gen_z_branch = tree->GetBranch("is_gen_z");
		if (is_gen_z_branch) {is_gen_z_branch->SetAddress(&is_gen_z_);}
	}
	is_gen_tt_branch = 0;
	if (tree->GetBranch("is_gen_tt") != 0) {
		is_gen_tt_branch = tree->GetBranch("is_gen_tt");
		if (is_gen_tt_branch) {is_gen_tt_branch->SetAddress(&is_gen_tt_);}
	}
	is_gen_ee_branch = 0;
	if (tree->GetBranch("is_gen_ee") != 0) {
		is_gen_ee_branch = tree->GetBranch("is_gen_ee");
		if (is_gen_ee_branch) {is_gen_ee_branch->SetAddress(&is_gen_ee_);}
	}
	is_gen_mm_branch = 0;
	if (tree->GetBranch("is_gen_mm") != 0) {
		is_gen_mm_branch = tree->GetBranch("is_gen_mm");
		if (is_gen_mm_branch) {is_gen_mm_branch->SetAddress(&is_gen_mm_);}
	}
	is_gen_acc_num_branch = 0;
	if (tree->GetBranch("is_gen_acc_num") != 0) {
		is_gen_acc_num_branch = tree->GetBranch("is_gen_acc_num");
		if (is_gen_acc_num_branch) {is_gen_acc_num_branch->SetAddress(&is_gen_acc_num_);}
	}
	is_acc_den_branch = 0;
	if (tree->GetBranch("is_acc_den") != 0) {
		is_acc_den_branch = tree->GetBranch("is_acc_den");
		if (is_acc_den_branch) {is_acc_den_branch->SetAddress(&is_acc_den_);}
	}
	gen_lep1_id_branch = 0;
	if (tree->GetBranch("gen_lep1_id") != 0) {
		gen_lep1_id_branch = tree->GetBranch("gen_lep1_id");
		if (gen_lep1_id_branch) {gen_lep1_id_branch->SetAddress(&gen_lep1_id_);}
	}
	gen_lep1_charge_branch = 0;
	if (tree->GetBranch("gen_lep1_charge") != 0) {
		gen_lep1_charge_branch = tree->GetBranch("gen_lep1_charge");
		if (gen_lep1_charge_branch) {gen_lep1_charge_branch->SetAddress(&gen_lep1_charge_);}
	}
	gen_lep2_id_branch = 0;
	if (tree->GetBranch("gen_lep2_id") != 0) {
		gen_lep2_id_branch = tree->GetBranch("gen_lep2_id");
		if (gen_lep2_id_branch) {gen_lep2_id_branch->SetAddress(&gen_lep2_id_);}
	}
	gen_lep2_charge_branch = 0;
	if (tree->GetBranch("gen_lep2_charge") != 0) {
		gen_lep2_charge_branch = tree->GetBranch("gen_lep2_charge");
		if (gen_lep2_charge_branch) {gen_lep2_charge_branch->SetAddress(&gen_lep2_charge_);}
	}
	is_reco_ee_branch = 0;
	if (tree->GetBranch("is_reco_ee") != 0) {
		is_reco_ee_branch = tree->GetBranch("is_reco_ee");
		if (is_reco_ee_branch) {is_reco_ee_branch->SetAddress(&is_reco_ee_);}
	}
	is_reco_mm_branch = 0;
	if (tree->GetBranch("is_reco_mm") != 0) {
		is_reco_mm_branch = tree->GetBranch("is_reco_mm");
		if (is_reco_mm_branch) {is_reco_mm_branch->SetAddress(&is_reco_mm_);}
	}
	is_reco_acc_num_branch = 0;
	if (tree->GetBranch("is_reco_acc_num") != 0) {
		is_reco_acc_num_branch = tree->GetBranch("is_reco_acc_num");
		if (is_reco_acc_num_branch) {is_reco_acc_num_branch->SetAddress(&is_reco_acc_num_);}
	}
	reco_lep1_id_branch = 0;
	if (tree->GetBranch("reco_lep1_id") != 0) {
		reco_lep1_id_branch = tree->GetBranch("reco_lep1_id");
		if (reco_lep1_id_branch) {reco_lep1_id_branch->SetAddress(&reco_lep1_id_);}
	}
	reco_lep1_charge_branch = 0;
	if (tree->GetBranch("reco_lep1_charge") != 0) {
		reco_lep1_charge_branch = tree->GetBranch("reco_lep1_charge");
		if (reco_lep1_charge_branch) {reco_lep1_charge_branch->SetAddress(&reco_lep1_charge_);}
	}
	reco_lep2_id_branch = 0;
	if (tree->GetBranch("reco_lep2_id") != 0) {
		reco_lep2_id_branch = tree->GetBranch("reco_lep2_id");
		if (reco_lep2_id_branch) {reco_lep2_id_branch->SetAddress(&reco_lep2_id_);}
	}
	reco_lep2_charge_branch = 0;
	if (tree->GetBranch("reco_lep2_charge") != 0) {
		reco_lep2_charge_branch = tree->GetBranch("reco_lep2_charge");
		if (reco_lep2_charge_branch) {reco_lep2_charge_branch->SetAddress(&reco_lep2_charge_);}
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		run_isLoaded = false;
		ls_isLoaded = false;
		evt_isLoaded = false;
		sample_isLoaded = false;
		dataset_isLoaded = false;
		filename_isLoaded = false;
		is_real_data_isLoaded = false;
		scale1fb_isLoaded = false;
		scale1fb_cms2_isLoaded = false;
		lumi_isLoaded = false;
		xsec_isLoaded = false;
		filt_eff_isLoaded = false;
		is_gen_z_isLoaded = false;
		is_gen_tt_isLoaded = false;
		is_gen_ee_isLoaded = false;
		is_gen_mm_isLoaded = false;
		is_gen_acc_num_isLoaded = false;
		is_acc_den_isLoaded = false;
		gen_lep1_id_isLoaded = false;
		gen_lep1_charge_isLoaded = false;
		gen_lep2_id_isLoaded = false;
		gen_lep2_charge_isLoaded = false;
		is_reco_ee_isLoaded = false;
		is_reco_mm_isLoaded = false;
		is_reco_acc_num_isLoaded = false;
		reco_lep1_id_isLoaded = false;
		reco_lep1_charge_isLoaded = false;
		reco_lep2_id_isLoaded = false;
		reco_lep2_charge_isLoaded = false;
		gen_p4_isLoaded = false;
		gen_lep1_p4_isLoaded = false;
		gen_lep2_p4_isLoaded = false;
		reco_p4_isLoaded = false;
		reco_lep1_p4_isLoaded = false;
		reco_lep2_p4_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (run_branch != 0) run();
	if (ls_branch != 0) ls();
	if (evt_branch != 0) evt();
	if (sample_branch != 0) sample();
	if (dataset_branch != 0) dataset();
	if (filename_branch != 0) filename();
	if (is_real_data_branch != 0) is_real_data();
	if (scale1fb_branch != 0) scale1fb();
	if (scale1fb_cms2_branch != 0) scale1fb_cms2();
	if (lumi_branch != 0) lumi();
	if (xsec_branch != 0) xsec();
	if (filt_eff_branch != 0) filt_eff();
	if (is_gen_z_branch != 0) is_gen_z();
	if (is_gen_tt_branch != 0) is_gen_tt();
	if (is_gen_ee_branch != 0) is_gen_ee();
	if (is_gen_mm_branch != 0) is_gen_mm();
	if (is_gen_acc_num_branch != 0) is_gen_acc_num();
	if (is_acc_den_branch != 0) is_acc_den();
	if (gen_lep1_id_branch != 0) gen_lep1_id();
	if (gen_lep1_charge_branch != 0) gen_lep1_charge();
	if (gen_lep2_id_branch != 0) gen_lep2_id();
	if (gen_lep2_charge_branch != 0) gen_lep2_charge();
	if (is_reco_ee_branch != 0) is_reco_ee();
	if (is_reco_mm_branch != 0) is_reco_mm();
	if (is_reco_acc_num_branch != 0) is_reco_acc_num();
	if (reco_lep1_id_branch != 0) reco_lep1_id();
	if (reco_lep1_charge_branch != 0) reco_lep1_charge();
	if (reco_lep2_id_branch != 0) reco_lep2_id();
	if (reco_lep2_charge_branch != 0) reco_lep2_charge();
	if (gen_p4_branch != 0) gen_p4();
	if (gen_lep1_p4_branch != 0) gen_lep1_p4();
	if (gen_lep2_p4_branch != 0) gen_lep2_p4();
	if (reco_p4_branch != 0) reco_p4();
	if (reco_lep1_p4_branch != 0) reco_lep1_p4();
	if (reco_lep2_p4_branch != 0) reco_lep2_p4();
}

	int &run()
	{
		if (not run_isLoaded) {
			if (run_branch != 0) {
				run_branch->GetEntry(index);
			} else { 
				printf("branch run_branch does not exist!\n");
				exit(1);
			}
			run_isLoaded = true;
		}
		return run_;
	}
	int &ls()
	{
		if (not ls_isLoaded) {
			if (ls_branch != 0) {
				ls_branch->GetEntry(index);
			} else { 
				printf("branch ls_branch does not exist!\n");
				exit(1);
			}
			ls_isLoaded = true;
		}
		return ls_;
	}
	int &evt()
	{
		if (not evt_isLoaded) {
			if (evt_branch != 0) {
				evt_branch->GetEntry(index);
			} else { 
				printf("branch evt_branch does not exist!\n");
				exit(1);
			}
			evt_isLoaded = true;
		}
		return evt_;
	}
	string &sample()
	{
		if (not sample_isLoaded) {
			if (sample_branch != 0) {
				sample_branch->GetEntry(index);
			} else { 
				printf("branch sample_branch does not exist!\n");
				exit(1);
			}
			sample_isLoaded = true;
		}
		return *sample_;
	}
	string &dataset()
	{
		if (not dataset_isLoaded) {
			if (dataset_branch != 0) {
				dataset_branch->GetEntry(index);
			} else { 
				printf("branch dataset_branch does not exist!\n");
				exit(1);
			}
			dataset_isLoaded = true;
		}
		return *dataset_;
	}
	string &filename()
	{
		if (not filename_isLoaded) {
			if (filename_branch != 0) {
				filename_branch->GetEntry(index);
			} else { 
				printf("branch filename_branch does not exist!\n");
				exit(1);
			}
			filename_isLoaded = true;
		}
		return *filename_;
	}
	bool &	is_real_data()
	{
		if (not is_real_data_isLoaded) {
			if (is_real_data_branch != 0) {
				is_real_data_branch->GetEntry(index);
			} else { 
				printf("branch is_real_data_branch does not exist!\n");
				exit(1);
			}
			is_real_data_isLoaded = true;
		}
		return is_real_data_;
	}
	double &	scale1fb()
	{
		if (not scale1fb_isLoaded) {
			if (scale1fb_branch != 0) {
				scale1fb_branch->GetEntry(index);
			} else { 
				printf("branch scale1fb_branch does not exist!\n");
				exit(1);
			}
			scale1fb_isLoaded = true;
		}
		return scale1fb_;
	}
	double &	scale1fb_cms2()
	{
		if (not scale1fb_cms2_isLoaded) {
			if (scale1fb_cms2_branch != 0) {
				scale1fb_cms2_branch->GetEntry(index);
			} else { 
				printf("branch scale1fb_cms2_branch does not exist!\n");
				exit(1);
			}
			scale1fb_cms2_isLoaded = true;
		}
		return scale1fb_cms2_;
	}
	double &	lumi()
	{
		if (not lumi_isLoaded) {
			if (lumi_branch != 0) {
				lumi_branch->GetEntry(index);
			} else { 
				printf("branch lumi_branch does not exist!\n");
				exit(1);
			}
			lumi_isLoaded = true;
		}
		return lumi_;
	}
	double &	xsec()
	{
		if (not xsec_isLoaded) {
			if (xsec_branch != 0) {
				xsec_branch->GetEntry(index);
			} else { 
				printf("branch xsec_branch does not exist!\n");
				exit(1);
			}
			xsec_isLoaded = true;
		}
		return xsec_;
	}
	double &	filt_eff()
	{
		if (not filt_eff_isLoaded) {
			if (filt_eff_branch != 0) {
				filt_eff_branch->GetEntry(index);
			} else { 
				printf("branch filt_eff_branch does not exist!\n");
				exit(1);
			}
			filt_eff_isLoaded = true;
		}
		return filt_eff_;
	}
	bool &	is_gen_z()
	{
		if (not is_gen_z_isLoaded) {
			if (is_gen_z_branch != 0) {
				is_gen_z_branch->GetEntry(index);
			} else { 
				printf("branch is_gen_z_branch does not exist!\n");
				exit(1);
			}
			is_gen_z_isLoaded = true;
		}
		return is_gen_z_;
	}
	bool &	is_gen_tt()
	{
		if (not is_gen_tt_isLoaded) {
			if (is_gen_tt_branch != 0) {
				is_gen_tt_branch->GetEntry(index);
			} else { 
				printf("branch is_gen_tt_branch does not exist!\n");
				exit(1);
			}
			is_gen_tt_isLoaded = true;
		}
		return is_gen_tt_;
	}
	bool &	is_gen_ee()
	{
		if (not is_gen_ee_isLoaded) {
			if (is_gen_ee_branch != 0) {
				is_gen_ee_branch->GetEntry(index);
			} else { 
				printf("branch is_gen_ee_branch does not exist!\n");
				exit(1);
			}
			is_gen_ee_isLoaded = true;
		}
		return is_gen_ee_;
	}
	bool &	is_gen_mm()
	{
		if (not is_gen_mm_isLoaded) {
			if (is_gen_mm_branch != 0) {
				is_gen_mm_branch->GetEntry(index);
			} else { 
				printf("branch is_gen_mm_branch does not exist!\n");
				exit(1);
			}
			is_gen_mm_isLoaded = true;
		}
		return is_gen_mm_;
	}
	bool &	is_gen_acc_num()
	{
		if (not is_gen_acc_num_isLoaded) {
			if (is_gen_acc_num_branch != 0) {
				is_gen_acc_num_branch->GetEntry(index);
			} else { 
				printf("branch is_gen_acc_num_branch does not exist!\n");
				exit(1);
			}
			is_gen_acc_num_isLoaded = true;
		}
		return is_gen_acc_num_;
	}
	bool &	is_acc_den()
	{
		if (not is_acc_den_isLoaded) {
			if (is_acc_den_branch != 0) {
				is_acc_den_branch->GetEntry(index);
			} else { 
				printf("branch is_acc_den_branch does not exist!\n");
				exit(1);
			}
			is_acc_den_isLoaded = true;
		}
		return is_acc_den_;
	}
	int &gen_lep1_id()
	{
		if (not gen_lep1_id_isLoaded) {
			if (gen_lep1_id_branch != 0) {
				gen_lep1_id_branch->GetEntry(index);
			} else { 
				printf("branch gen_lep1_id_branch does not exist!\n");
				exit(1);
			}
			gen_lep1_id_isLoaded = true;
		}
		return gen_lep1_id_;
	}
	int &gen_lep1_charge()
	{
		if (not gen_lep1_charge_isLoaded) {
			if (gen_lep1_charge_branch != 0) {
				gen_lep1_charge_branch->GetEntry(index);
			} else { 
				printf("branch gen_lep1_charge_branch does not exist!\n");
				exit(1);
			}
			gen_lep1_charge_isLoaded = true;
		}
		return gen_lep1_charge_;
	}
	int &gen_lep2_id()
	{
		if (not gen_lep2_id_isLoaded) {
			if (gen_lep2_id_branch != 0) {
				gen_lep2_id_branch->GetEntry(index);
			} else { 
				printf("branch gen_lep2_id_branch does not exist!\n");
				exit(1);
			}
			gen_lep2_id_isLoaded = true;
		}
		return gen_lep2_id_;
	}
	int &gen_lep2_charge()
	{
		if (not gen_lep2_charge_isLoaded) {
			if (gen_lep2_charge_branch != 0) {
				gen_lep2_charge_branch->GetEntry(index);
			} else { 
				printf("branch gen_lep2_charge_branch does not exist!\n");
				exit(1);
			}
			gen_lep2_charge_isLoaded = true;
		}
		return gen_lep2_charge_;
	}
	bool &	is_reco_ee()
	{
		if (not is_reco_ee_isLoaded) {
			if (is_reco_ee_branch != 0) {
				is_reco_ee_branch->GetEntry(index);
			} else { 
				printf("branch is_reco_ee_branch does not exist!\n");
				exit(1);
			}
			is_reco_ee_isLoaded = true;
		}
		return is_reco_ee_;
	}
	bool &	is_reco_mm()
	{
		if (not is_reco_mm_isLoaded) {
			if (is_reco_mm_branch != 0) {
				is_reco_mm_branch->GetEntry(index);
			} else { 
				printf("branch is_reco_mm_branch does not exist!\n");
				exit(1);
			}
			is_reco_mm_isLoaded = true;
		}
		return is_reco_mm_;
	}
	bool &	is_reco_acc_num()
	{
		if (not is_reco_acc_num_isLoaded) {
			if (is_reco_acc_num_branch != 0) {
				is_reco_acc_num_branch->GetEntry(index);
			} else { 
				printf("branch is_reco_acc_num_branch does not exist!\n");
				exit(1);
			}
			is_reco_acc_num_isLoaded = true;
		}
		return is_reco_acc_num_;
	}
	int &reco_lep1_id()
	{
		if (not reco_lep1_id_isLoaded) {
			if (reco_lep1_id_branch != 0) {
				reco_lep1_id_branch->GetEntry(index);
			} else { 
				printf("branch reco_lep1_id_branch does not exist!\n");
				exit(1);
			}
			reco_lep1_id_isLoaded = true;
		}
		return reco_lep1_id_;
	}
	int &reco_lep1_charge()
	{
		if (not reco_lep1_charge_isLoaded) {
			if (reco_lep1_charge_branch != 0) {
				reco_lep1_charge_branch->GetEntry(index);
			} else { 
				printf("branch reco_lep1_charge_branch does not exist!\n");
				exit(1);
			}
			reco_lep1_charge_isLoaded = true;
		}
		return reco_lep1_charge_;
	}
	int &reco_lep2_id()
	{
		if (not reco_lep2_id_isLoaded) {
			if (reco_lep2_id_branch != 0) {
				reco_lep2_id_branch->GetEntry(index);
			} else { 
				printf("branch reco_lep2_id_branch does not exist!\n");
				exit(1);
			}
			reco_lep2_id_isLoaded = true;
		}
		return reco_lep2_id_;
	}
	int &reco_lep2_charge()
	{
		if (not reco_lep2_charge_isLoaded) {
			if (reco_lep2_charge_branch != 0) {
				reco_lep2_charge_branch->GetEntry(index);
			} else { 
				printf("branch reco_lep2_charge_branch does not exist!\n");
				exit(1);
			}
			reco_lep2_charge_isLoaded = true;
		}
		return reco_lep2_charge_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_p4()
	{
		if (not gen_p4_isLoaded) {
			if (gen_p4_branch != 0) {
				gen_p4_branch->GetEntry(index);
			} else { 
				printf("branch gen_p4_branch does not exist!\n");
				exit(1);
			}
			gen_p4_isLoaded = true;
		}
		return *gen_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_lep1_p4()
	{
		if (not gen_lep1_p4_isLoaded) {
			if (gen_lep1_p4_branch != 0) {
				gen_lep1_p4_branch->GetEntry(index);
			} else { 
				printf("branch gen_lep1_p4_branch does not exist!\n");
				exit(1);
			}
			gen_lep1_p4_isLoaded = true;
		}
		return *gen_lep1_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_lep2_p4()
	{
		if (not gen_lep2_p4_isLoaded) {
			if (gen_lep2_p4_branch != 0) {
				gen_lep2_p4_branch->GetEntry(index);
			} else { 
				printf("branch gen_lep2_p4_branch does not exist!\n");
				exit(1);
			}
			gen_lep2_p4_isLoaded = true;
		}
		return *gen_lep2_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_p4()
	{
		if (not reco_p4_isLoaded) {
			if (reco_p4_branch != 0) {
				reco_p4_branch->GetEntry(index);
			} else { 
				printf("branch reco_p4_branch does not exist!\n");
				exit(1);
			}
			reco_p4_isLoaded = true;
		}
		return *reco_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_lep1_p4()
	{
		if (not reco_lep1_p4_isLoaded) {
			if (reco_lep1_p4_branch != 0) {
				reco_lep1_p4_branch->GetEntry(index);
			} else { 
				printf("branch reco_lep1_p4_branch does not exist!\n");
				exit(1);
			}
			reco_lep1_p4_isLoaded = true;
		}
		return *reco_lep1_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_lep2_p4()
	{
		if (not reco_lep2_p4_isLoaded) {
			if (reco_lep2_p4_branch != 0) {
				reco_lep2_p4_branch->GetEntry(index);
			} else { 
				printf("branch reco_lep2_p4_branch does not exist!\n");
				exit(1);
			}
			reco_lep2_p4_isLoaded = true;
		}
		return *reco_lep2_p4_;
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern DrellYan dy_obj;
#endif

namespace dy {
	const int &run();
	const int &ls();
	const int &evt();
	const string &sample();
	const string &dataset();
	const string &filename();
	const bool &is_real_data();
	const double &scale1fb();
	const double &scale1fb_cms2();
	const double &lumi();
	const double &xsec();
	const double &filt_eff();
	const bool &is_gen_z();
	const bool &is_gen_tt();
	const bool &is_gen_ee();
	const bool &is_gen_mm();
	const bool &is_gen_acc_num();
	const bool &is_acc_den();
	const int &gen_lep1_id();
	const int &gen_lep1_charge();
	const int &gen_lep2_id();
	const int &gen_lep2_charge();
	const bool &is_reco_ee();
	const bool &is_reco_mm();
	const bool &is_reco_acc_num();
	const int &reco_lep1_id();
	const int &reco_lep1_charge();
	const int &reco_lep2_id();
	const int &reco_lep2_charge();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_lep1_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_lep2_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_lep1_p4();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_lep2_p4();
}
#endif
