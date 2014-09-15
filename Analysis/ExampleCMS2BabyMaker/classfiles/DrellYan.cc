#include "DrellYan.h"
DrellYan dy_obj;
namespace dy {
	const int &run() { return dy_obj.run(); }
	const int &ls() { return dy_obj.ls(); }
	const int &evt() { return dy_obj.evt(); }
	const string &sample() { return dy_obj.sample(); }
	const string &dataset() { return dy_obj.dataset(); }
	const string &filename() { return dy_obj.filename(); }
	const bool &is_real_data() { return dy_obj.is_real_data(); }
	const double &scale1fb() { return dy_obj.scale1fb(); }
	const double &scale1fb_cms2() { return dy_obj.scale1fb_cms2(); }
	const double &lumi() { return dy_obj.lumi(); }
	const double &xsec() { return dy_obj.xsec(); }
	const double &filt_eff() { return dy_obj.filt_eff(); }
	const bool &is_gen_z() { return dy_obj.is_gen_z(); }
	const bool &is_gen_tt() { return dy_obj.is_gen_tt(); }
	const bool &is_gen_ee() { return dy_obj.is_gen_ee(); }
	const bool &is_gen_mm() { return dy_obj.is_gen_mm(); }
	const bool &is_gen_acc_num() { return dy_obj.is_gen_acc_num(); }
	const bool &is_acc_den() { return dy_obj.is_acc_den(); }
	const int &gen_lep1_id() { return dy_obj.gen_lep1_id(); }
	const int &gen_lep1_charge() { return dy_obj.gen_lep1_charge(); }
	const int &gen_lep2_id() { return dy_obj.gen_lep2_id(); }
	const int &gen_lep2_charge() { return dy_obj.gen_lep2_charge(); }
	const bool &is_reco_ee() { return dy_obj.is_reco_ee(); }
	const bool &is_reco_mm() { return dy_obj.is_reco_mm(); }
	const bool &is_reco_acc_num() { return dy_obj.is_reco_acc_num(); }
	const int &reco_lep1_id() { return dy_obj.reco_lep1_id(); }
	const int &reco_lep1_charge() { return dy_obj.reco_lep1_charge(); }
	const int &reco_lep2_id() { return dy_obj.reco_lep2_id(); }
	const int &reco_lep2_charge() { return dy_obj.reco_lep2_charge(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_p4() { return dy_obj.gen_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_lep1_p4() { return dy_obj.gen_lep1_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &gen_lep2_p4() { return dy_obj.gen_lep2_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_p4() { return dy_obj.reco_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_lep1_p4() { return dy_obj.reco_lep1_p4(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > &reco_lep2_p4() { return dy_obj.reco_lep2_p4(); }
}
