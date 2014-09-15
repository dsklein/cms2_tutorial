#include "Analysis/ExampleCMS2BabyMaker/interface/CMS2BabyMaker.h"

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
#include "Math/VectorUtil.h"

// CMS2
#include "CMS2/NtupleMacrosHeader/interface/CMS2.h"
#include "CMS2/NtupleMacrosCore/interface/eventSelections.h"
#include "CMS2/NtupleMacrosCore/interface/trackSelections.h"
#include "CMS2/NtupleMacrosCore/interface/mcSelections.h"
#include "CMS2/NtupleMacrosCore/interface/electronSelections.h"
#include "CMS2/NtupleMacrosCore/interface/MITConversionUtilities.h"
#include "CMS2/NtupleMacrosCore/interface/susySelections.h"

// stuff from packages 
#include "Packages/LooperTools/interface/GoodRun.h"
#include "Packages/LooperTools/interface/DorkyEventIdentifier.h"

// -------------------------------------------------//
// Class to hold the ntuple information 
// -------------------------------------------------//

TreeInfo::TreeInfo()
    : run                  ( -999999    ) 
    , ls                   ( -999999    ) 
    , evt                  ( -999999    ) 
    , sample               ( ""         ) 
    , dataset              ( ""         ) 
    , filename             ( ""         ) 
    , is_real_data         ( false      ) 
    , scale1fb             ( 1.0        ) 
    , scale1fb_cms2        ( 1.0        ) 
    , lumi                 ( 1.0        ) 
    , xsec                 ( -999999    ) 
    , kfactor              ( -999999    ) 
    , filt_eff             ( -999999    ) 
	, is_gen_z             ( false      ) 
	, is_gen_tt            ( false      ) 
    , is_gen_ee            ( false      ) 
    , is_gen_mm            ( false      ) 
	, is_gen_acc_num       ( false      )
	, is_acc_den           ( false      )
    , gen_p4               ( 0, 0, 0, 0 ) 
    , gen_lep1_p4          ( 0, 0, 0, 0 ) 
    , gen_lep1_id          ( -999999    ) 
    , gen_lep1_charge      ( -999999    ) 
    , gen_lep2_p4          ( 0, 0, 0, 0 ) 
    , gen_lep2_id          ( -999999    ) 
    , gen_lep2_charge      ( -999999    ) 
    , is_reco_ee           ( false      ) 
    , is_reco_mm           ( false      ) 
	, is_reco_acc_num      (false       )
    , reco_p4              ( 0, 0, 0, 0 ) 
    , reco_lep1_p4         ( 0, 0, 0, 0 ) 
    , reco_lep1_id         ( -999999    ) 
    , reco_lep1_charge     ( -999999    ) 
    , reco_lep2_p4         ( 0, 0, 0, 0 ) 
    , reco_lep2_id         ( -999999    ) 
    , reco_lep2_charge     ( -999999    ) 
{
}

void TreeInfo::Reset()
{
    run                  = -999999; 
    ls                   = -999999; 
    evt                  = -999999; 
    sample               = ""; 
    dataset              = ""; 
    filename             = ""; 
    is_real_data         = false; 
    scale1fb             = 1.0; 
    scale1fb_cms2        = 1.0; 
    lumi                 = 1.0; 
    xsec                 = -999999; 
    kfactor              = -999999; 
    filt_eff             = -999999; 
    is_gen_z             = false; 
    is_gen_tt            = false; 
    is_gen_ee            = false; 
    is_gen_mm            = false; 
	is_gen_acc_num       = false;
	is_acc_den           = false;
    gen_p4               = LorentzVector(0, 0, 0, 0); 
    gen_lep1_p4          = LorentzVector(0, 0, 0, 0); 
    gen_lep1_id          = -999999; 
    gen_lep1_charge      = -999999; 
    gen_lep2_p4          = LorentzVector(0, 0, 0, 0); 
    gen_lep2_id          = -999999; 
    gen_lep2_charge      = -999999; 
    is_reco_ee            = false; 
    is_reco_mm            = false; 
	is_reco_acc_num       = false;
    reco_p4               = LorentzVector(0, 0, 0, 0); 
    reco_lep1_p4          = LorentzVector(0, 0, 0, 0); 
    reco_lep1_id          = -999999; 
    reco_lep1_charge      = -999999; 
    reco_lep2_p4          = LorentzVector(0, 0, 0, 0); 
    reco_lep2_id          = -999999; 
    reco_lep2_charge      = -999999; 
}
    
void TreeInfo::SetBranches(TTree& tree)
{
    tree.Branch("run"                  , &run                  );
    tree.Branch("ls"                   , &ls                   );
    tree.Branch("evt"                  , &evt                  );
    tree.Branch("sample"               , &sample               );
    tree.Branch("dataset"              , &dataset              );
    tree.Branch("filename"             , &filename             );
    tree.Branch("is_real_data"         , &is_real_data         );
    tree.Branch("scale1fb"             , &scale1fb             );
    tree.Branch("scale1fb_cms2"        , &scale1fb_cms2        );
    tree.Branch("lumi"                 , &lumi                 );
    tree.Branch("xsec"                 , &xsec                 );
    tree.Branch("filt_eff"             , &filt_eff             );
    tree.Branch("is_gen_z"             , &is_gen_z             );
    tree.Branch("is_gen_tt"            , &is_gen_tt            );
    tree.Branch("is_gen_ee"            , &is_gen_ee            );
    tree.Branch("is_gen_mm"            , &is_gen_mm            );
	tree.Branch("is_gen_acc_num"       , &is_gen_acc_num       );
	tree.Branch("is_acc_den"           , &is_acc_den           );
    tree.Branch("gen_lep1_id"          , &gen_lep1_id          );
    tree.Branch("gen_lep1_charge"      , &gen_lep1_charge      );
    tree.Branch("gen_lep2_id"          , &gen_lep2_id          );
    tree.Branch("gen_lep2_charge"      , &gen_lep2_charge      );
    tree.Branch("is_reco_ee"           , &is_reco_ee           );
    tree.Branch("is_reco_mm"           , &is_reco_mm           );
	tree.Branch("is_reco_acc_num"      , &is_reco_acc_num      );
    tree.Branch("reco_lep1_id"         , &reco_lep1_id         );
    tree.Branch("reco_lep1_charge"     , &reco_lep1_charge     );
    tree.Branch("reco_lep2_id"         , &reco_lep2_id         );
    tree.Branch("reco_lep2_charge"     , &reco_lep2_charge     );

    tree.Branch("gen_p4"       , "LorentzVector" , &gen_p4     );
    tree.Branch("gen_lep1_p4"  , "LorentzVector" , &gen_lep1_p4);
    tree.Branch("gen_lep2_p4"  , "LorentzVector" , &gen_lep2_p4);
    tree.Branch("reco_p4"      , "LorentzVector" , &reco_p4     );
    tree.Branch("reco_lep1_p4" , "LorentzVector" , &reco_lep1_p4);
    tree.Branch("reco_lep2_p4" , "LorentzVector" , &reco_lep2_p4);
}

std::ostream& operator<< (std::ostream& out, const TreeInfo& info)
{
    out << "run                  = " << info.run                  << std::endl;
    out << "ls                   = " << info.ls                   << std::endl;
    out << "evt                  = " << info.evt                  << std::endl;
    out << "sample               = " << info.sample               << std::endl;
    out << "dataset              = " << info.dataset              << std::endl;
    out << "filename             = " << info.filename             << std::endl;
    out << "is_real_data         = " << info.is_real_data         << std::endl;
    out << "scale1fb             = " << info.scale1fb             << std::endl;
    out << "scale1fb_cms2        = " << info.scale1fb_cms2        << std::endl;
    out << "lumi                 = " << info.lumi                 << std::endl;
    out << "xsec                 = " << info.xsec                 << std::endl;
    out << "kfactor              = " << info.kfactor              << std::endl;
    out << "filt_eff             = " << info.filt_eff             << std::endl;
    out << "is_gen_z             = " << info.is_gen_z             << std::endl;
    out << "is_gen_tt            = " << info.is_gen_tt            << std::endl;
	out << "is_gen_ee            = " << info.is_gen_ee            << std::endl;
    out << "is_gen_mm            = " << info.is_gen_mm            << std::endl;
	out << "is_gen_acc_num       = " << info.is_gen_acc_num       << std::endl;
	out << "is_acc_den           = " << info.is_acc_den           << std::endl;
    out << "gen_p4.mass()        = " << info.gen_p4.mass()        << std::endl;
    out << "gen_lep1_p4.pt()     = " << info.gen_lep1_p4.pt()     << std::endl;
    out << "gen_lep1_id          = " << info.gen_lep1_id          << std::endl;
    out << "gen_lep1_charge      = " << info.gen_lep1_charge      << std::endl;
    out << "gen_lep2_p4.pt()     = " << info.gen_lep2_p4.pt()     << std::endl;
    out << "gen_lep2_id          = " << info.gen_lep2_id          << std::endl;
    out << "gen_lep2_charge      = " << info.gen_lep2_charge      << std::endl;
    out << "is_reco_ee           = " << info.is_reco_ee           << std::endl;
    out << "is_reco_mm           = " << info.is_reco_mm           << std::endl;
	out << "is_reco_acc_num      = " << info.is_reco_acc_num      << std::endl;
    out << "reco_p4.mass()       = " << info.reco_p4.mass()       << std::endl;
    out << "reco_lep1_p4.pt()    = " << info.reco_lep1_p4.pt()    << std::endl;
    out << "reco_lep1_id         = " << info.reco_lep1_id         << std::endl;
    out << "reco_lep1_charge     = " << info.reco_lep1_charge     << std::endl;
    out << "reco_lep2_p4.pt()    = " << info.reco_lep2_p4.pt()    << std::endl;
    out << "reco_lep2_id         = " << info.reco_lep2_id         << std::endl;
    out << "reco_lep2_charge     = " << info.reco_lep2_charge     << std::endl;
    return out;
}

// ------------------------------------ //
//  Looper class 
// ------------------------------------ //

// construct:
CMS2BabyMaker::CMS2BabyMaker
(
    const std::string& sample_name, 
    const std::string& output_filename,
    const std::string& runlist_filename,
    const double lumi,
    const double nevts_corr,
    const bool verbose
)
    : m_sample_name(sample_name)
    , m_output_filename(output_filename)
    , m_runlist_filename(runlist_filename)
    , m_nevts_corr(nevts_corr)
    , m_lumi(lumi)
    , m_verbose(verbose)
    , m_info()
    , m_file(*TFile::Open(output_filename.c_str(), "RECREATE"))
    , m_tree(*new TTree("tree", "CMS2 Example Baby TTree"))
{
}

// ------------------------------------ //
// Stuff to do before job starts
// ------------------------------------ //

void CMS2BabyMaker::BeginJob()
{
    // Set all the branches to the TTree
    m_info.SetBranches(m_tree);
}

// ------------------------------------ //
// simple lepton information 
// ------------------------------------ //

struct LepInfo
{
    LorentzVector p4;
    int id;
    int charge;
    int mom_id;
};

// simple dilepton hypothesis 
struct DilHyp
{
    LepInfo lep1;
    LepInfo lep2;
};

// compare the hypothesis
bool CompareDilHyp(const DilHyp& hyp1, const DilHyp& hyp2)
{
    int hyp1_type = -1;
    switch (hyp1.lep1.id*hyp1.lep2.id)
    {
        case -11*11: hyp1_type = 2; break; // ee
        case -13*13: hyp1_type = 1; break; // mu mu
    	case -15*15: hyp1_type = 3; break; // tau tau
    }
    int hyp2_type = -1;
    switch (hyp2.lep1.id*hyp2.lep2.id)
    {
        case -11*11: hyp2_type = 2; break; // ee
        case -13*13: hyp2_type = 1; break; // mu mu
    	case -15*15: hyp2_type = 3; break; // tau tau
    }

    // choose mm over ee over tau tau
    if (hyp1_type != hyp2_type)
    {
        return (hyp1_type < hyp2_type);
    }
    else
    {
        // choose one closer to pdg value of mass of Z-boson
        static const double Mz = 91.1876;
        const double dm1 = fabs((hyp1.lep1.p4 + hyp1.lep2.p4).mass() - Mz);
        const double dm2 = fabs((hyp2.lep1.p4 + hyp2.lep2.p4).mass() - Mz);
        return (dm1 < dm2);
    }
    return true;
}

// ------------------------------------ //
// Helper functions
// ------------------------------------ //

double LeptonD0(const int lep_id, const int lep_idx)
{
    const int vtxidx = firstGoodVertex();
    if (vtxidx < 0)
    {
        std::cout << "[LeptonD0] WARNING - first good vertex index < 0.  Returning bogus value 999999" << std::endl;
        return 999999.0;
    }
    if (abs(lep_id)==13)
    {
        const int trkidx = tas::mus_trkidx().at(lep_idx);
        if (trkidx >= 0)
        {
            return trks_d0_pv(trkidx, vtxidx).first;
        }
    }
    else if (abs(lep_id)==11)
    {
        const int gsfidx = tas::els_gsftrkidx().at(lep_idx);
        if (gsfidx >= 0) 
        {
            return gsftrks_d0_pv(gsfidx, vtxidx).first;
        }
        const int trkidx = tas::els_trkidx().at(lep_idx);
        if (trkidx >= 0)
        {
            return trks_d0_pv(trkidx, vtxidx).first;
        }
    }

    // return bogus for non electon/muon
    return -999999.0;
}

double LeptonDz(const int lep_id, const int lep_idx)
{
    const int vtxidx = firstGoodVertex();
    if (vtxidx < 0)
    {
        std::cout << "[LeptonDz] WARNING - first good vertex index < 0.  Returning bogus value 999999" << std::endl;
        return 999999.0;
    }
    if (abs(lep_id)==13)
    {
        const int trkidx = tas::mus_trkidx().at(lep_idx);
        if (trkidx >= 0)
        {
            return trks_dz_pv(trkidx, vtxidx).first;
        }
    }
    else if (abs(lep_id)==11)
    {
        const int gsfidx = tas::els_gsftrkidx().at(lep_idx);
        if (gsfidx >= 0)
        {
            return gsftrks_dz_pv(gsfidx, vtxidx).first;
        }
        const int trkidx = tas::els_trkidx().at(lep_idx);
        if (trkidx >= 0)
        {
            return trks_dz_pv(trkidx, vtxidx).first;
        }
    }

    // return bogus for non electon/muon
    return -999999.0;
}

double electronIsolationPF2012_cone03(const int el_idx)
{
    // electron pT
    const double pt = tas::els_p4().at(el_idx).pt();

    // get effective area
    const double AEff = fastJetEffArea03_v2( fabs(tas::els_etaSC().at(el_idx)) );

    // pf iso
    const double pfiso_ch = tas::els_iso03_pf2012ext_ch().at(el_idx);
    const double pfiso_em = tas::els_iso03_pf2012ext_em().at(el_idx);
    const double pfiso_nh = tas::els_iso03_pf2012ext_nh().at(el_idx);

    // rho
    const double rhoPrime = std::max(tas::evt_kt6pf_foregiso_rho(), 0.0f);
    const double pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, 0.0);
    const double pfiso = (pfiso_ch + pfiso_n) / pt;

    return pfiso;
}

double muonIsoValuePF2012(const int mu_idx)
{
    const double chiso = tas::mus_isoR04_pf_ChargedHadronPt().at(mu_idx);
    const double nhiso = tas::mus_isoR04_pf_NeutralHadronEt().at(mu_idx);
    const double emiso = tas::mus_isoR04_pf_PhotonEt().at(mu_idx);
    const double deltaBeta = tas::mus_isoR04_pf_PUPt().at(mu_idx);
    const double pt = tas::mus_p4().at(mu_idx).pt();
    const double absiso = chiso + max(0.0, nhiso + emiso - 0.5 * deltaBeta);
    return (absiso / pt);
}

bool lepIsKinematic(const LorentzVector p4, const int lepid)
{
  double abseta = fabs(p4.eta());

  if( abs(lepid)==13 ) {
	if( p4.pt()<25 ) return false;
	if( abseta > 2.1 ) return false;
	return true;
  }
  else if( abs(lepid)==11 ) {
	if( p4.Et()<25 ) return false;
	if( abseta>2.5 ) return false;
	if( abseta>1.4442 && abseta<1.566 ) return false;
	return true;
  }
  else if( abs(lepid)==15 ) {
	if( p4.pt()<25 ) return false;
	if( abseta > 2.5 ) return false;
	return true;
  }
  else return false;
}

// ------------------------------------ //
// Stuff to do on each event 
// ------------------------------------ //

void CMS2BabyMaker::Analyze(const long event, const std::string& current_filename)
{
    if (m_verbose)
    {
        std::cout << "\n[DrellYanNtupleMaker] Running on run, ls, event: " 
            << tas::evt_run()       << ", "
            << tas::evt_lumiBlock() << ", "
            << tas::evt_event()     << std::endl;
    }

    // reset the TTree variables
    // ---------------------- // 
    m_info.Reset();

    // event information 
    // ---------------------- // 

    m_info.run          = tas::evt_run();
    m_info.ls           = tas::evt_lumiBlock();
    m_info.evt          = tas::evt_event();
    m_info.sample       = m_sample_name; 
    m_info.dataset      = tas::evt_dataset().front().Data();
    m_info.filename     = current_filename; 
    m_info.is_real_data = tas::evt_isRealData(); 
    if (!tas::evt_isRealData())
    {
        m_info.scale1fb_cms2 = tas::evt_scale1fb();
        m_info.scale1fb      = m_info.scale1fb_cms2 * m_nevts_corr; 
        m_info.lumi          = m_lumi; 
        m_info.xsec          = tas::evt_xsec_excl();
        m_info.kfactor       = tas::evt_kfactor();
        m_info.filt_eff      = tas::evt_filt_eff();
    }

    // gen information 
    // ---------------------- // 

    if (!tas::evt_isRealData())
    {
        // loop over gen particles and create a list of status 3 leptons 
        std::vector<LepInfo> gen_infos;
        for (size_t gen_idx = 0; gen_idx < tas::genps_id().size(); ++gen_idx)
        {
            // only keep status 3
            if (cms2.genps_status().at(gen_idx) != 3) {continue;}
        
            // only keep charged leptons
            const int id            = tas::genps_id().at(gen_idx);
			const int charge        = -1*id/abs(id); // e.g. 11 == e^- and -11 == e+
			const LorentzVector& p4 = tas::genps_p4().at(gen_idx);
			const int mom_id        = tas::genps_id_mother().at(gen_idx);
			const LepInfo gen_info{p4, id, charge, mom_id}; 

            if (abs(id) == 11 || abs(id) == 13)
            {
                gen_infos.push_back(gen_info);
            }
            if (abs(id) == 15)
            {
                // Only consider tau --> nu_tau + nu_lep + lep events.
                // We count neutrinos because that guarantees that 
                // there is a corresponding lepton and that it comes from
                // a leptonic tau decay. You can get electrons from converted photons
                // which are radiated by charged pions from the tau decay, but that's
                // hadronic and we don't care for those.
                int nu_count = 0;
                for (const int& d_id : tas::genps_lepdaughter_id().at(gen_idx))
                {
                    if (abs(d_id)==12 || abs(d_id)==14) ++nu_count;
                }
                if (nu_count < 1) {continue;}

                // now find the lepton
                for (size_t d_idx = 0; d_idx != tas::genps_lepdaughter_id().at(gen_idx).size(); ++d_idx)
                {
                    const int d_id = tas::genps_lepdaughter_id().at(gen_idx).at(d_idx);
                    if (abs(d_id)==11 || abs(d_id)==13)
                    {
                        // const int d_charge        = -1*d_id/abs(d_id); // e.g. 11 == e^- and -11 == e+
                        // const LorentzVector& d_p4 = tas::genps_lepdaughter_p4().at(gen_idx).at(d_idx);
                        // const int d_mom_id        = id; 
                        // const LepInfo daughter_info{d_p4, d_id, d_charge, d_mom_id}; 
                        // gen_infos.push_back(daughter_info);
						gen_infos.push_back(gen_info);
                        break;
                    }
                }
            }
        }

        // form the gen hypothesis pairs from the LepInfos 
        std::vector<DilHyp> gen_hyps;
        for (size_t idx1 = 0; idx1 < gen_infos.size(); idx1++)
        {
            for (size_t idx2 = idx1 + 1; idx2 < gen_infos.size(); idx2++)
            {
                // sort by pt
                const LepInfo& gen_l1 = (gen_infos.at(idx1).p4.pt() > gen_infos.at(idx2).p4.pt() ? gen_infos.at(idx1) : gen_infos.at(idx2));
                const LepInfo& gen_l2 = (gen_infos.at(idx1).p4.pt() > gen_infos.at(idx2).p4.pt() ? gen_infos.at(idx2) : gen_infos.at(idx1));
                const DilHyp gen_hyp{gen_l1, gen_l2};

                // ensure opposite sign
                if (gen_hyp.lep1.charge == gen_hyp.lep2.charge) {continue;}

                // ensure same flavor
                if (abs(gen_hyp.lep1.id) != abs(gen_hyp.lep2.id)) {continue;}

                gen_hyps.push_back(gen_hyp);
            }
        }

        // sort the gen hyps in order of precedence
        std::sort(gen_hyps.begin(), gen_hyps.end(), CompareDilHyp);

        // fill the tree info with gen info
        if (m_verbose) {std::cout << "number of gen hyps = " << gen_hyps.size() << std::endl;}
        if (!gen_hyps.empty())
        {
            const DilHyp& gen_hyp = gen_hyps.front();

			//dytt baby should exclude all gen ee/mm events
			if( m_sample_name=="dytt" && (abs(gen_hyp.lep1.id)==11 || abs(gen_hyp.lep1.id)==13)  ) return;

            // assign info 
            m_info.is_gen_z        = (gen_hyp.lep1.mom_id==23 && gen_hyp.lep2.mom_id==23);
			m_info.is_gen_tt       = (gen_hyp.lep1.id*gen_hyp.lep2.id == -15*15);
            m_info.is_gen_ee       = (gen_hyp.lep1.id*gen_hyp.lep2.id == -11*11);
            m_info.is_gen_mm       = (gen_hyp.lep1.id*gen_hyp.lep2.id == -13*13);
            m_info.gen_p4          = (gen_hyp.lep1.p4 + gen_hyp.lep2.p4);
            m_info.gen_lep1_p4     = gen_hyp.lep1.p4;
            m_info.gen_lep1_id     = gen_hyp.lep1.id;
            m_info.gen_lep1_charge = gen_hyp.lep1.charge;
            m_info.gen_lep2_p4     = gen_hyp.lep2.p4;
            m_info.gen_lep2_id     = gen_hyp.lep2.id;
            m_info.gen_lep2_charge = gen_hyp.lep2.charge;

			m_info.is_acc_den      = (m_info.is_gen_z && !m_info.is_gen_tt && m_info.gen_p4.mass() > 60 && m_info.gen_p4.mass() < 120 );
			m_info.is_gen_acc_num  = (m_info.is_acc_den && lepIsKinematic(gen_hyp.lep1.p4, 15) && lepIsKinematic(gen_hyp.lep2.p4, 15));
        }
    }


    // reco information 
    // ---------------------- //

	bool lepdebug = false;
	//if( m_info.run==190782 && m_info.ls==195 && m_info.evt==229129850 ) lepdebug = true;
	//if( m_info.run==190895 && m_info.ls==827 && m_info.evt==826363592 ) lepdebug = true;
	if( lepdebug ) std::cout << event << "\t" << m_info.run << "\t" << m_info.ls << "\t" << m_info.evt << std::endl;

	std::vector<LepInfo> reco_infos;

	// Perform the electron selections
    // ----------------------------------------------------- // 

	bool passes_ee_trigger = true;
	bool passes_mm_trigger = true;

	// Check which triggers the event passes
	if( tas::evt_isRealData() ) {
	  if( !passUnprescaledHLTTriggerPattern("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") ) {
		if( lepdebug ) std::cout << "This event doesn't pass the electron trigger." << std::endl;
		passes_ee_trigger = false;
	  }
	  if( !passUnprescaledHLTTriggerPattern("HLT_Mu17_Mu8_v") &&
		  !passUnprescaledHLTTriggerPattern("HLT_Mu17_TkMu8_v") ) {
		if( lepdebug ) std::cout << "This event doesn't pass the muon triggers." << std::endl;
		passes_mm_trigger = false;
	  }
	}

	// Reco-level ID and iso cuts
	for( size_t el_idx = 0; el_idx<tas::els_p4().size(); ++el_idx ) {

	  if( lepdebug ) {
		std::cout << "\nElectron " << el_idx+1 << " of " << tas::els_p4().size() << ":" << std::endl;
	  }

	  const LorentzVector& el_p4 = tas::els_p4().at(el_idx);

	  bool is_barrel = false;
	  bool is_endcap = false;

	  //Transverse energy cut
	  if( lepdebug ) std::cout << "ET = " << el_p4.Et() << std::endl;
	  if( el_p4.Et() < 25. ) continue;

	  //Eta cut and classification
	  const double el_etaSC      = tas::els_etaSC().at(el_idx);
	  if( lepdebug ) std::cout << "etaSC = " << el_etaSC << std::endl;
	  if( fabs(el_etaSC) > 2.5 ) continue;
	  else if( fabs(el_etaSC) < 1.4442)  is_barrel = true; //barrel electron
	  else if( fabs(el_etaSC) > 1.566 )  is_endcap = true; //endcap electron
	  else continue; //electron is in the crack

	  //PF Isolation cut
	  if( lepdebug ) std::cout << "iso = " << electronIsolationPF2012_cone03(el_idx) << std::endl;
	  if( electronIsolationPF2012_cone03(el_idx) > 0.15 ) continue;

	  //Sigma_ietaieta cut
	  const double el_sieie = tas::els_sigmaIEtaIEta().at(el_idx);
	  if( lepdebug ) std::cout << "sieie = " << el_sieie << std::endl;
	  if( is_barrel && (el_sieie > 0.01) ) continue;
	  else if( is_endcap && (el_sieie > 0.03) ) continue;

	  //Delta phi_in cut
	  const double el_dphi = tas::els_dPhiIn().at(el_idx);
	  if( lepdebug ) std::cout << "dphi = " << el_dphi << std::endl;
	  if( is_barrel && (fabs(el_dphi)>0.06) ) continue;
	  else if( is_endcap && (fabs(el_dphi)>0.03) ) continue;

	  //Delta eta_in cut
	  const double el_deta = tas::els_dEtaIn().at(el_idx);
	  if( lepdebug ) std::cout << "deta = " << el_deta << std::endl;
	  if( is_barrel && (fabs(el_deta)>0.004) ) continue;
	  else if( is_endcap && (fabs(el_deta)>0.007) ) continue;

	  // |1/E - 1/p| cut
	  if( lepdebug ) std::cout << "1/E - 1/p = " << fabs((1. - tas::els_eOverPIn().at(el_idx)) / tas::els_ecalEnergy().at(el_idx)) << std::endl;
	  if( fabs((1. - tas::els_eOverPIn().at(el_idx)) / tas::els_ecalEnergy().at(el_idx)) > 0.05 ) continue; 

	  // H over E cut
	  const double HoverE = tas::els_hOverE().at(el_idx);
	  if( lepdebug ) std::cout << "HoverE = " << HoverE << std::endl;
	  if( is_barrel && (HoverE > 0.12) ) continue;
	  else if( is_endcap && (HoverE > 0.10) ) continue;

	  // (charge and pdgid bookkeeping)
	  const int el_charge = tas::els_charge().at(el_idx);
	  const int el_pdgid = el_charge * -11;

	  //D0,PV cut
	  if( lepdebug ) std::cout << "d0 = " << LeptonD0(el_pdgid, el_idx) << std::endl;
	  if( fabs(LeptonD0(el_pdgid, el_idx)) > 0.02 ) continue;

	  //DZ,PV cut
	  if( lepdebug ) std::cout << "dZ = " << LeptonDz(el_pdgid, el_idx) << std::endl;
	  if( fabs(LeptonDz(el_pdgid, el_idx)) > 0.1 ) continue;

	  //Missing hits cut
	  if( lepdebug ) std::cout << "missing hits = " << tas::els_exp_innerlayers().at(el_idx) << std::endl;
	  if( tas::els_exp_innerlayers().at(el_idx) > 1 ) continue;

	  //Conversion fit cut
	  if( lepdebug ) std::cout << "conversion fit = " << isMITConversion(el_idx, 0, 1e-6, 2.0, true, false) << std::endl;
	  if( isMITConversion(el_idx, 0, 1e-6, 2.0, true, false) ) continue;

	  //////////////////////////////////////////////////////
	  /////// If we've gotten this far, the electron passes! 

	  if( lepdebug ) std::cout << "This electron passes." << std::endl;

	  const LepInfo el_info{ el_p4, el_pdgid, el_charge, int(el_idx) };
	  reco_infos.push_back( el_info );

	}

	// Perform the muon selections
    // ----------------------------------------------------- //

	for( size_t mu_idx = 0; mu_idx<tas::mus_p4().size(); ++mu_idx ) {

	  if( lepdebug ) std::cout << "\nMuon " << mu_idx+1 << " of " << tas::mus_p4().size() << ":" << std::endl;

	  const LorentzVector& mu_p4 = tas::mus_p4().at(mu_idx);

	  //pT cut
	  if( lepdebug ) std::cout << "pT = " << mu_p4.pt() << std::endl;
	  if( mu_p4.pt() < 25. ) continue;

	  //eta (fiducial) cut
	  if( lepdebug ) std::cout << "eta = " << mu_p4.eta() << std::endl;
	  if( fabs(mu_p4.eta()) > 2.1 ) continue;

	  //PF isolation cut
	  if( lepdebug ) std::cout << "iso = " << muonIsoValuePF2012(mu_idx) << std::endl;
	  if( muonIsoValuePF2012(mu_idx) > 0.12 ) continue;

	  //Global muon ID requirement
	  if( lepdebug ) std::cout << "globalmuonID = " << (tas::mus_type().at(mu_idx) & (1<<1)) << std::endl;
	  if( (tas::mus_type().at(mu_idx) & (1<<1)) == 0 ) continue;

	  //PF muon ID requirement
	  if( lepdebug ) std::cout << "PFmuonID = " << (tas::mus_type().at(mu_idx) & (1<<5)) << std::endl;
	  if( (tas::mus_type().at(mu_idx) & (1<<5)) == 0 ) continue;

	  //chi2/ndof cut
	  if( lepdebug ) std::cout << "chi2/ndof = " << (tas::mus_gfit_chi2().at(mu_idx) / tas::mus_gfit_ndof().at(mu_idx)) << std::endl;
	  if( (tas::mus_gfit_chi2().at(mu_idx) / tas::mus_gfit_ndof().at(mu_idx)) > 10. ) continue;

	  //Valid standalone muon hits cut
	  if( lepdebug ) std::cout << "vstahits = " << tas::mus_gfit_validSTAHits().at(mu_idx) << std::endl;
	  if( tas::mus_gfit_validSTAHits().at(mu_idx) < 1 ) continue;

	  //Matched stations cut
	  if( lepdebug ) std::cout << "matched stations = " << tas::mus_numberOfMatchedStations().at(mu_idx) << std::endl;
	  if( tas::mus_numberOfMatchedStations().at(mu_idx) < 2 ) continue;

	  // (some track bookkeeping)
	  const int mu_trk_idx = tas::mus_trkidx().at(mu_idx);
	  if( mu_trk_idx < 0 ) continue;

	  //Tracker layers cut
	  if( lepdebug ) std::cout << "tracker layers = " << tas::trks_nlayers().at(mu_trk_idx) << std::endl;
	  if( tas::trks_nlayers().at(mu_trk_idx) < 6 ) continue;

	  //Valid pixel hits cut
	  if( lepdebug ) std::cout << "valid pixel hits = " << tas::trks_valid_pixelhits().at(mu_trk_idx) << std::endl;
	  if( tas::trks_valid_pixelhits().at(mu_trk_idx) < 1 ) continue;

	  // (charge and pdgid bookkeeping)
	  const int mu_charge = tas::mus_charge().at(mu_idx);
	  const int mu_pdgid  = mu_charge * -13;

	  //D0,PV cut
	  if( lepdebug ) std::cout << "d0 = " << LeptonD0(mu_pdgid, mu_idx) << std::endl;
	  if( fabs(LeptonD0(mu_pdgid, mu_idx)) > 0.02 ) continue;

	  //DZ,PV cut
	  if( lepdebug ) std::cout << "dZ = " << LeptonDz(mu_pdgid, mu_idx) << std::endl;
	  if( fabs(LeptonDz(mu_pdgid, mu_idx)) > 0.5 ) continue;


	  //////////////////////////////////////////////////////
	  /////// If we've gotten this far, the muon passes! 

	  if( lepdebug ) std::cout << "This muon passes." << std::endl;

	  const LepInfo mu_info{ mu_p4, mu_pdgid, mu_charge, int(mu_idx) };
	  reco_infos.push_back( mu_info );
	}

	// Match the leptons into hyps, and find the best hyp
	// -------------------------------------------------- //

	std::vector<DilHyp> reco_hyps;
	for (size_t idx1 = 0; idx1 < reco_infos.size(); idx1++) {
	  for (size_t idx2 = idx1 + 1; idx2 < reco_infos.size(); idx2++) {
		// sort by pt
		const LepInfo& reco_l1 = (reco_infos.at(idx1).p4.pt() > reco_infos.at(idx2).p4.pt() ? reco_infos.at(idx1) : reco_infos.at(idx2));
		const LepInfo& reco_l2 = (reco_infos.at(idx1).p4.pt() > reco_infos.at(idx2).p4.pt() ? reco_infos.at(idx2) : reco_infos.at(idx1));
		const DilHyp reco_hyp{reco_l1, reco_l2};

		// ensure opposite sign
		if (reco_hyp.lep1.charge == reco_hyp.lep2.charge) {continue;}

		// ensure same flavor
		if (abs(reco_hyp.lep1.id) != abs(reco_hyp.lep2.id)) {continue;}

		// ensure inside the mass window
		if( (reco_hyp.lep1.p4 + reco_hyp.lep2.p4).M() < 60 ||
			(reco_hyp.lep1.p4 + reco_hyp.lep2.p4).M() > 120 ) continue;

		// ensure appropriate trigger is passed
		if( abs(reco_hyp.lep1.id)==11 && !passes_ee_trigger ) continue;
		if( abs(reco_hyp.lep1.id)==13 && !passes_mm_trigger ) continue;

		reco_hyps.push_back(reco_hyp);

		// If we've gotten this far, this event can potentially go into the reco acceptance numerator
		if( !m_info.is_acc_den ) continue;

		m_info.is_reco_acc_num = true;

	  }
	}

	std::sort(reco_hyps.begin(), reco_hyps.end(), CompareDilHyp);


    // Fill the treeinfo member of the CMS2BabyMaker class
    // --------------------------------------------------- // 

	if (m_verbose) {std::cout << "number of reco hyps = " << reco_hyps.size() << std::endl;}
	if (!reco_hyps.empty())
	  {
		const DilHyp& reco_hyp = reco_hyps.front();

		// assign info 
		m_info.is_reco_ee       = (reco_hyp.lep1.id*reco_hyp.lep2.id == -11*11);
		m_info.is_reco_mm       = (reco_hyp.lep1.id*reco_hyp.lep2.id == -13*13);
		m_info.reco_p4          = (reco_hyp.lep1.p4 + reco_hyp.lep2.p4);
		m_info.reco_lep1_p4     = reco_hyp.lep1.p4;
		m_info.reco_lep1_id     = reco_hyp.lep1.id;
		m_info.reco_lep1_charge = reco_hyp.lep1.charge;
		m_info.reco_lep2_p4     = reco_hyp.lep2.p4;
		m_info.reco_lep2_id     = reco_hyp.lep2.id;
		m_info.reco_lep2_charge = reco_hyp.lep2.charge;

		if( lepdebug ) {
		  if( m_info.is_reco_ee ) std::cout << "We picked electrons " << reco_hyp.lep1.mom_id+1 << " and " << reco_hyp.lep2.mom_id+1 << std::endl;
		  if( m_info.is_reco_mm ) std::cout << "We picked muons "     << reco_hyp.lep1.mom_id+1 << " and " << reco_hyp.lep2.mom_id+1 << std::endl;
		}
	  }

    // fill the tree
    // ---------------------- // 
    if (m_verbose) {std::cout << m_info << std::endl;}

    m_tree.Fill();

    // done
    return;
}

// ------------------------------------ //
// Stuff to do after job finishes
// ------------------------------------ //

void CMS2BabyMaker::EndJob()
{
    std::cout << "[CMS2BabyMaker] Saving hists to output file: " << m_output_filename << std::endl;
    m_file.cd(); 
    m_tree.Write();
    m_file.Close(); 
    return;
}

// ------------------------------------ //
// loop over event 
// ------------------------------------ //

void CMS2BabyMaker::ScanChain(TChain& chain, const long num_events)
{
    // Benchmark
    TBenchmark bmark;
    bmark.Start("benchmark");

    // Stuff to do before job starts
	BeginJob();

    //~-~-~-~-~-~-~-~-~-~-~-//
    //Set json here for data//
    //~-~-~-~-~-~-~-~-~-~-~-//

    // set the "good run" list 
    if (!m_runlist_filename.empty())
    {
        set_goodrun_file(m_runlist_filename.c_str());
    }
	if( m_sample_name == "data" ) std::cout << "Using good run list: " << m_runlist_filename << std::endl;

    //~-~-~-~-~-~-~-~-~-~-~-~-~-~-//
    // Loop over events to Analyze//
    //~-~-~-~-~-~-~-~-~-~-~-~-~-~-//

    // number of events to run on
    const unsigned long num_events_chain = (num_events < 0 or num_events > chain.GetEntries() ? chain.GetEntries() : num_events);

    // count the total and duplicates events
    int i_permille_old             = 0;
    unsigned long num_events_total = 0;
    unsigned long num_duplicates   = 0;
    unsigned long num_bad_events   = 0;

    // files list --> TChain doesn't propoerly with CMS2.h so 
    // we have to looper manually over the files :(
    TObjArray* const list_of_files = chain.GetListOfFiles();
    TIter file_iter(list_of_files);
    TFile* current_file = NULL;

    // loop over files in the chain
    while ((current_file = (TFile*)file_iter.Next()))
    {
        // get the trees in each file
        TFile* const file = TFile::Open(current_file->GetTitle());
        TTree* const tree = dynamic_cast<TTree*>(file->Get("Events"));

        // initialize the chain
        cms2.Init(tree);
        TTreeCache::SetLearnEntries(10);
        chain.SetCacheSize(128*1024*1024);

        // Loop over Events in current file
        if (num_events_total >= num_events_chain) {break;}
        const long num_events_tree = tree->GetEntriesFast();

        // loop over events
        for (long event = 0; event != num_events_tree; ++event)
        {
            // quit if the total is > the number in the chain
            if (num_events_total >= num_events_chain) {continue;}

            // Get Event Content
            tree->LoadTree(event);
            cms2.GetEntry(event);
            ++num_events_total;

            //parse events from json
            if (tas::evt_isRealData())
            {
                // for read data, check that the run is good
                if (!m_runlist_filename.empty())
                {
                    if (!goodrun_json(tas::evt_run(), tas::evt_lumiBlock()))
                    {
                        num_bad_events++;
                        continue;
                    }
                }

                // check for duplicates
                const DorkyEventIdentifier id = {tas::evt_run(), tas::evt_event(), tas::evt_lumiBlock()};
                if (is_duplicate(id))
                {
                    // cout << "\t! ERROR: found duplicate." << endl;
                    num_duplicates++;
                    continue;
                }
            }

            // Progress
            const int i_permille = floor(1000 * num_events_total/ static_cast<float>(num_events_chain));
            if (i_permille != i_permille_old)
			  {
                printf("  \015\033[32m ---> \033[1m\033[31m%4.1f%%" "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                fflush(stdout);
                i_permille_old = i_permille;
			  }

            //~-~-~-~-~-~-~-//
            // Analysis Code//
            //~-~-~-~-~-~-~-//
            Analyze(event, current_file->GetTitle());

        } //end event loop

        // cleanup
        file->Close();

    } // end file loop

    // Stuff to do after job finishes
	EndJob();

    // return
    bmark.Stop("benchmark");
    std::cout << std::endl;
    std::cout << num_events_total << " Events Processed" << std::endl;
    std::cout << num_duplicates   << " Duplicate Events" << std::endl;
    std::cout << num_bad_events   << " Bad Events"       << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "CPU  Time:	" << Form("%.01f", bmark.GetCpuTime("benchmark") ) << std::endl;
    std::cout << "Real Time:	" << Form("%.01f", bmark.GetRealTime("benchmark")) << std::endl;
    std::cout << std::endl;

    // done
    return;
}
