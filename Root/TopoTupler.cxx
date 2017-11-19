#include "TopoStudies/TopoTupler.h"

//std/stl
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

//ROOT
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

//xAOD/EDM
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"

//TOLS
#include "TrigConfxAOD/xAODConfigTool.h"

//ANA
using namespace topo;

#define print( ARG ) std::cout << std::setw(20) << std::left << "TopoTupler::" << __FUNCTION__ << "    " << ARG;
#define SET_DUAL_TOOL( TOOLHANDLE, TOOLTYPE, TOOLNAME ) \
    ASG_SET_ANA_TOOL_TYPE(TOOLHANDLE, TOOLTYPE); \
    TOOLHANDLE.setName(TOOLNAME);

const float gev = 1e-3;

et_greater_jet by_jet_et;
et_greater_emtau by_ele_et;
pt_greater_mu by_muo_pt;

//////////////////////////////////////////////////////////////////////////////
TopoTupler::TopoTupler() :
    m_dbg(false),
    m_isMC(false),
    m_dsid(-1),
    m_output_setup(false),
    m_filename(""),
    m_mu_filter(false),
    m_em_filter(false),
    m_all_filter(false),
    m_event(new xAOD::TEvent(xAOD::TEvent::kAthenaAccess) ),
    m_l1_met_ex(-1.0),
    m_l1_met_ey(-1.0),
    m_l1_met_sum(-1.0),
    m_l1_met_phi(-1.0),
    m_output_file(nullptr),
    m_output_tree(nullptr)
{
    print("TopoTupler wakes up\n");
}
//////////////////////////////////////////////////////////////////////////////
float TopoTupler::dphi_mpi_pi(float dphi)
{
    if(dphi >= M_PI) dphi -= 2*M_PI;
    if(dphi < -M_PI) dphi += 2*M_PI;
    return dphi;
}
//////////////////////////////////////////////////////////////////////////////
float TopoTupler::delta_r(float dphi, float deta)
{
    dphi = dphi_mpi_pi(dphi);
    return sqrt( dphi * dphi + deta * deta );
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::Init(TTree *tree)
{
    m_store.clear();
    m_event->readFrom(tree);

#ifdef R21ANA
    if(!m_tdt.isUserConfigured()) {
        m_trigConfTool.setTypeAndName("TrigConf::xAODConfigTool/xAODConfigTool");
        CHECK( m_trigConfTool.retrieve() );

        m_tdt.setTypeAndName("Trig::TrigDecisionTool/TrigDecisionTool");
        CHECK( m_tdt.setProperty("ConfigTool", m_trigConfTool.getHandle()) );
        CHECK( m_tdt.setProperty("TrigDecisionKey", "xTrigDecision") );
        CHECK( m_tdt.retrieve() );
    }
#else
    string toolName = "";
    if(!m_tdt.isUserConfigured()) {
        toolName = "xAODConfigTool";
        SET_DUAL_TOOL(m_trigConfTool, TrigConf::xAODConfigTool, toolName);
        CHECK( m_trigConfTool.retrieve() );

        toolName = "TrigDecisionTool";
        SET_DUAL_TOOL(m_tdt, Trig::TrigDecisionTool, toolName);
        CHECK( m_tdt.setProperty("ConfigTool", m_trigConfTool.getHandle()) );
        CHECK( m_tdt.setProperty("TrigDecisionKey", "xTrigDecision") );
        CHECK( m_tdt.retrieve() );
    }

#endif

    get_metadata();
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::setup_output_tree()
{
    // assume this is called after we have established MC/DATA and DSID

    stringstream ofn;
    ofn << "topotuple_";
    if(is_mc()) ofn << "mc";
    else ofn << "data";
    ofn << "_" << dsid();
    if(em_filter()) ofn << "_EMfilter";
    else if(mu_filter()) ofn << "_MUfilter";
    else if(all_filter()) ofn << "_ALLfilter";
    ofn << ".root";

    m_output_file = new TFile(ofn.str().c_str(), "RECREATE");
    m_output_tree = new TTree("trig", "trig");

    // BRANCHES
    m_output_tree->Branch("isMC", &m_isMC);
    m_output_tree->Branch("dsid", &m_dsid);

    // STANDARD TRIGGERS
    m_output_tree->Branch("L1_2EM15VH", &m_pass_L1_2EM15VH);
    m_output_tree->Branch("L1_MU20", &m_pass_L1_MU20);
    m_output_tree->Branch("L1_2MU10", &m_pass_L1_2MU10);
    m_output_tree->Branch("L1_EM15VH_MU10", &m_pass_L1_EM15VH_MU10);
    m_output_tree->Branch("L1_EM22VHI", &m_pass_L1_EM22VHI);

    // MU STUFF
    m_output_tree->Branch("n_muons", &m_n_muons);
    m_output_tree->Branch("muon_pt", &m_muon_pt);
    m_output_tree->Branch("muon_phi", &m_muon_phi);
    m_output_tree->Branch("muon_eta", &m_muon_eta);

    // ELECTRON STUFF
    m_output_tree->Branch("n_electrons", &m_n_electrons);
    m_output_tree->Branch("electron_pt", &m_electron_pt);
    m_output_tree->Branch("electron_phi", &m_electron_phi);
    m_output_tree->Branch("electron_eta", &m_electron_eta);
    m_output_tree->Branch("electron_emIso", &m_electron_emIso);
    m_output_tree->Branch("electron_hadIso", &m_electron_hadIso);

    // JET STUFF
    m_output_tree->Branch("n_jets", &m_n_jets);
    m_output_tree->Branch("jet_pt", &m_jet_pt);
    m_output_tree->Branch("jet_phi", &m_jet_phi);
    m_output_tree->Branch("jet_eta", &m_jet_eta);
    m_output_tree->Branch("jet_pt6", &m_jet_pt6);
    m_output_tree->Branch("jet_pt8", &m_jet_pt8);

    // MET STUFF
    m_output_tree->Branch("met_ex", &m_l1_met_ex);
    m_output_tree->Branch("met_ey", &m_l1_met_ey);
    m_output_tree->Branch("met_sum", &m_l1_met_sum);
    m_output_tree->Branch("met_phi", &m_l1_met_phi);


    // DANGLE
    m_output_tree->Branch("dphi_muo_met", &m_dphi_muo_met);
    m_output_tree->Branch("dphi_ele_met", &m_dphi_ele_met);
    m_output_tree->Branch("dr_muo_jet", &m_dr_muo_jet);
    m_output_tree->Branch("dr_ele_jet", &m_dr_ele_jet);

    // DILEPTON EE
    m_output_tree->Branch("ee_mll", &m_ee_mll);
    m_output_tree->Branch("ee_dphi", &m_ee_dphi);
    m_output_tree->Branch("ee_dr", &m_ee_dr);

    // DILEPTON MM
    m_output_tree->Branch("mm_mll", &m_mm_mll);
    m_output_tree->Branch("mm_dphi", &m_mm_dphi);
    m_output_tree->Branch("mm_dr", &m_mm_dr);

    // DILEPTON EM
    m_output_tree->Branch("em_mll", &m_em_mll);
    m_output_tree->Branch("em_dphi", &m_em_dphi);
    m_output_tree->Branch("em_dr", &m_em_dr);
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::get_metadata()
{
    //stringstream sx;

    //// data or MC?
    //m_isMC = event()->eventType( xAOD::EventInfo::IS_SIMULATION );
    //sx << "Treating sample as " ? (is_mc() ? "MC" : "DATA");
    //print(sx.str()); sx.str("");

}
//////////////////////////////////////////////////////////////////////////////
Bool_t TopoTupler::Notify()
{
    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::Begin(TTree* /*tree*/)
{
    if(dbg()) print("\n");
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::SlaveBegin(TTree* /*tree*/)
{
    if(dbg()) print("\n");
    m_timer.Start();
}
//////////////////////////////////////////////////////////////////////////////
Bool_t TopoTupler::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo") );

    if(dbg() || chain_entry % 5000 == 0) {
        stringstream sx;
        sx << " **** Processing entry " << setw(6) << chain_entry << " ****\n";
        print(sx.str()); sx.str("");
    }

    m_isMC = ei->eventType( xAOD::EventInfo::IS_SIMULATION );
    if(is_mc()) {
        m_dsid = ei->mcChannelNumber();
    }
    else {
        m_dsid = ei->runNumber();
    }

    if(!m_output_setup) {

        m_output_setup = true;
        setup_output_tree();
    }

    if(mu_filter() || em_filter() || all_filter()) {
        if(!pass_filter()) return kTRUE; 
    }


    fill_l1_roi_objects();

    fill_tree_muon();
    fill_tree_electron();
    fill_tree_jet();
    fill_tree_standard_triggers();
    fill_tree_dangle();
    fill_tree_dilepton_ee();
    fill_tree_dilepton_mm();
    fill_tree_dilepton_em();


    m_output_tree->Fill();
    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_l1_roi_objects()
{

    ///////////////////////////////////////
    // MUON ROI
    ///////////////////////////////////////
    m_l1_muons.clear();
    const xAOD::MuonRoIContainer* l1_muon = 0;
    CHECK( event()->retrieve( l1_muon, "LVL1MuonRoIs" ));
    for(auto & muon : * l1_muon) {
        m_l1_muons.push_back(muon);
    }
    std::sort(m_l1_muons.begin(), m_l1_muons.end(), by_muo_pt);

    ///////////////////////////////////////
    // EM TAU ELECTRONS
    ///////////////////////////////////////
    m_l1_electrons.clear();
    const xAOD::EmTauRoIContainer* l1_emtau = 0;
    CHECK( event()->retrieve( l1_emtau, "LVL1EmTauRoIs") );
    for(auto & emtau : * l1_emtau) {
        if(emtau->roiType()==2) continue; // skip Tau roi
        m_l1_electrons.push_back(emtau);
    }
    std::sort(m_l1_electrons.begin(), m_l1_electrons.end(), by_ele_et);

    ///////////////////////////////////////
    // JET ROI
    ///////////////////////////////////////
    m_l1_jets.clear();
    const xAOD::JetRoIContainer* l1_jet = 0;
    CHECK( event()->retrieve( l1_jet, "LVL1JetRoIs") );
    for(auto & jet : * l1_jet) {
        m_l1_jets.push_back(jet);
    }
    std::sort(m_l1_jets.begin(), m_l1_jets.end(), by_jet_et);

    ///////////////////////////////////////
    // MET
    ///////////////////////////////////////
    const xAOD::EnergySumRoI* l1_e_sum = 0;
    CHECK( event()->retrieve( l1_e_sum, "LVL1EnergySumRoI") );
    float ex = l1_e_sum->exMiss();
    float ey = l1_e_sum->eyMiss();
    m_l1_met_ex = ex * gev;
    m_l1_met_ey = ey * gev;
    m_l1_met_sum = sqrt( ex * ex + ey * ey ) * gev;
    m_l1_met_phi = atan2(ey, ex);

}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_muon()
{

    m_n_muons = m_l1_muons.size();
    m_muon_pt.clear();
    m_muon_phi.clear();
    m_muon_eta.clear();
    for(auto & mu : m_l1_muons) {
        m_muon_pt.push_back(mu->thrValue() * gev);
        m_muon_phi.push_back(mu->phi());
        m_muon_eta.push_back(mu->eta());
    }

}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_electron()
{
    m_n_electrons = m_l1_electrons.size();
    m_electron_pt.clear();
    m_electron_phi.clear();
    m_electron_eta.clear();
    m_electron_emIso.clear();
    m_electron_hadIso.clear();
    for(auto & el : m_l1_electrons) {
        m_electron_pt.push_back(el->eT() * gev);
        m_electron_phi.push_back(el->phi());
        m_electron_eta.push_back(el->eta());
        m_electron_emIso.push_back(el->emIsol() * gev);
        m_electron_hadIso.push_back(el->hadIsol() * gev);
    }
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_jet()
{
    m_n_jets = m_l1_jets.size();
    m_jet_pt.clear();
    m_jet_phi.clear();
    m_jet_eta.clear();
    m_jet_pt6.clear();
    m_jet_pt8.clear();
    for(auto & jet : m_l1_jets) {
        m_jet_pt.push_back(jet->et4x4() * gev);
        m_jet_phi.push_back(jet->phi());
        m_jet_eta.push_back(jet->eta());
        m_jet_pt6.push_back(jet->et6x6() * gev);
        m_jet_pt8.push_back(jet->et8x8() * gev);
    }

}
//////////////////////////////////////////////////////////////////////////////
bool TopoTupler::pass_filter()
{
    string trigger = "";
    bool pass = false;
    if(em_filter()) {
        trigger = "HLT_noalg_L1EM10VH";
        //trigger = "HLT_noalg_L1EM20VH";
        //auto allChains = m_tdt->getChainGroup("HLT_noalg_.*");
        //stringstream sx;
        //cout << "--------------------------------" << endl;
        //for(auto & trig : allChains->getListOfTriggers()) {
        //    if(m_tdt->isPassed(trig))  sx << " " << trig << "\n";
        //    //sx << "  " << trig << " [" << (m_tdt->isPassed(trig) ? "1" : "0") << "]";
        //}
        //cout << sx.str() << endl;
        if(m_tdt->isPassed(trigger)) {
            pass = true;
            if(dbg()) cout << "PASS EM FILTER " << trigger << endl;
        }
    }
    else if(mu_filter()) {
        trigger = "HLT_noalg_L1MU4";
        if(m_tdt->isPassed(trigger)) {
            pass = true;
            if(dbg()) cout << "PASS MU FILTER " << trigger << endl;
        }
    }
    else if(all_filter()) {
        trigger = "HLT_noalg_.*";
        auto allChains = m_tdt->getChainGroup(trigger);
        for(auto & trig : allChains->getListOfTriggers()) {
            if(m_tdt->isPassed(trig)) { pass = true; break; }
        }
    }
    return pass;
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_standard_triggers()
{
    string trigger = "";

    // check mumu L1 seed
    trigger = "L1_MU20";
    bool pass_L1_MU20 = false;
    if(m_tdt->isPassed(trigger)) { pass_L1_MU20 = true; }
    m_pass_L1_MU20 = (pass_L1_MU20 ? 1 : 0);

    trigger = "L1_2MU10";
    bool pass_L1_2MU10 = false;
    if(m_tdt->isPassed(trigger)) { pass_L1_2MU10 = true; }
    m_pass_L1_2MU10 = (pass_L1_2MU10 ? 1 : 0);

    // check ee L1 seed
    trigger = "L1_2EM15VH";
    bool pass_L1_2EM15VH = false;
    if(m_tdt->isPassed(trigger)) { pass_L1_2EM15VH = true; }
    m_pass_L1_2EM15VH = (pass_L1_2EM15VH ? 1 : 0);

    // check emu L1 seed
    trigger = "L1_EM15VH_MU10";
    bool pass_L1_EM15VH_MU10 = false;
    if(m_tdt->isPassed(trigger)) { pass_L1_EM15VH_MU10 = true; }
    m_pass_L1_EM15VH_MU10 = (pass_L1_EM15VH_MU10 ? 1: 0);

    // single e
    trigger = "L1_EM22VHI";
    bool pass_L1_EM22VHI = false;
    if(m_tdt->isPassed(trigger)) { pass_L1_EM22VHI = true; }
    m_pass_L1_EM22VHI = (pass_L1_EM22VHI ? 1 : 0);
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_dangle()
{
    // deltas for the muons
    m_dphi_muo_met.clear();
    m_dr_muo_jet.clear();

    for(int imu = 0; imu < (int)m_l1_muons.size(); imu++) {
        for(int ijet = 0; ijet < (int)m_l1_jets.size(); ijet++) {
            float dphi = m_l1_muons.at(imu)->phi() - m_l1_jets.at(ijet)->phi();
            float deta = m_l1_muons.at(imu)->eta() - m_l1_jets.at(ijet)->eta();
            float dr = delta_r(dphi, deta);
            m_dr_muo_jet.push_back(dr);
        } // ijet
        float dphi = m_l1_muons.at(imu)->phi() - m_l1_met_phi;
        dphi = std::abs(dphi_mpi_pi(dphi));
        m_dphi_muo_met.push_back(dphi);
    } // imu
    std::sort(m_dphi_muo_met.begin(), m_dphi_muo_met.end());
    std::sort(m_dr_muo_jet.begin(), m_dr_muo_jet.end());
    //cout << "******************************" << endl;
    //cout << "muo dr jet before = ";
    //for(auto dr : m_dr_muo_jet) {
    //    cout << dr << " ";
    //}
    //cout << endl;
    std::reverse(m_dr_muo_jet.begin(), m_dr_muo_jet.end());
    //cout << "muo dr jet after = ";
    //for(auto dr : m_dr_muo_jet) {
    //    cout << dr << " ";
    //}
    //cout << endl;

    // delta for the electrons
    m_dphi_ele_met.clear();
    m_dr_ele_jet.clear();
    for(int iel = 0; iel < (int)m_l1_electrons.size(); iel++) {
        for(int ijet = 0; ijet < (int)m_l1_jets.size(); ijet++) {
            float dphi = m_l1_electrons.at(iel)->phi() - m_l1_jets.at(ijet)->phi();     
            float deta = m_l1_electrons.at(iel)->eta() - m_l1_jets.at(ijet)->eta();
            float dr = delta_r(dphi, deta);
            m_dr_ele_jet.push_back(dr);
        } // ijet
        float dphi = m_l1_electrons.at(iel)->phi() - m_l1_met_phi;
        dphi = std::abs(dphi_mpi_pi(dphi));
        m_dphi_ele_met.push_back(dphi);
    } // iel
    std::sort(m_dphi_ele_met.begin(), m_dphi_ele_met.end());
    std::sort(m_dr_ele_jet.begin(), m_dr_ele_jet.end());
    std::reverse(m_dr_ele_jet.begin(), m_dr_ele_jet.end());
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_dilepton_ee()
{
    m_ee_mll.clear();
    m_ee_dphi.clear();
    m_ee_dr.clear();

    for(int iel = 0; iel < (int)m_l1_electrons.size(); iel++) {
        for(int jel = iel+1; jel < (int)m_l1_electrons.size(); jel++) {
            float dphi = m_l1_electrons.at(iel)->phi() - m_l1_electrons.at(jel)->phi();
            float deta = m_l1_electrons.at(iel)->eta() - m_l1_electrons.at(jel)->eta();
            float dr = delta_r(dphi, deta);
            dphi = std::abs(dphi_mpi_pi(dphi));
            m_ee_dphi.push_back(dphi);
            m_ee_dr.push_back(dr);
            float cosh_deta = cosh(deta);
            float cos_dphi = cos(dphi);
            float e0 = m_l1_electrons.at(iel)->eT();
            float e1 = m_l1_electrons.at(jel)->eT();
            float m = 2. * e0 * e1 * (cosh_deta - cos_dphi);
            m = sqrt(m) * gev;
            m_ee_mll.push_back(m);
        } // jel
    } // iel
    std::sort(m_ee_dphi.begin(), m_ee_dphi.end());
    std::sort(m_ee_dr.begin(), m_ee_dr.end());
    std::sort(m_ee_mll.begin(), m_ee_mll.end());
    std::reverse(m_ee_mll.begin(), m_ee_mll.end());
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_dilepton_mm()
{
    m_mm_mll.clear();
    m_mm_dphi.clear();
    m_mm_dr.clear();

    for(int imu = 0; imu < (int)m_l1_muons.size(); imu++) {
        for(int jmu = imu+1; jmu < (int)m_l1_muons.size(); jmu++) {
            float dphi = m_l1_muons.at(imu)->phi() - m_l1_muons.at(jmu)->phi();
            float deta = m_l1_muons.at(imu)->eta() - m_l1_muons.at(jmu)->eta();
            float dr = delta_r(dphi, deta);
            dphi = std::abs(dphi_mpi_pi(dphi));
            m_mm_dphi.push_back(dphi);
            m_mm_dr.push_back(dr);
            float cosh_deta = cosh(deta);
            float cos_dphi = cos(dphi);
            float e0 = m_l1_muons.at(imu)->thrValue();
            float e1 = m_l1_muons.at(jmu)->thrValue();
            float m = 2. * e0 * e1 * (cosh_deta - cos_dphi);
            m = sqrt(m) * gev;
            m_mm_mll.push_back(m);
        } // jmu
    } // imu
    std::sort(m_mm_dphi.begin(), m_mm_dphi.end());
    std::sort(m_mm_dr.begin(), m_mm_dr.end());
    std::sort(m_mm_mll.begin(), m_mm_mll.end());
    std::reverse(m_mm_mll.begin(), m_mm_mll.end());
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::fill_tree_dilepton_em()
{
    m_em_mll.clear();
    m_em_dphi.clear();
    m_em_dr.clear();

    for(int iel = 0; iel < (int)m_l1_electrons.size(); iel++) {
        for(int imu = 0; imu < (int)m_l1_muons.size(); imu++) {
            float dphi = m_l1_electrons.at(iel)->phi() - m_l1_muons.at(imu)->phi();
            float deta = m_l1_electrons.at(iel)->eta() - m_l1_muons.at(imu)->eta();
            float dr = delta_r(dphi, deta);
            dphi = std::abs(dphi_mpi_pi(dphi));
            m_em_dphi.push_back(dphi);
            m_em_dr.push_back(dr);
            float cosh_deta = cosh(deta);
            float cos_dphi = cos(dphi);
            float e0 = m_l1_electrons.at(iel)->eT();
            float e1 = m_l1_muons.at(imu)->thrValue();
            float m = 2. * e0 * e1 * (cosh_deta - cos_dphi);
            m = sqrt(m) * gev;
            m_em_mll.push_back(m);
        } // imu
    } // iel
    std::sort(m_em_dphi.begin(), m_em_dphi.end());
    std::sort(m_em_dr.begin(), m_em_dr.end());
    std::sort(m_em_mll.begin(), m_em_mll.end());
    std::reverse(m_em_mll.begin(), m_em_mll.end());
}
//////////////////////////////////////////////////////////////////////////////
void TopoTupler::Terminate()
{
    stringstream sx;
    m_output_file->cd();
    m_output_tree->Write(0, TObject::kOverwrite);
    sx << "----------------------------------------------------\n";
    print(sx.str()); sx.str("");
    sx << "Output tree saved to : " << m_output_file->GetName() << "\n";
    print(sx.str()); sx.str("");
    sx << "----------------------------------------------------\n";
    print(sx.str()); sx.str("");
    m_output_file->Close();

}
//////////////////////////////////////////////////////////////////////////////
TopoTupler::~TopoTupler()
{
    print(" ~TopoTupler sleeps now~");
}
