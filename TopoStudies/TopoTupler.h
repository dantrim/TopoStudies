#ifndef TOPO_TOPOSTUDIES_H
#define TOPO_TOPOSTUDIES_H

//std/stl
#include <vector>
#include <string>

//infrastructure
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "AsgTools/AnaToolHandle.h"
#include "xAODTrigger/EnergySumRoI.h"
#include "xAODTrigger/EnergySumRoIAuxInfo.h"

//TOOLS
namespace TrigConf {
    class ITrigConfigTool;
}
#include "TrigDecisionTool/TrigDecisionTool.h"


//ROOT
#include "TSelector.h"
#include "TStopwatch.h"
#include "TFile.h"

#undef CHECK
#define CHECK( ARG ) \
    do { \
        const bool result = ARG; \
        if( ! result ) { \
            ::Error("TopoTupler", "Failed to execute: \"%s\"", \
                #ARG ); \
            exit(1); \
        } \
    } while (false)

namespace topo {

    struct et_greater_jet {
        bool operator()(const xAOD::JetRoI* a, const xAOD::JetRoI* b) { return a->et4x4() > b->et4x4(); }
    };

    struct et_greater_emtau {
        bool operator()(const xAOD::EmTauRoI* a, const xAOD::EmTauRoI* b) { return a->eT() > b->eT(); }
    };

    struct pt_greater_mu {
        bool operator()(const xAOD::MuonRoI* a, const xAOD::MuonRoI* b) { return a->thrValue() > b->thrValue(); }
    };

    class TopoTupler : public TSelector
    {
        public :
            TopoTupler();
            virtual ~TopoTupler();

            void set_dbg(bool doit) { m_dbg = doit; }
            bool dbg() { return m_dbg; }

            bool is_mc() { return m_isMC; }
            int dsid() { return m_dsid; }
            void set_filename(std::string name) { m_filename = name; }
            std::string filename() { return m_filename; }

            void set_em_filter(bool doit) { m_em_filter = doit; }
            bool em_filter() { return m_em_filter; }
            void set_mu_filter(bool doit) { m_mu_filter = doit; }
            bool mu_filter() { return m_mu_filter; }

            void get_metadata();
            void setup_output_tree();
            float dphi_mpi_pi(float dphi);
            float delta_r(float dphi, float deta);

            xAOD::TEvent* event() { return m_event; }

            bool pass_filter();

            void fill_l1_roi_objects();
            void fill_tree_muon();
            void fill_tree_electron();
            void fill_tree_jet();
            void fill_tree_standard_triggers();
            void fill_tree_dangle();
            void fill_tree_dilepton_ee();
            void fill_tree_dilepton_mm();
            void fill_tree_dilepton_em();

            // TSELECTOR
            virtual Int_t Version() const { return 2; }
            virtual void Init(TTree* tree);
            virtual Bool_t Notify();
            virtual void Begin(TTree* tree);
            virtual void SlaveBegin(TTree* tree);
            virtual void Terminate();
            virtual Bool_t Process(Long64_t entry);

        private :
            bool m_dbg;
            bool m_isMC;
            int m_dsid;
            bool m_output_setup;
            std::string m_filename;
            bool m_mu_filter;
            bool m_em_filter;
            TStopwatch m_timer;
            xAOD::TEvent* m_event;
            xAOD::TStore m_store;

            // event-wise object containers
            std::vector<xAOD::MuonRoI*> m_l1_muons;
            std::vector<xAOD::EmTauRoI*> m_l1_electrons;
            std::vector<xAOD::JetRoI*> m_l1_jets;
            float m_l1_met_ex; // gev
            float m_l1_met_ey; // gev
            float m_l1_met_sum; // gev
            float m_l1_met_phi; // gev

            asg::AnaToolHandle<TrigConf::ITrigConfigTool> m_trigConfTool;
            asg::AnaToolHandle<Trig::TrigDecisionTool> m_tdt;

            // OUTPUT
            TFile* m_output_file;
            TTree* m_output_tree;

            // BRANCHES
            int m_n_muons;
            std::vector<int> m_muon_pt;
            std::vector<float> m_muon_phi;
            std::vector<float> m_muon_eta;

            int m_n_electrons;
            std::vector<float> m_electron_pt;
            std::vector<float> m_electron_phi;
            std::vector<float> m_electron_eta;
            std::vector<float> m_electron_emIso;
            std::vector<float> m_electron_hadIso;

            int m_n_jets;
            std::vector<float> m_jet_pt;
            std::vector<float> m_jet_phi;
            std::vector<float> m_jet_eta;
            std::vector<float> m_jet_pt6;
            std::vector<float> m_jet_pt8;

            int m_pass_L1_2EM15VH;
            int m_pass_L1_MU20;
            int m_pass_L1_2MU10;
            int m_pass_L1_EM15VH_MU10;
            int m_pass_L1_EM22VHI;

            std::vector<float> m_dphi_muo_met;
            std::vector<float> m_dphi_ele_met;
            std::vector<float> m_dr_muo_jet;
            std::vector<float> m_dr_ele_jet;

            // dilepton
            std::vector<float> m_ee_mll;
            std::vector<float> m_ee_dphi;
            std::vector<float> m_ee_dr;

            std::vector<float> m_mm_mll;
            std::vector<float> m_mm_dphi;
            std::vector<float> m_mm_dr;

            std::vector<float> m_em_mll;
            std::vector<float> m_em_dphi;
            std::vector<float> m_em_dr;

    }; // class TopoTupler

} // namespace topo



#endif
