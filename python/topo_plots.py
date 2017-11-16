#!/usr/bin/env python

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from optparse import OptionParser
import sys
import os

class Sample :
    def __init__(self, name, filename, color, signal) :
        self._name = name
        self._filename = filename
        self._color = color
        self._signal = signal
        self._tree = self.get_tree(filename)

    @property
    def name(self) :
        return self._name
    @property
    def filename(self) :
        return self._filename
    @property
    def color(self) :
        return self._color
    @property
    def signal(self) :
        return self._signal
    @property
    def tree(self) :
        return self._tree

    def get_tree(self, filename) :
        chain = r.TChain("trig")
        if not os.path.isfile(filename) :
            print "ERROR file %s does not exist" % filename
            sys.exit()
        chain.AddFile(filename)
        return chain

def make_plots(samples) :

    variables = ["L1_2EM15VH", "L1_MU20", "L1_2MU10", "L1_EM15VH_MU10", "L1_EM22VHI",
            "n_muons", "muon_pt[0]", "abs(muon_phi[0])", "abs(muon_eta[0])",
            "n_electrons", "electron_pt[0]", "abs(electron_phi[0])", "abs(electron_eta[0])",
            "n_jets", "jet_pt[0]", "abs(jet_phi[0])", "abs(jet_eta[0])",
            "met_sum", "met_phi", 
            "dphi_muo_met[0]",
            "dphi_ele_met[0]",
            "dr_muo_jet[0]",
            "dr_ele_jet[0]",
            "ee_mll[0]",
            "ee_dphi[0]",
            "ee_dr[0]",
            "mm_mll[0]",
            "mm_dphi[0]",
            "mm_dr[0]",
            "em_mll[0]",
            "em_dphi[0]",
            "em_dr[0]"]

    bounds = {}
    bounds["L1_2EM15VH"] = {1, 0, 1}
    bounds["L1_MU20"] = {1,0,1}
    bounds["L1_2MU10"] = {1,0,1}
    bounds["L1_EM15VH_MU10"] = {1,0,1}
    bounds["L1_EM22VHI"] = {1,0,1}
    bounds["n_muons"] = {1,0,6}
    bounds["muon_pt[0]"] = {1, 0, 25}
    bounds["abs(muon_phi[0])"] = {0.1, 0, 3.2}
    bounds["abs(muon_eta[0])"] = {0.1, 0, 2.5}
    bounds["n_electrons"] = {1,0,12}
    bounds["electron_pt[0]"] = {5, 0, 60}
    bounds["abs(electron_phi[0])"] = {0.1, 0, 3.2}
    bounds["abs(electron_eta[0])"] = {0.1, 0, 2.5}

def main() :

    sample_dir = "../samples/"
    name_reg = "topotuple_mc_"
    ttbar = Sample("ttbar", "%s%s410009.root" % (sample_dir, name_reg), r.kBlack, False)
    hhsm = Sample("hhSM", "%s%s342053.root" % (sample_dir, name_reg), r.kRed, True)

    samples = [ttbar, hhsm]
    make_plots(samples)
    

if __name__ == "__main__" :
    main()
