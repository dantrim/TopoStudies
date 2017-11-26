#!/usr/bin/env python

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gStyle.SetOptStat(False)
r.gROOT.SetBatch(True)

import os
import sys

from optparse import OptionParser

import numpy as np


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

def make_ac_plot(sample) :

    c = r.TCanvas("c_%s" % sample.name, "", 800, 600)
    c.SetGrid(1,1)
    c.cd()


    lep_pts = list(np.arange(8,32,2))
    n_bins = len(lep_pts)
    bin_low = lep_pts[0]
    bin_high = lep_pts[-1]

    h = r.TH2F("hlep_pt_%s" % sample.name, ";Sub-lead lepton p_{T};Lead lepton p_{T}", n_bins, bin_low, bin_high, n_bins, bin_low, bin_high)

    el_cut_initial = "n_jets>=2 && jet_pt[0]>=12 && jet_pt[1]>=12 && n_electrons>=1 && dr_ele_jet[0]>=2"
    #el_fail_standard = "1"
    #mu_fail_standard = "L1_MU20==0 && L1_2MU10==0 && L1_EM15VH_MU10==0"
    el_fail_standard = "L1_2EM15VH==0 && L1_EM15VH_MU10==0 && L1_EM22VHI==0"
    sel = r.TCut("1")

    for iptcut, ipt in enumerate(lep_pts) :
        for jptcut, jpt in enumerate(lep_pts) :
            if jpt >= ipt : continue
            lead_pt = ipt
            sublead_pt = jpt
            el_cut = el_cut_initial + " && electron_pt[0]>%d && electron_pt[1]>%d" % (ipt, jpt) 
            el_cut = "(%s)" % el_cut
            el_cut = "%s && (%s)" % (el_cut, el_fail_standard)
            el_cut = "(%s)" % el_cut
           # print el_cut
           # sys.exit()
            el_cut = r.TCut(el_cut)

            h_den = r.TH1F("h_den_%s_%d_%d" % (sample.name, iptcut, jptcut), "", 4, 0, 4)
            h_num = r.TH1F("h_num_%s_%d_%d" % (sample.name, iptcut, jptcut), "", 4, 0, 4)

            cmd_den = "isMC>>%s" % h_den.GetName()
            cmd_num = "isMC>>%s" % h_num.GetName()

            sample.tree.Draw(cmd_den, sel, "goff")
            den_counts = h_den.Integral()

            sample.tree.Draw(cmd_num, el_cut * sel, "goff")
            num_counts = h_num.Integral()

            fraction = float(num_counts) / float(den_counts)
            h.Fill(iptcut+1, jptcut+1, fraction)

    r.gStyle.SetPalette(r.kBird)

    h.Draw("colz")
    c.Update()


    c.SaveAs("test_ac.pdf")

    sys.exit()
    
    

def make_ac_plots(samples) :

    for isample, sample in enumerate(samples) :
        if not sample.signal : continue
        make_ac_plot(sample)

def main() :

    sample_dir = "/data/uclhc/uci/user/dantrim/L1TopoStudies/TopoStudies/samples/" 
    name_reg = "topotuple_mc_"
    ttbar = Sample("ttbar", "%s%s410009.root" % (sample_dir, name_reg), r.kBlack, False)
    hhsm = Sample("hhSM", "%s%s342053.root" % (sample_dir, name_reg), r.kRed, False)
    hh300 = Sample("hh300", "%s%s343766.root" % (sample_dir, name_reg), r.kBlue, True)
    hh600 = Sample("hh600", "%s%s343772.root" % (sample_dir, name_reg), r.kMagenta, True)
    hh800 = Sample("hh800", "%s%s343775.root" % (sample_dir, name_reg), 8, True)

    samples = [ttbar, hhsm, hh300, hh600, hh800]

    make_ac_plots(samples)


if __name__ == "__main__" :
    main()
