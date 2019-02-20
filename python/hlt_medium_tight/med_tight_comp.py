#!/usr/bin/env python

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gStyle.SetOptStat(False)
r.gROOT.SetBatch(True)

import sys
import os
import glob

import numpy as np
n_em_single = 0
n_em_asym = 0
n_em_sym = 0


sampledir = "/data/uclhc/uci/user/dantrim/n0234val/" 
dsid_to_name = {}
dsid_to_name["342053"] = "hhSM"
dsid_to_name["343764"] = "X260"
dsid_to_name["343766"] = "X300"
dsid_to_name["343772"] = "X600"
dsid_to_name["410009"] = "ttbar"

class Sample :
    def __init__(self, dsid = "", name = "", files = []) :
        self._dsid = dsid
        self._name = name
        self._files = files

        self._medium_tree = None
        self._tight_tree = None

        self.get_trees(files)

    @property
    def dsid(self) :
        return self._dsid
    @property
    def name(self) :
        return self._name
    @property
    def files(self) :
        return self._files
    @property
    def medium_tree(self) :
        return self._medium_tree
    @property
    def tight_tree(self) :
        return self._tight_tree

    def get_trees(self, list_of_files) :

        for f in list_of_files :
            c = r.TChain("superNt")
            c.AddFile(f)
            if "tight" in f.lower() :
                self._tight_tree = c
            elif "medium" in f.lower() :
                self._medium_tree = c

def get_samples() :

    out = []

    for dsid, name in dsid_to_name.iteritems() :
        files = glob.glob("%s/CENTRAL*%s*.root" % (sampledir, dsid) )
        s = Sample(dsid, name, files)
        out.append(s)

    return out

def get_trigger_cut(ele_id) :

    sE15 = "trig_e24_lhmedium_L1EM20VHI==1 || trig_e60_lhmedium==1"
    dE15 = "trig_2e12_lhloose_L12EM10VH==1"
    
    sM15 = "trig_mu20_iloose_L1MU15==1 || trig_mu40==1"
    dM15 = "trig_mu18_mu8noL1==1"
    
    m15 = "trig_e17_lhloose_mu14==1 || trig_e24_lhmedium_L1EM20VHI_mu8noL1==1 || trig_e7_lhmedium_mu24==1"
    
    sE16 = "trig_e60_lhmedium_nod0==1 || trig_e60_lhmedium==1"
    dE16 = "trig_2e17_lhvloose_nod0==1"
    
    sM16 = "trig_mu26_ivarmedium==1 || trig_mu50==1"
    dM16 = "trig_mu22_mu8noL1==1"
    
    m16 = "trig_e17_lhloose_nod0_mu14==1 || trig_e26_lhmedium_nod0_L1EM22VHI_mu8noL1==1 || trig_e7_lhmedium_nod0_mu24==1"

    if "tight" in ele_id.lower() :
        sE16 += " || trig_e26_lhtight_nod0_ivarloose==1"

#    cut = "(year==2015 && ( %s || %s || %s || %s || %s ))" % ( sE15, dE15, sM15, dM15, m15 )
#    cut += " || "
#    cut += " (year==2016 && ( %s || %s || %s || %s || %s ))" % ( sE16, dE16, sM16, dM16, m16 )
#    cut = "( %s )" % cut

    cut = "(year==2015 && ( %s || %s || %s ))" % (  dE15, dM15, m15 )
    cut += " || "
    cut += " (year==2016 && (  %s || %s || %s ))" % (  dE16,  dM16, m16 )
    cut = "( %s )" % cut

    return cut


def make_lep_pt_2d_for_ele_id(sample, ele_id) :

    c = r.TCanvas('c_%s_%s' % (sample.name, ele_id), '', 800, 600)
    c.SetGrid(1,1)
    c.cd()

    lep_pts = list(np.arange(6,52,2))
    n_bins = len(lep_pts)
    binlow = lep_pts[0]
    binhigh = lep_pts[-1]

    h = r.TH2F("h_%s_%s" % (sample.name, ele_id), "%s %s;Lead lepton p_{T} [GeV];Sub-lead lepton p_{T} [GeV]" % (sample.name, ele_id), n_bins, binlow, binhigh, n_bins, binlow, binhigh)

    xax = h.GetXaxis()
    yax = h.GetYaxis()

    xax.SetTitleFont(42)
    xax.SetLabelFont(42)

    yax.SetTitleFont(42)
    yax.SetLabelFont(42)

    trig_cut_string = get_trigger_cut(ele_id)

    trig_cut = r.TCut(trig_cut_string)
    sel = r.TCut("1")

    draw_cmd = "l1_pt:l0_pt>>%s" % h.GetName()
    if "medium" in ele_id.lower() :
        sample.medium_tree.Draw(draw_cmd, trig_cut * sel, "goff")
    elif "tight" in ele_id.lower() :
        sample.tight_tree.Draw(draw_cmd, trig_cut * sel, "goff")

    h.Scale(1/h.Integral())

    h.Draw("colz")
    c.Update()

    # save

    #c.SaveAs("./med_tight_plots/lep_pt_2d_%s_%s.pdf" % (sample.name, ele_id))
    #c.SaveAs("./med_tight_plots/lep_pt_2d_%s_%s_nosingle.pdf" % (sample.name, ele_id))
    #c.SaveAs("./med_tight_plots/lep_pt_2d_%s_%s_nosingle_nosmartasym.pdf" % (sample.name, ele_id))
    c.SaveAs("./med_tight_plots/lep_pt_2d_%s_%s_nosmartasym.pdf" % (sample.name, ele_id))

def make_lep_pt_1d_for_ele_id(sample, ele_id) :

    c = r.TCanvas("c_%s_%s" % (sample.name, ele_id), "", 800, 600)
    c.SetGrid(1,1)
    c.cd()

    upper = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
    lower = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)
    up_height = 0.75
    dn_height = 0.30
    upper.SetPad(0.0, 1.0 - up_height, 1.0, 1.0)
    lower.SetPad(0.0, 0.0, 1.0, dn_height)
    upper.Draw()
    lower.Draw()
    c.Update()

    upper.cd()

    lep_pts = list(np.arange(10, 100, 5))
    n_bins = len(lep_pts)
    binlow = lep_pts[0]
    binhigh = lep_pts[-1]

    do_sublead = False
    xlabel = "Lead lepton p_{T} [GeV]"
    if do_sublead :
        xlabel = "Sub-lead lepton p_{T} [GeV]"

    h0 = r.TH1F("h_pt_%s_%s" % (sample.name, ele_id), "%s %s; %s;Entries" % (sample.name, ele_id, xlabel), n_bins, binlow, binhigh)
    h0.SetLineColor(r.kBlack)
    htrig = r.TH1F("h_pt_trig_%s_%s" % (sample.name, ele_id), "%s %s; %s;Entries" % (sample.name, ele_id, xlabel), n_bins, binlow, binhigh)
    htrig.SetLineColor(r.kBlue)
    for h in [h0, htrig] :
        xax = h.GetXaxis()
        yax = h.GetYaxis()

        xax.SetTitleFont(42)
        xax.SetLabelFont(42)

        yax.SetTitleFont(42)
        yax.SetLabelFont(42)

    trig_cut_string = get_trigger_cut(ele_id)
    trig_cut = r.TCut(trig_cut_string)
    sel = r.TCut("1")

    draw_cmd = "l0_pt>>%s" % h0.GetName()
    draw_cmd_trig = "l0_pt>>%s" % htrig.GetName()
    if do_sublead :
        draw_cmd = draw_cmd.replace("l0_pt", "l1_pt")
        draw_cmd_trig = draw_cmd_trig.replace("l0_pt", "l1_pt")
    if "medium" in ele_id.lower() :
        sample.medium_tree.Draw(draw_cmd, sel, "goff")
        sample.medium_tree.Draw(draw_cmd_trig, sel * trig_cut, "goff")
    elif "tight" in ele_id.lower() :
        sample.tight_tree.Draw(draw_cmd, sel, "goff")
        sample.tight_tree.Draw(draw_cmd_trig, sel * trig_cut, "goff")

    maxy = h0.GetMaximum()
    h0.SetMaximum(1.2 * maxy)
    htrig.SetMaximum(1.2 * maxy)

    h0.Draw("hist")
    htrig.Draw("hist same")
    upper.Update()
    c.Update()

    lower.cd()

    hratio = htrig.Clone("hratio_%s" % htrig.GetName())
    hratio.Divide(h0)
    hratio.SetMaximum(1.2)
    hratio.SetMinimum(0.0)
    hratio.SetLineColor(r.kBlue)
    hratio.Draw("hist")
    lower.Update()
    c.Update()

    line = r.TLine(binlow, 1.0, binhigh, 1.0)
    line.SetLineStyle(2)
    line.SetLineColor(r.kRed)
    line.Draw()
    lower.Update()
    c.Update()

    save_label = "lead"
    if do_sublead :
        save_label = "sublead"
    c.SaveAs("./med_tight_plots/lep_pt_1d_%s_%s_%s_nosingle.pdf" % (save_label, sample.name, ele_id))
    
            

def make_lep_pt_2d(sample) :

    ele_id = ["medium", "tight"]

    for el in ele_id :
        make_lep_pt_2d_for_ele_id(sample, el)

def make_lep_pt_1d(sample) :
    ele_id = ["medium", "tight"]
    for el in ele_id :
        make_lep_pt_1d_for_ele_id(sample, el)

def pass_ee_logic(lead_pt, sub_pt, evt, ele_id) :

    year = evt.year

    if year == 2015 :
        if (lead_pt >= 26) and (evt.trig_e24_lhmedium_L1EM20VHI==1) :
            return True
        else :
            if (lead_pt>=14) and (sub_pt>=14) and (evt.trig_2e12_lhloose_L12EM10VH==1) :
                return True

    elif year == 2016 :
        if (lead_pt >= 28) and (evt.trig_e26_lhtight_nod0_ivarloose==1) and (ele_id.lower()=="tight") :
            return True
        else :
            if(lead_pt>=19) and (sub_pt>=19) and (evt.trig_2e17_lhvloose_nod0==1) :
                return True

    return False

def pass_mm_logic(lead_pt, sub_pt, evt) :

    year = evt.year

    if year == 2015 :
        if (lead_pt >= 22) and (evt.trig_mu20_iloose_L1MU15==1) :
            return True
        else :
            if (lead_pt>=20) and (sub_pt>=10) and (evt.trig_mu18_mu8noL1==1) :
                return True
    elif year == 2016 :
        if (lead_pt >= 28) and (evt.trig_mu26_ivarmedium==1) :
            return True
        else :
            if (lead_pt>=24) and (sub_pt>=10) and (evt.trig_mu22_mu8noL1==1) :
                return True

    return False

def pass_em_logic(lead_pt, sub_pt, evt, ele_id) :

    global n_em_single, n_em_asym, n_em_sym

    year = evt.year

    if year == 2015 :
        if (lead_pt >= 26) and (evt.trig_e24_lhmedium_L1EM20VHI==1) :
            n_em_single += 1
            return True
        elif (lead_pt >= 26) and (sub_pt >= 10) and (evt.trig_e24_lhmedium_L1EM20VHI_mu8noL1==1) :
            n_em_asym += 1
            return True
        else :
            if (lead_pt >= 20) and (sub_pt >= 17) and (evt.trig_e17_lhloose_mu14==1) :
                n_em_sym += 1
                return True

    elif year == 2016 :
        if ele_id.lower() == "tight" :
            if (lead_pt >= 28) and (evt.trig_e26_lhtight_nod0_ivarloose==1) :
                n_em_single += 1
                return True
        elif ele_id.lower() == "medium" :
            if (lead_pt >= 62) and (evt.trig_e60_lhmedium_nod0==1) :
                n_em_single += 1
                return True
        elif (lead_pt >= 28) and (sub_pt >= 10) and (evt.trig_e26_lhmedium_nod0_L1EM22VHI_mu8noL1==1) :
            n_em_asym += 1
            return True
        else :
            if (lead_pt >= 20) and (sub_pt >= 17) and (evt.trig_e17_lhloose_mu14==1) :
                n_em_syn += 1
                return True

    return False

def pass_me_logic(lead_pt, sub_pt, evt) :

    year = evt.year

    if year == 2015 :
        if (lead_pt >= 22) and (evt.trig_mu20_iloose_L1MU15==1) :
            return True
        elif (lead_pt >= 26) and (sub_pt >= 10) and (evt.trig_e7_lhmedium_mu24==1) :
            return True
        elif (lead_pt >= 17) and (sub_pt >= 20) and (evt.trig_e17_lhloose_mu14==1) :
            return True
        

    elif year == 2016 :
        if (lead_pt >= 28) and (evt.trig_mu26_ivarmedium==1) :
            return True
        elif (lead_pt >= 26) and (sub_pt >= 10) and (evt.trig_e7_lhmedium_nod0_mu24==1) :
            return True
        elif (lead_pt >= 17) and (sub_pt >= 20) and (evt.trig_e17_lhloose_mu14==1) :
            return True

    return False

def pass_simple_OR(evt, ele_id) :

    year = evt.year

    if year == 2015 :
        if ((evt.trig_e24_lhmedium_L1EM20VHI==1) or
            (evt.trig_e60_lhmedium==1) or
            (evt.trig_2e12_lhloose_L12EM10VH==1) or
            (evt.trig_mu20_iloose_L1MU15==1) or
            (evt.trig_mu40==1) or
            (evt.trig_mu18_mu8noL1==1) or
            (evt.trig_e17_lhloose_mu14==1) or
            (evt.trig_e24_lhmedium_L1EM20VHI_mu8noL1==1) or
            (evt.trig_e7_lhmedium_mu24==1)) :
            return True
    elif year == 2016 :

        if ele_id == "medium" :
            if ((evt.trig_e60_lhmedium_nod0==1) or
                (evt.trig_e60_lhmedium==1) or
                (evt.trig_2e17_lhvloose_nod0==1) or
                (evt.trig_mu26_ivarmedium==1) or
                (evt.trig_mu50==1) or
                (evt.trig_mu22_mu8noL1) or
                (evt.trig_e17_lhloose_nod0_mu14==1) or
                (evt.trig_e26_lhmedium_nod0_L1EM22VHI_mu8noL1==1) or
                (evt.trig_e7_lhmedium_nod0_mu24==1)) :
                return True 
        elif ele_id == "tight" :
            if ((evt.trig_e60_lhmedium_nod0==1) or
                (evt.trig_e26_lhtight_nod0_ivarloose==1) or
                (evt.trig_e60_lhmedium==1) or
                (evt.trig_2e17_lhvloose_nod0==1) or
                (evt.trig_mu26_ivarmedium==1) or
                (evt.trig_mu50==1) or
                (evt.trig_mu22_mu8noL1) or
                (evt.trig_e17_lhloose_nod0_mu14==1) or
                (evt.trig_e26_lhmedium_nod0_L1EM22VHI_mu8noL1==1) or
                (evt.trig_e7_lhmedium_nod0_mu24==1)) :
                return True 

    return False

def get_smart_counts_for_ele_id(sample, ele_id) :

    n_den_ee = 0
    n_pass_ee = 0

    n_den_mm = 0
    n_pass_mm = 0

    n_den_em = 0
    n_pass_em = 0

    n_den_me = 0
    n_pass_me = 0

    n_pass_simple = 0

    n_den_current = 0
    n_pass_current = 0

    tree = None
    if "medium" in ele_id.lower() :
        tree = sample.medium_tree
    elif "tight" in ele_id.lower() :
        tree = sample.tight_tree

    for ievent, event in enumerate(tree) :

        is_ee = (event.isEE == 1)
        is_mm = (event.isMM == 1)
        is_em = (event.isEM == 1)
        is_me = (event.isME == 1)

        if is_ee :
            n_den_ee += 1
        elif is_mm :
            n_den_mm += 1
        elif is_em :
            n_den_em += 1
        elif is_me :
            n_den_me += 1

        lead_pt = event.l0_pt
        sublead_pt = event.l1_pt

        year = event.year

        pass_current_pt = ( (lead_pt>=25 and sublead_pt>=20) )
        n_den_current += 1
        if year == 2015 :
            if (event.trig_pass2015 == 1) and pass_current_pt :
                n_pass_current += 1
        elif year == 2016 :
            if (event.trig_pass2016update == 1) and pass_current_pt:
                n_pass_current += 1

        if is_ee and pass_ee_logic(lead_pt, sublead_pt, event, ele_id) :
            n_pass_ee += 1
        elif is_mm and pass_mm_logic(lead_pt, sublead_pt, event) :
            n_pass_mm += 1
        elif is_em and pass_em_logic(lead_pt, sublead_pt, event, ele_id) :
            n_pass_em += 1
        elif is_me and pass_me_logic(lead_pt, sublead_pt, event) :
            n_pass_me += 1

        if pass_simple_OR(event, ele_id) :
            n_pass_simple += 1


    print 50 * '-'
    print " %s %s " % (sample.name, ele_id)
    print "     n_den_ee        : %d" % n_den_ee
    print "     n_pass_ee       : %d  (%.2f)" % ( n_pass_ee, ( (float(n_pass_ee) / float(n_den_ee)) * 100. ) ) 
    print "     n_den_mm        : %d" % n_den_mm
    print "     n_pass_mm       : %d  (%.2f)" % ( n_pass_mm, ( (float(n_pass_mm) / float(n_den_mm)) * 100. ) )
    print "     n_den_em        : %d" % n_den_em
    print "     n_pass_em       : %d  (%.2f)" % ( n_pass_em, ( (float(n_pass_em) / float(n_den_em)) * 100. ) )
    print "     n_den_me        : %d" % n_den_me
    print "     n_pass_me       : %d  (%.2f)" % ( n_pass_me, ( (float(n_pass_me) / float(n_den_me)) * 100. ) )
    print 25 * '- '
    n_den_total = n_den_ee + n_den_mm + n_den_em + n_den_me
    n_pass_total = n_pass_ee + n_pass_mm + n_pass_em + n_pass_me
    print "     n_den_total     : %d" % ( n_den_total ) 
    print "     n_pass_total    : %d  (%.2f)" % ( n_pass_total, ( (float(n_pass_total) / float(n_den_total)) * 100. ) )
    print "     n_den_simple OR : %d" % ( n_den_total )
    print "     n_pass_simple OR: %d  (%.2f)" % ( n_pass_simple, ( (float(n_pass_simple) / float(n_den_total)) * 100. ) )
    print 25 * "= "
    print "     n_den_current   : %d" % n_den_current
    print "     n_pass_current  : %d  (%.2f)" % ( n_pass_current, ( (float(n_pass_current) / float(n_den_current)) * 100. ) )
            
    global n_em_single, n_em_asym, n_em_sym
    print ""
    print 50 * "="
    print " n_pass em_single    : %d" % n_em_single
    print " n_pass em_asym      : %d" % n_em_asym
    print " n_pass em_sym       : %d" % n_em_sym

def get_smart_counts(sample) :

    ele_id = ["medium"]#, "tight"]
    ele_id = ["tight"]#, "tight"]
    #ele_id = ["medium", "tight"]
    for el in ele_id :
        get_smart_counts_for_ele_id(sample, el)

def main() :

    samples = get_samples()

    for sample in samples :
        #make_lep_pt_2d(sample)
        #make_lep_pt_1d(sample)
        get_smart_counts(sample)

    

#________________________________
if __name__ == "__main__" :
    main()
