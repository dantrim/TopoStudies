#!/usr/bin/env python

import ROOT as r
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)
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

def make_plot(samples, var, bounds, var_name) :

    #selection_cut = "L1_MU20==1"
    #selection_name = "L1_MU20"
    #selection_name = "No Selection"
    #selection_cut = "n_jets>=2 && jet_pt[0]>=12 && jet_pt[1]>=12 && n_muons>=1 && muon_pt[0]>=6 && dr_muo_jet[0]>=2"
    #selection_name = "2J12 and MU6 and DRMUJET_2"
    #selection_cut = "L1_EM22VHI" 
    #selection_name = "L1_EM22VHI"
    #selection_cut = "n_jets>=2 && jet_pt[0]>=12 && jet_pt[1]>=12 && n_electrons>=1 && electron_pt[0]>=20 && dr_ele_jet[0]>=2"
    #selection_name = "2J12 and EL20 and DRELJET_2"

    twoj12 = "n_jets>=2 && jet_pt[0]>=12 && jet_pt[1]>=12"
    muX = "n_muons>=1 && muon_pt[0]>=6"
    elX = "n_electrons>=1 && electron_pt[0]>=28"

    selection_cut = twoj12
    selection_name = "2J12"

    var_name_ok = var.replace("[","").replace("]","").replace("abs(", "").replace(")","")

    histos = []
    n_bins = (bounds[2] - bounds[1]) / bounds[0]
    n_bins = int(n_bins)
    maxy = -1
    if "L1_" in var :
        print 65 * '-'
    for sample in samples :
        h = r.TH1F("h_%s_%s" % (sample.name, var_name_ok), ";%s;a.u." % var_name,n_bins, bounds[1], bounds[2])
        h.GetXaxis().SetTitleFont(42)
        h.GetXaxis().SetLabelFont(42)
        h.GetYaxis().SetTitleFont(42)
        h.GetYaxis().SetLabelFont(42)
        if sample.signal :
            h.SetLineStyle(2)
        elif sample.name == "ttbar" :
            h.SetFillColorAlpha(40, 0.2)
        h.SetLineWidth(2)
        h.SetLineColor(sample.color)

        cmd = "%s>>%s" % ( var, h.GetName() )
        cut = r.TCut(selection_cut)
        sel = r.TCut("1")
        sample.tree.Draw(cmd, cut * sel, "goff")

        counts = h.Integral()
        if "L1_" in var :
            print "%s : %.2f" % ( sample.name, counts )

        if h.Integral() == 0 :
            print "ho integral is zero for sample %s and variable %s" % (samples[ih].name, var)
            continue
        h.Scale(1/h.Integral())
        if h.GetMaximum() > maxy : maxy = h.GetMaximum()
    
        histos.append(h)

    c = r.TCanvas("c_%s" % var_name_ok, "", 800, 600)
    c.SetGrid(1,1)
    c.cd()

    leg = r.TLegend(0.67, 0.67, 0.85, 0.85)
    #leg.SetBorderSize(0)

    for ihist, hist in enumerate(histos) :
        hist.SetMaximum(1.3*maxy)
        if ihist == 0 :
            hist.Draw("hist")
        else :
            hist.Draw("hist same")
        c.Update()
        leg.AddEntry(hist, samples[ihist].name, "l")

    leg.Draw()
    c.Update()

    # text
    text = r.TLatex()
    text.SetTextFont(42)
    text.SetTextSize(0.75*text.GetTextSize())
    text.DrawLatexNDC(0.16, 0.85, selection_name)
    text.DrawLatexNDC(0.16, 0.80, selection_cut)

    # save
    save_name = "topo_plot_%s_%s.pdf" % ( ''.join(selection_name.strip().split()), var_name_ok)
    save_out = " >>> Saving plot %s : " % save_name
    save_name = "./topo_roi_plots/%s" % save_name
    save_out += os.path.abspath(save_name)
    c.SaveAs(save_name)
    print save_out
    

def make_eff_plot(samples, lep_trig_den, var, bounds, var_name, cut_left) :
    
    var_name_ok = var.replace("(","").replace(")","").replace("abs(","").replace(")","")
    c = r.TCanvas("c_eff_%s_%s" % (lep_trig_den, var_name_ok), "", 800, 600)
    c.cd()
    upper = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
    #upper.SetGrid(1,1)
    lower = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)

    up_height = 0.75
    dn_height = 0.30
    upper.SetPad(0.0, 1.0 - up_height, 1.0, 1.0)
    lower.SetPad(0.0, 0.0, 1.0, dn_height)

    upper.SetTickx(0)
    upper.SetGridx(1)
    upper.SetGridy(1)
    lower.SetGridx(1)
    lower.SetTicky(0)

    upper.SetFrameFillColor(0)
    upper.SetFillColor(0)

    upper.SetRightMargin(0.05)
    lower.SetRightMargin(0.05)

    upper.SetLeftMargin(0.14)
    lower.SetLeftMargin(0.14)

    upper.SetTopMargin(0.7 * upper.GetTopMargin())

    upper.SetBottomMargin(0.09)
    lower.SetBottomMargin(0.4)

    upper.Draw()
    lower.Draw()
    c.Update()

    # draw stuff on top
    upper.cd()

    selection_tcut = "n_jets>=2 && jet_pt[0]>=12 && jet_pt[1]>=12"
    selection_name = ""

    if "MU" in lep_trig_den :
        selection_tcut += " && n_muons>=1 && muon_pt[0]>=10"
        selection_name = "2J12 and MU10"
    elif "EM" in lep_trig_den :
        selection_tcut += " && n_electrons>=1 && electron_pt[0]>=24"
        selection_name = "2J12 and EM24"

    histos = []
    den_counts = []
    n_bins = (bounds[2] - bounds[1]) / bounds[0]
    n_bins = int(n_bins)
    maxy = -1
    title_offset = 0
    label_offset = 0

    for sample in samples :
        h = r.TH1F("h_%s_%s" % (sample.name, var_name_ok), ";%s;a.u." % var_name, n_bins, bounds[1], bounds[2])
        title_offset = h.GetXaxis().GetTitleOffset()
        label_offset = h.GetXaxis().GetLabelOffset()
        h.GetXaxis().SetTitleFont(42)
        h.GetXaxis().SetLabelFont(42)
        h.GetXaxis().SetLabelOffset(999)
        h.GetXaxis().SetTitleOffset(999)
        h.GetYaxis().SetTitleFont(42)
        h.GetYaxis().SetLabelFont(42)
        h.GetYaxis().SetTitleSize(1.15*h.GetYaxis().GetTitleSize())
        if sample.signal :
            h.SetLineStyle(2)
        elif sample.name == "ttbar" :
            h.SetFillColorAlpha(40, 0.17)
        h.SetLineWidth(2)
        h.SetLineColor(sample.color)

        cmd = "%s>>%s" % (var, h.GetName())
        cut = r.TCut(selection_tcut)
        sel = r.TCut("1")
        sample.tree.Draw(cmd, cut * sel, "goff")

        hd = r.TH1F("h_den_%s_%s" % (sample.name, var_name_ok), "", 20, 0, 20)
        trig = lep_trig_den
        cmd = "%s>>%s" % (lep_trig_den, hd.GetName())
        trig_cut = "%s==1" % lep_trig_den
#        trig_cut = "(%s==1 || L1_EM22VHI==1)" % lep_trig_den
        cut = r.TCut(trig_cut)
#        cut = r.TCut("%s==1" % lep_trig_den)
        sample.tree.Draw(cmd, cut * sel, "goff")
        trig_counts = hd.Integral()
        #print "TRIGGER COUNTS %s for sample %s = %.f" % (lep_trig_den, sample.name, trig_counts)
        #print "    > topo counts = %.2f" % h.Integral()
        den_counts.append(trig_counts)

        #h.Scale(1/h.Integral())
        #if h.GetMaximum() > maxy : maxy = h.GetMaximum()
        histos.append(h)



    ############
    # lower
    lower.cd()

    eff_histos = []
    for isample, sample in enumerate(samples) :
        he = histos[isample].Clone("eff_%s" % histos[isample].GetName())
        he.Reset("ICE")
        he.SetMaximum(2)
        #if "EM" in lep_trig_den :
        #    he.SetMaximum(3.5)
        he.SetMinimum(0)
        he.SetFillStyle(0)
        he.GetXaxis().SetTitleOffset(1.2 * title_offset)
        he.GetXaxis().SetLabelOffset(1.1*label_offset)
        he.GetYaxis().SetTitleSize(3 * he.GetYaxis().GetTitleSize())
        he.GetYaxis().SetLabelSize(2.6 * he.GetYaxis().GetLabelSize())

        arrow = "#downarrow"
        if cut_left :
            arrow = "#uparrow"
        title = "#frac{X}{%s} %s" % (lep_trig_den, arrow)
        he.GetYaxis().SetTitle(title)
        he.GetYaxis().SetTitleOffset(0.33)#* he.GetYaxis().GetTitleOffset())

        he.GetYaxis().SetNdivisions(4)
        he.GetXaxis().SetTitleFont(42)
        he.GetXaxis().SetLabelFont(42)
        he.GetYaxis().SetTitleFont(42)
        he.GetYaxis().SetLabelFont(42)
        he.GetXaxis().SetTitleSize(3.2 * he.GetXaxis().GetTitleSize())
        he.GetXaxis().SetLabelSize(2.6 * he.GetXaxis().GetLabelSize())

        nbins = histos[isample].GetNbinsX()
        for ibin in xrange(nbins) :
            #bc = histos[isample].GetBinContent(ibin+1)
            integral_eff = 0
            if cut_left :
                integral_eff = histos[isample].Integral(0, ibin+1)
            else :
                integral_eff = histos[isample].Integral(ibin+1, -1)
            gain = float(integral_eff) / float(den_counts[isample])
            he.SetBinContent(ibin+1, gain)
        eff_histos.append(he)

    for ihist, hist in enumerate(eff_histos) :
        if ihist == 0:
            hist.Draw("hist")
        else :
            hist.Draw("hist same")
        c.Update()

    line = r.TLine(bounds[1], 1.0, bounds[2], 1.0)
    line.SetLineStyle(3)
    line.SetLineColor(r.kRed)
    line.SetLineWidth(2)
    line.Draw()
    lower.Update()


    # DRAW UPPER
    upper.cd()
    for h in histos :
        h.Scale(1/h.Integral())
        if h.GetMaximum() > maxy : maxy = h.GetMaximum()

    leg = r.TLegend(0.65, 0.65, 0.88, 0.88)
    for ihist, hist in enumerate(histos) :
        hist.SetMaximum(1.3 * maxy)
        if ihist == 0 :
            hist.Draw("hist")
        else :
            hist.Draw("hist same")
        c.Update()
        if samples[ihist].name == "ttbar" :
            leg.AddEntry(hist, samples[ihist].name, "fl")
        else :
            leg.AddEntry(hist, samples[ihist].name, "l")

    leg.Draw()
    c.Update()

    # text
    text = r.TLatex()
    text.SetTextFont(42)
    text.SetTextSize(0.75 * text.GetTextSize())
    text.DrawLatexNDC(0.16, 0.85, selection_name)
    c.Update()

    #########################
    # save
    save_name = "topo_plot_eff_%s_%s_%s.pdf" % ( ''.join(selection_name.strip().split()), lep_trig_den, var_name_ok)
    save_out = " >> Saving plot %s : " % save_name
    save_name = "./topo_roi_eff_plots/%s" % save_name
    save_out += os.path.abspath(save_name)
    c.SaveAs(save_name)
    #leg = r.TLegend(0.65, 0.65, 0.88, 0.88)
    #for ihist, hist in enumerate(histos) :
    #    hist.SetMaximum(1.3 * maxy)
    #    if ihist == 0 :
    #        hist.Draw("hist")
    #    else :
    #        hist.Draw("hist same")
    #    c.Update()
    #    if samples[ihist].name == "ttbar" :
    #        leg.AddEntry(hist, samples[ihist].name, "fl")
    #    else :
    #        leg.AddEntry(hist, samples[ihist].name, "l")

    #leg.Draw()
    #c.Update()
    print save_out

def make_eff_or_plot(samples, lep_trig_den, var, do_left) :

    lep_trig_den_ok = 'or'.join(lep_trig_den)
    var_name_ok = var.replace("(","").replace(")","").replace("abs(","").replace(")","")

    c = r.TCanvas("c_or_eff_%s_%s" % ( lep_trig_den_ok, var_name_ok ), "", 800, 600)
    c.cd()
    upper = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
    middle = r.TPad("middle", "middle", 0.0, 0.0, 1.0, 1.0)
    lower = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)

    up_height = 0.9
    mid_height_low = 0.25
    mid_height_high = 0.40
    dn_height = 0.25
    upper.SetPad(0.0, mid_height_high, 1.0, 1.0)
    middle.SetPad(0.0, mid_height_low, 1.0, mid_height_high)
    lower.SetPad(0.0, 0.0, 1.0, mid_height_low)

    upper.SetTickx(0)
    upper.SetGridx(1)
    upper.SetGridy(1)
    middle.SetGridx(1)
    middle.SetTicky(0)
    lower.SetGridx(1)
    lower.SetTicky(0)

    upper.SetFrameFillColor(0)
    upper.SetFillColor(0)

    right_margin = 0.05
    upper.SetRightMargin(right_margin)
    middle.SetRightMargin(right_margin)
    lower.SetRightMargin(right_margin)

    left_margin = 0.14
    upper.SetLeftMargin(left_margin)
    middle.SetLeftMargin(left_margin)
    lower.SetLeftMargin(left_margin)

    upper.SetBottomMargin(0.04)
    middle.SetBottomMargin(0.15)
    lower.SetBottomMargin(0.47)

    upper.SetTopMargin(0.09)
    middle.SetTopMargin(0.05)
    lower.SetTopMargin(0.02)

    upper.Draw()
    middle.Draw()
    lower.Draw()
    c.Update()

    variables = {}
    nice_name = ""
    bounds = []
    if "dphi" in var :
        variables['e'] = "abs(dphi_ele_met[0])"
        variables['m'] = "abs(dphi_muo_met[0])"
        nice_name = "|#Delta#phi(l,-#sum E_{T})|"
        bounds = [0.1, 0, 3.2]
    elif "dr" in var :
        variables['e'] = "dr_ele_jet[0]"
        variables['m'] = "dr_muo_jet[0]"
        nice_name = "#DeltaR(l,Jj)"
        bounds = [0.2, 0, 7]


    selection_tcut = "n_jets>=2 && jet_pt[0]>=12 && jet_pt[1]>=12"
    ele_pt_threshold = 22
    muo_pt_threshold = 6
    selection_tcut += " && ( (n_electrons>=1 && electron_pt[0]>=%d) || (n_muons>=1 && muon_pt[0]>=%d) )" % (ele_pt_threshold, muo_pt_threshold)
    selection_name = "2J12 and MU%d or EM%d" % (muo_pt_threshold, ele_pt_threshold)

    histos = []
    den_counts = []
    n_bins = (bounds[2] - bounds[1]) / bounds[0]
    n_bins = int(n_bins)
    maxy = -1
    title_offset = 0
    label_offset = 0

    for sample in samples :
        h = r.TH1F("h_%s_%s" % ( sample.name, var_name_ok), ";%s;a.u." % nice_name, n_bins, bounds[1], bounds[2]) 
        title_offset = h.GetXaxis().GetTitleOffset() 
        label_offset = h.GetXaxis().GetLabelOffset()
        h.GetXaxis().SetTitleFont(42)
        h.GetXaxis().SetLabelFont(42)
        h.GetXaxis().SetLabelOffset(999)
        h.GetXaxis().SetTitleOffset(999)
        h.GetYaxis().SetTitleFont(42)
        h.GetYaxis().SetLabelFont(42)

        h.GetYaxis().SetTitleSize(1.15 * h.GetYaxis().GetTitleSize())
        if sample.signal :
            h.SetLineStyle(2)
        elif sample.name == "ttbar" :
            h.SetFillColorAlpha(40, 0.17)
        h.SetLineWidth(2)
        h.SetLineColor(sample.color)

        for ievent, event in enumerate(sample.tree) :
            n_jets = event.n_jets
            if not n_jets >= 2 : continue
            lead_jet_pt = event.jet_pt[0]
            sublead_jet_pt = event.jet_pt[1]
            pass_lj = (lead_jet_pt >= 12)
            pass_sj = (sublead_jet_pt >= 12)
            if not (pass_lj and pass_sj) : continue

            n_ele = event.n_electrons
            n_muo = event.n_muons

            pass_trigger = False
            # pass ele
            if n_ele >= 1 :
                lep_pt = event.electron_pt[0]
                if lep_pt >= ele_pt_threshold :
                    dangle = 0.0
                    if "dphi" in var : dangle = abs(event.dphi_ele_met[0])
                    else : dangle = event.dr_ele_jet[0]
                    h.Fill(dangle)
            # pass muo
            if n_muo >= 1 :
                lep_pt = event.muon_pt[0]
                if lep_pt >= muo_pt_threshold :
                    dangle = 0.0
                    if "dphi" in var : dangle = abs(event.dphi_muo_met[0])
                    else : dangle = event.dr_muo_jet[0]
                    h.Fill(dangle)

        histos.append(h)

        hd = r.TH1F("h_den_%s_%s" % ( sample.name, var_name_ok), "", 20, 0, 20)
        #trig_cut = "(%s==1 || %s==1)" % (lep_trig_den[0], lep_trig_den[1])
        trig_cut = "(%s==1 || %s==1 || L1_2EM15VH==1 || L1_EM15VH_MU10==1)" % (lep_trig_den[0], lep_trig_den[1])
        trig_cut = r.TCut(trig_cut)
        cmd = "%s>>%s" % (lep_trig_den[0], hd.GetName())
        sel = r.TCut("1")
        sample.tree.Draw(cmd, trig_cut * sel, "goff")
        trig_counts = hd.Integral()
        den_counts.append(trig_counts)

    ##############################
    # middle
    middle.cd()
    eff_histos = []
    eff_histo_bkg = None
    for isample, sample in enumerate(samples) :
        he = histos[isample].Clone("eff_%s" % histos[isample].GetName())
        he.Reset("ICE")
        he.SetMaximum(2)
        he.SetMinimum(0)
        he.SetFillStyle(0)
        he.GetXaxis().SetTitleOffset(999) 
        he.GetXaxis().SetLabelOffset(999)
        he.GetYaxis().SetTitleSize(3.2 * he.GetYaxis().GetTitleSize())
        he.GetYaxis().SetLabelSize(4 * he.GetYaxis().GetLabelSize())

        arrow = "#downarrow"
        if do_left :
            arrow = "#uparrow"
        title = "#frac{X}{1L OR 2L} %s" % ( arrow ) 
        he.GetYaxis().SetTitle(title)
        he.GetYaxis().SetTitleOffset(0.33)

        he.GetYaxis().SetNdivisions(4)

        nbins = histos[isample].GetNbinsX()
        for ibin in xrange(nbins) :
            integral_eff = 0
            if do_left :
                integral_eff = histos[isample].Integral(0,ibin+1)
            else :
                integral_eff = histos[isample].Integral(ibin+1,-1)
            gain = float(integral_eff) / float(den_counts[isample])
            he.SetBinContent(ibin+1, gain)
        eff_histos.append(he)
        if "ttbar" in sample.name :
            eff_histo_bkg = he

    for ihist, hist in enumerate(eff_histos) :
        if ihist == 0 :
            hist.Draw("hist")
        else :
            hist.Draw("hist same")
        c.Update()
    line = r.TLine(bounds[1], 1.0, bounds[2], 1.0)
    line.SetLineStyle(3)
    line.SetLineColor(r.kRed)
    line.SetLineWidth(2)
    line.Draw()
    middle.Update()
    c.Update()

    # lower
    lower.cd()
    sob_histos = []
    for ihist, hist in enumerate(eff_histos) :
        if "ttbar" in hist.GetName() : continue
        hs = hist.Clone("sob_%s" % hist.GetName())
        hs.SetMinimum(0.98)
        hs.SetMaximum(2.2)
        hs.GetYaxis().SetTitle("#Delta S/B   ")
        hs.GetXaxis().SetTitleOffset(1.6*title_offset)
        hs.GetXaxis().SetLabelOffset(1.1*label_offset)
        #hs.GetYaxis().SetTitleOffset(0.85 * hs.GetYaxis().GetTitleOffset())
        hs.GetYaxis().SetLabelSize(0.65 * hs.GetYaxis().GetLabelSize())

        hs.GetXaxis().SetTitle(nice_name)
        hs.GetXaxis().SetTitleSize(3.2 * hs.GetXaxis().GetTitleSize())
        hs.GetXaxis().SetLabelSize(2.6 * hs.GetXaxis().GetLabelSize())
        

        hs.Divide(eff_histo_bkg)
        sob_histos.append(hs)

    for ihist, hist in enumerate(sob_histos) :
        if ihist == 0 :
            hist.Draw("hist")
        else :
            hist.Draw("hist same")
        c.Update()

    line2 = r.TLine(bounds[1], 1.0, bounds[2], 1.0)
    line2.SetLineStyle(3)
    line2.SetLineColor(r.kRed)
    line2.SetLineWidth(2)
    line2.Draw()
    middle.Update()

    # UPPER
    upper.cd()
    for h in histos :
        h.Scale(1/h.Integral())
        if h.GetMaximum() > maxy : maxy = h.GetMaximum()
    leg = r.TLegend(0.65, 0.65, 0.88, 0.88)
    for ihist, hist in enumerate(histos) :
        hist.SetMaximum(1.3 * maxy)
        if ihist == 0 :
            hist.Draw("hist")
        else :
            hist.Draw("hist same")
        c.Update()
        if samples[ihist].name == "ttbar" :
            leg.AddEntry(hist, samples[ihist].name, "fl")
        else :
            leg.AddEntry(hist, samples[ihist].name, 'l')
    leg.Draw()
    c.Update()

    # text
    text = r.TLatex()
    text.SetTextFont(42)
    text.SetTextSize(0.75 * text.GetTextSize())
    text.DrawLatexNDC(0.16, 0.85, selection_name)
    c.Update()

    # save
    save_name = "topo_plot_or_eff_%s_%s.pdf" % ( ''.join(selection_name.strip().split()), var_name_ok )
    save_out = " >> Saving plot %s : " % save_name
    save_name = "./topo_roi_or_eff_plots/%s" % save_name
    save_out = os.path.abspath(save_name)
    c.SaveAs(save_name)
    print save_out


def make_plots(samples, effa="", effb="") :

    variables = ["L1_2EM15VH", "L1_MU20", "L1_2MU10", "L1_EM15VH_MU10", "L1_EM22VHI",
            "n_muons", "muon_pt[0]", "abs(muon_phi[0])", "abs(muon_eta[0])",
            "n_electrons", "electron_pt[0]", "abs(electron_phi[0])", "abs(electron_eta[0])",
            "n_jets", "jet_pt[0]", "abs(jet_phi[0])", "abs(jet_eta[0])",
            "met_sum", "abs(met_phi)", 
            "abs(dphi_muo_met[0])",
            "abs(dphi_ele_met[0])",
            "dr_muo_jet[0]",
            "dr_ele_jet[0]",
            "ee_mll[0]",
            "abs(ee_dphi[0])",
            "ee_dr[0]",
            "mm_mll[0]",
            "abs(mm_dphi[0])",
            "mm_dr[0]",
            "em_mll[0]",
            "abs(em_dphi[0])",
            "em_dr[0]"]

    bounds = {}
    bounds["L1_2EM15VH"] = [1, 0, 2]
    bounds["L1_MU20"] = [1,0,2]
    bounds["L1_2MU10"] = [1,0,2]
    bounds["L1_EM15VH_MU10"] = [1,0,2]
    bounds["L1_EM22VHI"] = [1,0,2]
    bounds["n_muons"] = [1,0,6]
    bounds["muon_pt[0]"] = [1, 0, 25]
    bounds["abs(muon_phi[0])"] = [0.1, 0, 3.2]
    bounds["abs(muon_eta[0])"] = [0.1, 0, 2.5]
    bounds["n_electrons"] = [1,0,12]
    bounds["electron_pt[0]"] = [5, 0, 100]
    bounds["abs(electron_phi[0])"] = [0.1, 0, 3.2]
    bounds["abs(electron_eta[0])"] = [0.1, 0, 2.5]
    bounds["n_jets"] = [1, 0, 12]
    bounds["jet_pt[0]"] = [10, 0, 160]
    bounds["abs(jet_phi[0])"] = [0.1, 0, 3.2]
    bounds["abs(jet_eta[0])"] = [0.1, 0, 4.0]
    bounds["met_sum"] = [10, 0, 100]
    bounds["abs(met_phi)"] = [0.1, 0, 3.2]
    bounds["abs(dphi_muo_met[0])"] = [0.1, 0, 3.2]
    bounds["abs(dphi_ele_met[0])"] = [0.1, 0, 3.2]
    bounds["dr_muo_jet[0]"] = [0.2, 0, 7]
    bounds["dr_ele_jet[0]"] = [0.2, 0, 7]
    bounds["ee_mll[0]"] = [10, 0, 120]
    bounds["abs(ee_dphi[0])"] = [0.1, 0, 2]
    bounds["ee_dr[0]"] = [0.2, 0, 4] 
    bounds["mm_mll[0]"] = [10, 0, 120]
    bounds["abs(mm_dphi[0])"] = [0.1, 0, 3.2]
    bounds["mm_dr[0]"] = [0.2, 0, 4]
    bounds["em_mll[0]"] = [10, 0, 120]
    bounds["abs(em_dphi[0])"] = [0.1, 0, 3.2]
    bounds["em_dr[0]"] = [0.2, 0, 4]

    names = {} 
    names["L1_2EM15VH"] = "L1_2EM15VH"
    names["L1_MU20"] = "L1_MU20"
    names["L1_2MU10"] = "L1_2MU10"
    names["L1_EM15VH_MU10"] = "L1_EM15VH_MU10"
    names["L1_EM22VHI"] = "L1_EM22VHI"
    names["n_muons"] = "# muon ROI"
    names["muon_pt[0]"] = "muon ROI p_{T} threshold [GeV]"
    names["abs(muon_phi[0])"] = "|muon ROI #phi|"
    names["abs(muon_eta[0])"] = "|muon ROI #eta|"
    names["n_electrons"] = "# electron ROI"
    names["electron_pt[0]"] = "electron ROI E_{T} [GeV]"
    names["abs(electron_phi[0])"] = "|electron ROI #phi|"
    names["abs(electron_eta[0])"] = "|electron ROI #eta|"
    names["n_jets"] = "# jet ROI"
    names["jet_pt[0]"] = "jet ROI E_{T} [GeV]"
    names["abs(jet_phi[0])"] = "|jet ROI #phi|"
    names["abs(jet_eta[0])"] = "|jet ROI #eta|"
    names["met_sum"] = "-#sum E_{T} [GeV]"
    names["abs(met_phi)"] = "-#sum E_{T} #phi"
    names["abs(dphi_muo_met[0])"] = "|#Delta#phi(#mu, -#sum E_{T})|"
    names["abs(dphi_ele_met[0])"] = "|#Delta#phi(e, -#sum E_{T})|"
    names["dr_muo_jet[0]"] = "#DeltaR (#mu, Jj)"
    names["dr_ele_jet[0]"] = "#DeltaR (e, Jj)"
    names["ee_mll[0]"] = "m_{ee} [GeV]"
    names["abs(ee_dphi[0])"] = "|#Delta#phi_{ee}|"
    names["ee_dr[0]"] = "#DeltaR_{ee}"
    names["mm_mll[0]"] = "m_{#mu #mu} [GeV]"
    names["abs(mm_dphi[0])"] = "|#Delta#phi_{#mu #mu}|"
    names["mm_dr[0]"] = "#DeltaR_{#mu #mu}"
    names["em_mll[0]"] = "m_{e #mu} [GeV]"
    names["abs(em_dphi[0])"] = "|#Delta#phi_{e #mu}|"
    names["em_dr[0]"] = "#DeltaR_{e #mu}"

    n_vars = len(variables)
    effs = []
    if effa!="" :
        effs.append(effa)
    if effb!="" :
        effs.append(effb)

    if len(effs) == 0 :
        for ivar, var in enumerate(variables) :
            print "[%d/%d] %s" % (ivar+1, n_vars, var)
            bound = bounds[var]
            name = names[var]
            make_plot(samples, var, bound, name)
    elif len(effs) == 1 :
        lep_trig_den = effs[0]
        variables = [ "abs(dphi_muo_met[0])",
            "abs(dphi_ele_met[0])",
            "dr_muo_jet[0]",
            "dr_ele_jet[0]"]
        n_vars = len(variables)
        for ivar, var in enumerate(variables) :
            print "[%d/%d] %s" % (ivar+1, n_vars, var)
            do_left = False
            if "dphi" in var :
                do_left = True
            bound = bounds[var]
            name = names[var]
            make_eff_plot(samples, lep_trig_den, var, bound, name, do_left)

    elif len(effs) == 2 :
        lep_trig_den = effs
        variables = [ "dphi_jet", "dr_jet" ]
        n_vars = len(variables)
        for ivar, var in enumerate(variables) :
            print "[%d/%d] %s" % (ivar+1, n_vars, var)
            do_left = False
            if "dphi" in var :
                do_left = True
            make_eff_or_plot(samples, lep_trig_den, var, do_left)

def main() :
    parser = OptionParser()
    parser.add_option("--eff-a", default="", help = "Set a trigger A to compare eff to")
    parser.add_option("--eff-b", default="", help = "Set a trigger B (OR'd with A)")
    (options, args) = parser.parse_args()
    eff_a = options.eff_a
    eff_b = options.eff_b

#    sample_dir = "../samples/"
    sample_dir = "/data/uclhc/uci/user/dantrim/L1TopoStudies/TopoStudies/samples/"
    name_reg = "topotuple_mc_"
    ttbar = Sample("ttbar", "%s%s410009.root" % (sample_dir, name_reg), r.kBlack, False)
    hhsm = Sample("hhSM", "%s%s342053.root" % (sample_dir, name_reg), r.kRed, False)
    hh300 = Sample("hh300", "%s%s343766.root" % (sample_dir, name_reg), r.kBlue, True)
    hh600 = Sample("hh600", "%s%s343772.root" % (sample_dir, name_reg), r.kMagenta, True)
    hh800 = Sample("hh800", "%s%s343775.root" % (sample_dir, name_reg), 8, True)

    samples = [ttbar, hhsm, hh300, hh600, hh800]


    if eff_a != "" and "L1_" not in eff_a :
        print "ERROR requesed eff-a does not have 'L1' in name"
        sys.exit()
    
    make_plots(samples, effa = eff_a, effb = eff_b)
    

if __name__ == "__main__" :
    main()
