###############################################
# Mt2dc analysis
# Purpose: To analysis the results of the m2dc
#    variable calculation to determine its properties,
#    it benefits, and its drawbacks.  Comparisons are
#    made with the original mt2 variable.  Final stylized
#    plots are not produced here but by mt2dc_makeFinalPlots.py.
#    This code executes the mt2dc calculation, performed 
#    by mt2dc.py.
# 
# 
# 
# Authors:  
#    Ewan Hill    <ewan.hill@utoronto.ca>
# 
# 
# Notes:
#    * The input and output files must be set in the 
#      code before it can be properly run.
# 
# Usage:
#    Before the code can be run, you must set up root.
#
#    Example of execution:
#       python mt2dc_analysis_v01.py
#           or if you want to save the output to a log file:
#       python mt2dc_analysis_v01.py >& out.log
#   
#   
#   
##############################################

import ROOT
import mt2dc as DC
import math



# With newer versions of root
# from ROOT import gROOT
# gROOT.LoadMacro("AtlasStyle.C")



####################################
# Define the input and output root files
####################################
f_inputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_03_March_07_skim_mg5_ttbar_jet_merged_001-716_ntuple_2l2b_v01.root", "read")

t = f_inputRoot.Get("variables")
type(t)

f_outputRoot = ROOT.TFile.Open("/Users/juliakim/Documents/2022_05_May_09_mt2dc_analysis_v01.root", "recreate")



####################################
# Define constants
####################################
m_W = 80.   # GeV - NEED TO LOOK UP MASS USED IN THE TTBAR GENERATOR !!!!!!!!!!!
m_t = 173.  # GeV - NEED TO LOOK UP MASS USED IN THE TTBAR GENERATOR !!!!!!!!!!!

####################################
# Define the plots to produce
####################################

h_ell1_pt = ROOT.TH1F("h_ell1_pt", "Pt of highest pt light lepton;Leading light lepton pt [GeV];Number of entries / 2 GeV", 100, 0, 200)
h_ell1_E = ROOT.TH1F("h_ell1_E", "E of highest pt light lepton;Leading light lepton Energy    [GeV];Number of entries / 2.5 GeV", 100, 0, 250)
h_bjet1_E = ROOT.TH1F("h_bjet1_E", "E of highest pt b-tagged jet;Leading b-tagged jet Energy    [GeV];Number of entries / 5 GeV", 100, 0, 500)
h_mT2_W = ROOT.TH1F("h_mT2_W", "mt2(ell1,ell2) = mt2(W);mt2(W)   [GeV];Number of entries / 1 GeV", 200, 0, 200)
h_mT2_t_11_22 = ROOT.TH1F("h_mT2_t_11_22", "mt2|11,22(t) = mt2(b1 ell1,b2 ell2);mt2(t|11,22)   [GeV];Number of entries / 1 GeV", 300, 0, 300)
h_mT2_t_12_21 = ROOT.TH1F("h_mT2_t_12_21", "mt2|12,21(t) = mt2(b1 ell2,b2 ell1);mt2(t|12,21)   [GeV];Number of entries / 1 GeV", 300, 0, 300)
h_mT2_t_min   = ROOT.TH1F("h_mT2_t_min",   "min( mt2(t|11,22), mt2(t|12,21) );mt2(t)_min   [GeV];Number of entries / 1 GeV",     300, 0, 300)

h_mT_ell1__forMt2Overlay = ROOT.TH1F("h_mT_ell1_forMt2Overlay",   "mt(ell1, EtMiss);mt(ell1)   [GeV];Number of entries / 1 GeV",     300, 0, 300)
# h_mT_ell2 = ROOT.TH1F("h_mT_ell2",   "mt(ell2, EtMiss);mt(ell2)   [GeV];Number of entries / 1 GeV",     300, 0, 300)
h_mT_ell1 = ROOT.TH1F("h_mT_ell1",   "mt(ell1, EtMiss);mt(ell1)   [GeV];Number of entries / 5 GeV",     200, 0, 1000)


# my own histogram 
h_ell2_pt = ROOT.TH1F("h_ell2_pt", "Pt of lowest pt light lepton; Second light lepton pt [GeV]; Number of entries / 2GeV", 100, 0, 200)


# Make the below a list with different alpha values !!!!!!!!!!!!!!!!
# h_mT2dc
h_muT2dc = ROOT.TH1F("h_muT2dc", "muT2DC(theta = 1.15 radians) = ( sin(1.15) mT2'(W) + cos(1.15) mT2'(t) ) / ( sin(1.15) m(W) + cos(1.15) m(t) );muT2DC   [unitless];Number of entries / 0.01", 300, 0, 3)
# histogram for transverse mass variable 


h_mT2prime_W = ROOT.TH1F("h_mT2prime_W", "mt2'(W);mt2'(W)   [GeV];Number of entries / 1 GeV", 200, 0, 200)
# h_mT2prime_t = ROOT.TH1F("h_mT2prime_t", "mt2'(t);mt2'(t)   [GeV];Number of entries / 1 GeV", 200, 0, 200)
# h_mT2dc_delta
# h_mT2dc_delta_norm


h_muT2_Wt = ROOT.TH1F("h_muT2_Wt", "muT2_Wt(theta = 1.15 radians) = ( sin(1.15) mT2(W) + cos(1.15) mT2(t) ) / ( sin(1.15) m(W) + cos(1.15) m(t) );muT2_Wt   [unitless];Number of entries / 0.01", 300, 0, 3)

# h2_alpha__vs__mT2dc


####################################
# Define other useful variables
####################################
nentries = t.GetEntries() # 60599 

# total # of events 


####################################
# Main analysis - loop over all events
# >>>>>>>>>>>>>>>> Main analysis
####################################



# Note: /Users/ewanhill/Documents/postdoc_ATLAS_UofT/fastDetectorSimulation/Delphes-3.5.0/classes/DelphesClasses.cc
# The mass is zero for electrons and muons in Delphes.  Not jets though.

for i in range(nentries):
    # Print status update to user every 1000 events
    if (   ( i % 1000 == 0 )   ): # every 1000th event 
       print(":: Processing entry ", i, " = ", i*1.0/nentries*100.0, "%.")    
       # Use *1.0 to force float division
       # percentage accomplished 
    # load tree
    if t.LoadTree(i) < 0:
       print("**could not load tree for entry #%s") % i
       break
    # load entry
    nb = t.GetEntry(i) #t. = access info. from tree 
    if nb <= 0:
       # no data
       continue

    # Fill 4-vector information for electron/muon 1 (highest pt light lepton)
    pt  = t.ell1_PT # single entry 
    eta = t.ell1_Eta
    phi = t.ell1_Phi
    m = 0
    p4_ell1 = ROOT.TLorentzVector() # 4-vec 
    p4_ell1.SetPtEtaPhiM(pt,eta,phi,m)

    h_ell1_pt.Fill( p4_ell1.Pt() ) # put info into ROOT histogram 
    h_ell1_E.Fill( p4_ell1.E() )
    
    
    # Fill 4-vector information for electron/muon (2nd highest p light lepton) 
    pt = t.ell2_PT
    eta = t.ell2_Eta
    phi = t.ell2_Phi
    m = 0
    p4_ell2 = ROOT.TLorentzVector() 
    p4_ell2.SetPtEtaPhiM(pt, eta, phi, m)
    
    h_ell1_pt.Fill(p4_ell1.Pt()) 
    h_ell1_E.Fill(p4_ell1.E())


    # Fill 4-vector information for bjet 1 (highest pt jet originating from a b-quark)
    pt  = t.bjet1_PT
    eta = t.bjet1_Eta
    phi = t.bjet1_Phi
    m = t.bjet1_Mass
    p4_bjet1 = ROOT.TLorentzVector()
    p4_bjet1.SetPtEtaPhiM(pt,eta,phi,m)

    h_bjet1_E.Fill( p4_bjet1.E() )



    # Fill 4-vector information for missing transverse energy
    pt  = t.EtMiss
    eta = 0
    phi = t.EtMiss_phi
    m = 0
    p4_met = ROOT.TLorentzVector()
    p4_met.SetPtEtaPhiM(pt,eta,phi,m)

    mt_ell1 = DC.mT_4vecCalc(p4_ell1, p4_met)
    h_mT_ell1.Fill( mt_ell1 )
    h_mT_ell1__forMt2Overlay.Fill( mt_ell1 )



    # Read in mt2 variables
    mt2_W = t.mt2_W_ell1ell2
    mt2_t_11_22 = t.mt2_t_bjet1ell1_bjet2ell2
    mt2_t_12_21 = t.mt2_t_bjet1ell2_bjet2ell1
    mt2_t_min = min(mt2_t_11_22, mt2_t_12_21)

   # Fill histos. with info
    h_mT2_W.Fill( mt2_W )
    h_mT2_t_11_22.Fill( mt2_t_11_22 )
    h_mT2_t_12_21.Fill( mt2_t_12_21 )
    h_mT2_t_min.Fill( mt2_t_min )


    # Calculate mt2dc stuff here............
    ###############
    # filler
    h_mT2prime_W.Fill( mt2_W + 10 )   # Fill histogram with junk for now.

    
    theta = 1.15
    maxMassTheta = math.cos(theta) * m_t + math.sin(theta) * m_W 
    if (   ( i % 500 == 0 )   ):
       print("maxMassTheta = ", maxMassTheta)
    muT2_Wt_theta = ( math.cos(theta) * mt2_t_min + math.sin(theta) * mt2_W ) / maxMassTheta 
    h_muT2_Wt.Fill( muT2_Wt_theta )
    h_muT2dc.Fill( muT2_Wt_theta + 0.1 )   # Fill histogram with junk for now.
    

# second hghest E electron 
# fix variable names 


test = DC.mt2dc([0.7])  # [0.5, 0.5]
print( test.get_alphaList() )
print( test.get_mT2dc_maxMass([10]) )
print( test.get_mT2dc() )
print( test.get_muT2dc([10]) )
[mySettings, myResults] = test.get_allSettingsAndResults() 
print(mySettings)
print(myResults)

# <<<<<<<<<<<<<<<< Main analysis


####################################
# Draw the histograms and save them.
####################################
c = ROOT.TCanvas()

h_ell1_pt.Draw("E")
c.SaveAs("h_ell1_PT.pdf")

h_ell1_E.Draw("E")
c.SaveAs("h_ell1_E.pdf")

h_bjet1_E.Draw("E")
c.SaveAs("h_bjet1_E.pdf")

h_mT2_W.Draw("E")
c.SaveAs("h_mT2_W.pdf")

h_mT2prime_W.Draw("E")
c.SaveAs("h_mT2prime_W.pdf")


h_mT2_t_11_22.Draw("E")
c.SaveAs("h_mT2_t_11_22.pdf")

h_mT2_t_12_21.Draw("E")
c.SaveAs("h_mT2_t_12_21.pdf")

h_mT2_t_min.Draw("E")
c.SaveAs("h_mT2_t_min.pdf")


h_muT2_Wt.Draw("E")
c.SaveAs("h_muT2_Wt.pdf")

h_muT2dc.Draw("E")
c.SaveAs("h_muT2dc.pdf")

h_mT_ell1.Draw("E")
c.SaveAs("h_mT_ell1.pdf")

h_mT_ell1__forMt2Overlay.Draw("E")
c.SaveAs("h_mT_ell1__forMt2Overlay.pdf")



# save to ROOT output files

h_ell1_pt.Write()
h_ell1_E.Write()
h_bjet1_E.Write()
h_mT2_W.Write()

h_mT2prime_W.Write()
h_mT2_t_11_22.Write()
h_mT2_t_12_21.Write()
h_mT2_t_min.Write()

h_muT2_Wt.Write()
h_muT2dc.Write()

h_mT_ell1.Write()
h_mT_ell1__forMt2Overlay.Write()


f_outputRoot.Close()


