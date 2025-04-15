#LC TTREE
from tqdm import tqdm
import numpy as np
import ROOT
#load data
chain = ROOT.TChain(f"DecayTree");
chain.Add("/sphenix/tg/tg01/hf/cdean/QM25_productions_TrackVertexMatch/pKpi_reco/*.root")
# Create a TTreeReader to read the TTree
reader = ROOT.TTreeReader(chain)

#for dalitz plot
p_PE = ROOT.TTreeReaderValue('Float_t')(reader, "track_3_pE")
p_PX = ROOT.TTreeReaderValue('Float_t')(reader, "track_3_px")
p_PY = ROOT.TTreeReaderValue('Float_t')(reader, "track_3_py")
p_PZ = ROOT.TTreeReaderValue('Float_t')(reader, "track_3_pz")
K_PE = ROOT.TTreeReaderValue('Float_t')(reader, "track_2_pE")
K_PX = ROOT.TTreeReaderValue('Float_t')(reader, "track_2_px")
K_PY = ROOT.TTreeReaderValue('Float_t')(reader, "track_2_py")
K_PZ = ROOT.TTreeReaderValue('Float_t')(reader, "track_2_pz")
pi_PE = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_pE")
pi_PX = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_px")
pi_PY = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_py")
pi_PZ = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_pz")
track_3_p = ROOT.TTreeReaderValue('Float_t')(reader, "track_3_p")
track_3_dEdx = ROOT.TTreeReaderValue('Float_t')(reader, "track_3_dEdx")

Lambda_cplus_mass = ROOT.TTreeReaderValue('Float_t')(reader,"Lambda_cplus_mass")
Lambda_cplus_DIRA = ROOT.TTreeReaderValue('Float_t')(reader,"Lambda_cplus_DIRA")
Lambda_cplus_decayLength = ROOT.TTreeReaderValue('Float_t')(reader,"Lambda_cplus_decayLength")
Lambda_cplus_pT = ROOT.TTreeReaderValue('Float_t')(reader,"Lambda_cplus_pT")

track_1_pT = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_pT")
track_2_pT = ROOT.TTreeReaderValue('Float_t')(reader, "track_2_pT")

track_1_INTT_nHits = ROOT.TTreeReaderValue('UInt_t')(reader, "track_1_INTT_nHits")
track_2_INTT_nHits = ROOT.TTreeReaderValue('UInt_t')(reader, "track_2_INTT_nHits")
track_3_INTT_nHits = ROOT.TTreeReaderValue('UInt_t')(reader, "track_3_INTT_nHits")

track_1_IP = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_IP")  
track_2_IP = ROOT.TTreeReaderValue('Float_t')(reader, "track_2_IP")  
track_3_IP = ROOT.TTreeReaderValue('Float_t')(reader, "track_3_IP")  

track_1_track_2_DCA = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_track_2_DCA") 
track_1_track_3_DCA = ROOT.TTreeReaderValue('Float_t')(reader, "track_1_track_3_DCA") 
track_2_track_3_DCA = ROOT.TTreeReaderValue('Float_t')(reader, "track_2_track_3_DCA") 

nPrimaryVertices = ROOT.TTreeReaderValue('Int_t')(reader, "nPrimaryVertices") 


h_M = ROOT.TH1D("h_M","Lc_M",140,2,2.7)
h_Lc_pK = ROOT.TH1D("h_Lc_pK","h_Lc_pK",2001,1,3)
h_Lc_pip = ROOT.TH1D("h_Lc_pip","h_Lc_pip",1501,1,2.5)
h_Lc_piK = ROOT.TH1D("h_Lc_piK","h_Lc_piK",2001,0,2)

#main event loop
for _ in tqdm(range(chain.GetEntries()), desc="Processing Lc TTree entries"):
    reader.Next()
    
    if not (track_3_p.__deref__() >=1.5 or (track_3_p.__deref__() >= 1.3 and track_3_dEdx.__deref__() > 450) or (track_3_p.__deref__() >= 1.0 and track_3_dEdx.__deref__() > 550) or (track_3_p.__deref__() >= 0.5 and track_3_dEdx.__deref__() > 800) or (track_3_p.__deref__() >= 0.7 and track_3_dEdx.__deref__() > 600) or (track_3_p.__deref__() < 0.4 and track_3_dEdx.__deref__() > 2000) or (track_3_p.__deref__() > 0.4 and track_3_dEdx.__deref__() > 1300)):
        continue
        
    if not (Lambda_cplus_DIRA.__deref__() > 0.90):
        continue
    
    if not (0.003 <= Lambda_cplus_decayLength.__deref__() < 0.2):
        continue
        
    if not (track_1_pT.__deref__() > 0.4 and track_2_pT.__deref__() > 0.2):
        continue
        
    if not (track_1_INTT_nHits.__deref__() >= 2 and track_2_INTT_nHits.__deref__() >= 2 and track_3_INTT_nHits.__deref__() >=2):
        continue
    
    if not (track_1_IP.__deref__() >= 0.005 and track_2_IP.__deref__() >= 0.005 and track_3_IP.__deref__() >=0.005):
        continue
        
    if not (track_1_track_2_DCA.__deref__() < 0.01 and track_1_track_3_DCA.__deref__() < 0.01 and track_2_track_3_DCA.__deref__() < 0.01):
        continue
        
    #if Lambda_cplus_pT.__deref__() > 2:
    #    continue
    
    Lc_pK = (np.sqrt((p_PE.__deref__() + K_PE.__deref__())**2 - (p_PX.__deref__() + K_PX.__deref__())**2 - (p_PY.__deref__() + K_PY.__deref__())**2 - (p_PZ.__deref__() + K_PZ.__deref__())**2))
    Lc_piK = (np.sqrt((pi_PE.__deref__() + K_PE.__deref__())**2 - (pi_PX.__deref__() + K_PX.__deref__())**2 - (pi_PY.__deref__() + K_PY.__deref__())**2 - (pi_PZ.__deref__() + K_PZ.__deref__())**2))
    Lc_pip =(np.sqrt((pi_PE.__deref__() + p_PE.__deref__())**2 - (pi_PX.__deref__() + p_PX.__deref__())**2 - (pi_PY.__deref__() + p_PY.__deref__())**2 - (pi_PZ.__deref__() + p_PZ.__deref__())**2))
    h_Lc_pK.Fill(Lc_pK)
    h_Lc_piK.Fill(Lc_piK)
    h_Lc_pip.Fill(Lc_pip)
    h_M.Fill(Lambda_cplus_mass.__deref__())
    

c_M = ROOT.TCanvas("c_M","c_M",4*400,1*400)
c_M.Divide(4,1)
c_M.cd(1)
h_M.Draw()
c_M.cd(2)
h_Lc_pK.Draw()
c_M.cd(3)
h_Lc_piK.Draw()
c_M.cd(4)
h_Lc_pip.Draw()
#c_M.cd(2)
#chain.Draw("Lambda_cplus_DIRA")
c_M.SaveAs("sig.png")



