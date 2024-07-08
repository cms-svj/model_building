import awkward as ak
import pickle
import numpy as np
import hist
import matplotlib as mpl
from coffea.nanoevents import NanoEventsFactory, DelphesSchema

def ET(vec):
    return np.sqrt(vec.px**2+vec.py**2+vec.mass**2)

def get_values(var,Events):
    return ak.flatten(Events[var],axis=None)

def fill_hist(var,nbins,bmin,bmax,label,Events):
    h = (
        hist.Hist.new
        .Reg(nbins, bmin, bmax, label=label)
        .Double()
    )
    h.fill(get_values(var,Events),weight=0.5)
    return h

def normalize_angle(angle):
    angle = np.mod(angle, 2 * np.pi)
    angle = np.where(angle >= np.pi, angle - 2 * np.pi, angle)
    return angle

def histogram(filename, helper):                    
    DelphesSchema.mixins["FatJet"] = "Jet"
    Events = NanoEventsFactory.from_root(
                file=filename,
                    treepath="Delphes",
                        schemaclass=DelphesSchema,
                        ).events()
    
    # require two jets
    events = Events
    mask = ak.num(events.FatJet)>=2
    events = events[mask]

    # Dijet
    events["Dijet"] = events.FatJet[:,0]+events.FatJet[:,1]

    # transverse mass calculation
    E1 = ET(events.Dijet)
    E2 = events.MissingET.MET
    MTsq = (E1+E2)**2-(events.Dijet.px+events.MissingET.px)**2-(events.Dijet.py+events.MissingET.py)**2
    events["MT"] = np.sqrt(MTsq, where=MTsq>=0)
            
    # 4-vectors for dijet
    events["Dijet_pt"] = events.Dijet.pt
    events["Dijet_eta"] = events.Dijet.eta
    events["Dijet_phi"] = events.Dijet.phi
    events["Dijet_mass"] = events.Dijet.mass

    events["MET"] = events.MissingET.MET

    ## For plotting individually for jet1 and jet2

    events["Jet1"] = events.FatJet[:,0]
    events["Jet2"] = events.FatJet[:,1]

    # 4-vectors for jet1 and jet2
    events["Jet1_pt"] = events["Jet1"].pt
    events["Jet2_pt"] = events["Jet2"].pt
    
    events["Jet1_eta"] = events["Jet1"].eta
    events["Jet2_eta"] = events["Jet2"].eta
    
    events["Jet1_phi"] = events["Jet1"].phi
    events["Jet2_phi"] = events["Jet2"].phi
    
    events["Jet1_mass"] = events["Jet1"].mass
    events["Jet2_mass"] = events["Jet2"].mass
    
    events["DeltaEta"] = np.abs(events["Jet2_eta"] - events["Jet1_eta"])
    events["DeltaPhi"] = np.abs(normalize_angle(events["Jet1_phi"] - events["Jet2_phi"]))
    
    events["DeltaPhi_MET_Jet1"] = np.abs(normalize_angle(events.MissingET.phi - events["Jet1_phi"]))
    events["DeltaPhi_MET_Jet2"] = np.abs(normalize_angle(events.MissingET.phi - events["Jet2_phi"]))

    # Stable inv frac
    dark_hadron_ids = helper.darkHadronIDs
    stable_particle_ids = helper.stableIDs
    pid = events.GenParticle["PID"]
    d1 = events.GenParticle["D1"]

    # Boolean array of whether a particle is dark
    
    is_dark = ak.zeros_like(pid)
    for dhid in dark_hadron_ids:
        is_dark = is_dark | (np.abs(pid)==dhid)
    is_dark = is_dark==1

    # PIDs of dark daughter
    dark_daughter = pid[d1[is_dark]] 
    is_dark_daughter = ak.zeros_like(dark_daughter)

    for dsid in stable_particle_ids:
        is_dark_daughter = is_dark_daughter | (np.abs(dark_daughter)==dsid)

    stability = [None]*len(events)
    for i in range(len(pid)): 
        if ak.sum(is_dark[i])!=0:
            stability[i] = ak.sum(is_dark_daughter[i])/ak.sum(is_dark[i])

    # Add the invisible fraction to the events
    events["stable_invisible_fraction"] = stability

    # store modified events array
    Events = events
    
    
    # Histogram

    # Creating hist objects

    hist_dict = {
        "Jet_mt": fill_hist("MT",25,0,1500,r"$m_{\text{T}}$ [GeV]",Events),
        "Dijet_pt": fill_hist("Dijet_pt",50,0,1000,r"$p_{\text{T}}$ [GeV]",Events),
        "Jet1_pt": fill_hist("Jet1_pt",50,0,1000,r"$Jet_1\_p_{\text{T}}$ [GeV]",Events),
        "Jet2_pt": fill_hist("Jet2_pt",50,0,1000,r"$Jet_2\_p_{\text{T}}$ [GeV]",Events),
        "MET": fill_hist("MET",50,0,1000,r"$p_{\text{T},miss}$ [GeV]",Events),
        "Dijet_eta": fill_hist("Dijet_eta",50,-10,10,r"$Jet_{\eta}$ [GeV]",Events),
        "Jet1_eta": fill_hist("Jet1_eta",50,-6,6,r"$\eta_{Jet1}$",Events),
        "Jet2_eta": fill_hist("Jet2_eta",50,-6,6,r"$\eta_{Jet2}$",Events),
        "Dijet_phi": fill_hist("Dijet_phi",25,-4,4,r"$Jet_{\phi}$",Events),
        "Jet1_phi": fill_hist("Jet1_phi",25,-4,4,r"$\phi_{Jet1}$",Events),
        "Jet2_phi": fill_hist("Jet2_phi",25,-4,4,r"$\phi_{Jet2}$",Events),
        "Dijet_mass": fill_hist("Dijet_mass",50,0,2300,r"$Jet_{mass}$ [GeV]",Events),
        "Jet1_mass": fill_hist("Jet1_mass",50,0,250,r"$Jet1_{mass}$ [GeV]",Events),
        "Jet2_mass": fill_hist("Jet2_mass",50,0,250,r"$Jet2_{mass}$ [GeV]",Events),
        "DeltaEta": fill_hist("DeltaEta",35,0,8.0,r"$\Delta\eta$",Events),
        "DeltaPhi": fill_hist("DeltaPhi",20,0,3.5,r"$\Delta\Phi$",Events),
        "DeltaPhi_MET_Jet1": fill_hist("DeltaPhi_MET_Jet1",25,0,3.5,r"$\Delta\phi\_P_{t,miss}\_Jet_1$",Events),
        "DeltaPhi_MET_Jet2": fill_hist("DeltaPhi_MET_Jet2",25,0,3.5,r"$\Delta\phi\_P_{t,miss}\_Jet_2$",Events),
        "stable_invisible_fraction": fill_hist("stable_invisible_fraction",25,0,1,r"stable_invisible_fraction",Events)
    }
    
    # Saving the histograms
    with open("Hists.pkl", "wb") as out:        
        pickle.dump(hist_dict, out)
    
