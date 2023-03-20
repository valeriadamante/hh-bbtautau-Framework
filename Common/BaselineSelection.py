import ROOT
import os
from .Utilities import *

print
initialized = False

ana_reco_object_collections = [ "Electron", "Muon", "Tau", "Jet", "FatJet", "boostedTau", "MET", "PuppiMET", "DeepMETResponseTune", "DeepMETResolutionTune"]

def Initialize(loadTF=False, loadHHBtag=False):
    global initialized
    if not initialized:
        headers_dir = os.path.dirname(os.path.abspath(__file__))
        header_path_RootExt = os.path.join(headers_dir, "RootExt.h")
        header_path_GenLepton = os.path.join(headers_dir, "GenLepton.h")
        header_path_Gen = os.path.join(headers_dir, "BaselineGenSelection.h")
        header_path_Reco = os.path.join(headers_dir, "BaselineRecoSelection.h")
        header_path_HHbTag = os.path.join(headers_dir, "HHbTagScores.h")
        ROOT.gInterpreter.Declare(f'#include "{header_path_RootExt}"')
        ROOT.gInterpreter.Declare(f'#include "{header_path_GenLepton}"')
        ROOT.gInterpreter.Declare(f'#include "{header_path_Gen}"')
        ROOT.gInterpreter.Declare(f'#include "{header_path_Reco}"')
        for wpcl in [WorkingPointsTauVSe,WorkingPointsTauVSmu,WorkingPointsTauVSjet,WorkingPointsbTag, WorkingPointsMuonID]:
            ROOT.gInterpreter.Declare(f'{generate_enum_class(wpcl)}')
        if(loadTF):
            import RunKit.includeCMSSWlibs as IncludeLibs
            IncludeLibs.includeLibTool("tensorflow")
        if(loadHHBtag):
            ROOT.gInterpreter.Declare(f'#include "{header_path_HHbTag}"')
            ROOT.gROOT.ProcessLine(f'HHBtagWrapper::Initialize("{os.environ["CMSSW_BASE"]}/src/HHTools/HHbtag/models/", 1)')
        initialized = True

leg_names = [ "Electron", "Muon", "Tau", "boostedTau" ]

channels = [ 'muMu', 'eMu', 'eE', 'muTau', 'eTau', 'tauTau' ] # in order of importance during the channel selection

channelLegs = {
    "eTau": [ "Electron", "Tau" ],
    "muTau": [ "Muon", "Tau" ],
    "tauTau": [ "Tau", "Tau" ],
    "muMu": [ "Muon", "Muon" ],
    "eMu": [ "Electron", "Muon" ],
    "eE": [ "Electron", "Electron" ],
}


def applyMETFlags(df, MET_flags):
    MET_flags_string = ' && '.join(MET_flags)
    return df.Filter(MET_flags_string, "MET filters")

def DefineGenObjects(df, isData=False, isHH=False, Hbb_AK4mass_mpv=125., p4_suffix='nano'):
    if isData:
        df = df.Define("genLeptons", "std::vector<reco_tau::gen_truth::GenLepton>()")
    else:
        df = df.Define("GenPart_daughters", "GetDaughters(GenPart_genPartIdxMother)")
        df = df.Define("genLeptons", """reco_tau::gen_truth::GenLepton::fromNanoAOD(GenPart_pt, GenPart_eta,
                                        GenPart_phi, GenPart_mass, GenPart_genPartIdxMother, GenPart_pdgId,
                                        GenPart_statusFlags, event)""")

    for lep in ["Electron", "Muon", "Tau"]:
        df = df.Define(f"{lep}_genMatchIdx",  f"MatchGenLepton({lep}_p4_{p4_suffix}, genLeptons, 0.2)")
        df = df.Define(f"{lep}_genMatch",  f"GetGenLeptonMatch({lep}_genMatchIdx, genLeptons)")
    if isData:
        return df

    if isHH:
        df = df.Define("genHttCand", """GetGenHTTCandidate(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags,
                                                       GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, false)""")
        df = df.Define("genHbbIdx", """GetGenHBBIndex(event, GenPart_pdgId, GenPart_daughters, GenPart_statusFlags)""")
        df = df.Define("genHbb_isBoosted", "GenPart_pt[genHbbIdx]>550")
    for var in ["GenJet", "GenJetAK8"]:
        df = df.Define(f"{var}_idx", f"CreateIndexes({var}_pt.size())")
        df = df.Define(f"{var}_p4", f"GetP4({var}_pt,{var}_eta,{var}_phi,{var}_mass, {var}_idx)")

    df = df.Define("GenJet_b_PF", "abs(GenJet_partonFlavour)==5")
    df = df.Define("GenJetAK8_b_PF", "abs(GenJetAK8_partonFlavour)==5")
    df = df.Define("GenJet_Hbb",f"FindTwoJetsClosestToMPV({Hbb_AK4mass_mpv}, GenJet_p4, GenJet_b_PF)")
    df = df.Define("GenJetAK8_Hbb", "FindGenJetAK8(GenJetAK8_mass, GenJetAK8_b_PF)")

    return df

def PassGenAcceptance(df):
    df = df.Filter("genHttCand", "GenHttCand present")
    return df.Filter("PassGenAcceptance(*genHttCand)", "GenHttCand Acceptance")

def GenJetSelection(df):
    df = df.Define("GenJet_B1","GenJet_pt > 20 && abs(GenJet_eta) < 2.5 && GenJet_Hbb")
    df = df.Define("GenJetAK8_B1","GenJetAK8_pt > 170 && abs(GenJetAK8_eta) < 2.5 && GenJetAK8_Hbb")
    return df.Filter("GenJet_idx[GenJet_B1].size()==2 || (GenJetAK8_idx[GenJetAK8_B1].size()==1 && genHbb_isBoosted)", "(One)Two b-parton (Fat)jets at least")

def GenJetHttOverlapRemoval(df):
    for var in ["GenJet", "GenJetAK8"]:
        df = df.Define(f"{var}_B2", f"RemoveOverlaps({var}_p4, {var}_B1,{{{{genHttCand->leg_p4[0], genHttCand->leg_p4[1]}},}}, 2, 0.5)" )
    return df.Filter("GenJet_idx[GenJet_B2].size()==2 || (GenJetAK8_idx[GenJetAK8_B2].size()==1 && genHbb_isBoosted)", "No overlap between genJets and genHttCands")

def RequestOnlyResolvedGenJets(df):
    return df.Filter("GenJet_idx[GenJet_B2].size()==2", "Resolved topology")


def SelectRecoP4(df, syst_name='nano'):
    for obj in ana_reco_object_collections:
        df = df.Define(f"{obj}_p4", f"{obj}_p4_{syst_name}")
    return df

def CreateRecoP4(df, suffix='nano'):
    if len(suffix) > 0:
        suffix = "_" + suffix
    if("TrigObj_pt" in df.GetColumnNames()):
        df = df.Define(f"TrigObj_idx", f"CreateIndexes(TrigObj_pt.size())")
        if "TrigObj_mass" not in df.GetColumnNames():
            df = df.Define("TrigObj_mass", "RVecF(TrigObj_pt.size(), 0.f)")
        df = df.Define(f"TrigObj_p4", f"GetP4(TrigObj_pt,TrigObj_eta,TrigObj_phi, TrigObj_mass, TrigObj_idx)")
    for obj in ana_reco_object_collections:
        if "MET" in obj:
            df = df.Define(f"{obj}_p4{suffix}", f"LorentzVectorM({obj}_pt, 0., {obj}_phi, 0.)")
        else:
            df = df.Define(f"{obj}_idx", f"CreateIndexes({obj}_pt.size())")
            df = df.Define(f"{obj}_p4{suffix}", f"GetP4({obj}_pt, {obj}_eta, {obj}_phi, {obj}_mass, {obj}_idx)")
    return df

def DefineMETCuts(met_thr, met_collections):
  cut = ' || '.join([f'{v}_pt > {met_thr}' for v in met_collections ])
  return f"( {cut} )"


def RecoLeptonsSelection(df, apply_filter=True):
    df = df.Define("Electron_B0", f"""
        v_ops::pt(Electron_p4) > 18 && abs(v_ops::eta(Electron_p4)) < 2.3 && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.045
        && (Electron_mvaIso_WP90 || (Electron_mvaNoIso_WP90 && Electron_pfRelIso03_all < 0.5))
    """)

    df = df.Define("Muon_B0", f"""
        v_ops::pt(Muon_p4) > 18 && abs(v_ops::eta(Muon_p4)) < 2.3 && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045
        && ( ((Muon_tightId || Muon_mediumId) && Muon_pfRelIso04_all < 0.5) || (Muon_highPtId && Muon_tkRelIso < 0.5) )
    """)
    df = df.Define("Tau_B0", f"""
        v_ops::pt(Tau_p4) > 15 && abs(v_ops::eta(Tau_p4)) < 2.5 && abs(Tau_dz) < 0.2 && Tau_decayMode != 5 && Tau_decayMode != 6
        && (    (    Tau_idDeepTau2017v2p1VSe >= {WorkingPointsTauVSe.VVLoose.value}
                  && Tau_idDeepTau2017v2p1VSmu >= {WorkingPointsTauVSmu.VLoose.value}
                  && Tau_idDeepTau2017v2p1VSjet >= {WorkingPointsTauVSjet.VVVLoose.value} )
             || (    Tau_idDeepTau2018v2p5VSe >= {WorkingPointsTauVSe.VVLoose.value}
                  && Tau_idDeepTau2018v2p5VSmu >= {WorkingPointsTauVSmu.VLoose.value}
                  && Tau_idDeepTau2018v2p5VSjet >= {WorkingPointsTauVSjet.VVVLoose.value} )
           )
    """)

    df = df.Define("boostedTau_B0", f"""
        v_ops::pt(boostedTau_p4) > 40 && abs(v_ops::eta(boostedTau_p4)) < 2.3 && abs(boostedTau_dz) < 0.2 && boostedTau_decayMode != 5
        && boostedTau_decayMode != 6 && boostedTau_idMVAnewDM2017v2 >= {WorkingPointsBoostedTauVSjet.VVLoose.value}
    """)

    df = df.Define("Electron_B0T", """
        Electron_B0 && (Electron_mvaIso_WP80
                        || (Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all < 0.15))
    """)

    df = df.Define("Muon_B0T", """
        Muon_B0 && ( ((Muon_tightId || Muon_mediumId) && Muon_pfRelIso04_all < 0.15)
                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)

    df = df.Define("Tau_B0T", f"""
        Tau_B0 && (
                    Tau_idDeepTau2017v2p1VSjet >= {WorkingPointsTauVSjet.Medium.value}
                   || Tau_idDeepTau2018v2p5VSjet >= {WorkingPointsTauVSjet.Loose.value} )
    """)

    df = df.Define("boostedTau_B0T", f"""
        boostedTau_B0 && boostedTau_idMVAnewDM2017v2 >= {WorkingPointsBoostedTauVSjet.Medium.value}
    """)

    met_cuts = DefineMETCuts(80, ["MET", "DeepMETResolutionTune", "DeepMETResponseTune", "PuppiMET"])

    ch_filters = []
    for leg1_idx in range(len(leg_names)):
        for leg2_idx in range(max(1, leg1_idx), len(leg_names)):
            leg1, leg2 = leg_names[leg1_idx], leg_names[leg2_idx]
            if leg1 == 'Tau' and leg2 == 'boostedTau': continue
            ch_filter = f"{leg1}{leg2}_B0"
            ch_filters.append(ch_filter)
            if leg1 == leg2:
                ch_filter_def = f"{leg1}_idx[{leg1}_B0].size() > 1 && {leg1}_idx[{leg1}_B0T].size() > 0"
            else:
                ch_filter_def = f"""
                    ({leg1}_idx[{leg1}_B0].size() > 0 && {leg2}_idx[{leg2}_B0T].size() > 0)
                    || ({leg1}_idx[{leg1}_B0T].size() > 0 && {leg2}_idx[{leg2}_B0].size() > 0)
                """
            df = df.Define(ch_filter, ch_filter_def)
        ch_filter = f"{leg1}MET_B0"
        ch_filters.append(ch_filter)
        ch_filter_def = f"{leg1}_idx[{leg1}_B0T].size() > 0 && {met_cuts}"
        df = df.Define(ch_filter, ch_filter_def)
    filter_expr = " || ".join(ch_filters)
    if apply_filter:
        return df.Filter(filter_expr, "Reco leptons requirements")
    else:
        return df, filter_expr

def RecoHttCandidateSelection(df, config):
    df = df.Define("Electron_iso", "Electron_pfRelIso03_all") \
           .Define("Muon_iso", "Muon_pfRelIso04_all") \
           .Define("Tau_iso", "-Tau_rawDeepTau2017v2p1VSjet")

    df = df.Define("Electron_B2_eTau_1", f"Electron_B0 && v_ops::pt(Electron_p4) > 20 && Electron_mvaIso_WP80")
    #df = df.Define("Electron_B2_eTau_1", f"Electron_B0 && v_ops::pt(Electron_p4) > 20 && Electron_mvaNoIso_WP80 && Electron_pfRelIso03_all<0.1")
    df = df.Define("Tau_B2_eTau_2", f"""
        Tau_B0 && v_ops::pt(Tau_p4) > 20
        && (Tau_idDeepTau2017v2p1VSe >= {getattr(WorkingPointsTauVSe, config["deepTauWPs"]["eTau"]["VSe"]).value})
        && (Tau_idDeepTau2017v2p1VSmu >= {getattr(WorkingPointsTauVSmu, config["deepTauWPs"]["eTau"]["VSmu"]).value})
    """)

    df = df.Define("Muon_B2_muTau_1", f"""
        Muon_B0 && v_ops::pt(Muon_p4) > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.15)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)
    df = df.Define("Tau_B2_muTau_2", f"""
        Tau_B0 && v_ops::pt(Tau_p4) > 20
        && (Tau_idDeepTau2017v2p1VSe >= {getattr(WorkingPointsTauVSe, config["deepTauWPs"]["muTau"]["VSe"]).value})
        && (Tau_idDeepTau2017v2p1VSmu >= {getattr(WorkingPointsTauVSmu, config["deepTauWPs"]["muTau"]["VSmu"]).value})
    """)

    df = df.Define("Tau_B2_tauTau_1", f"""
        Tau_B0 && v_ops::pt(Tau_p4) > 20
        && (Tau_idDeepTau2017v2p1VSe >= {getattr(WorkingPointsTauVSe, config["deepTauWPs"]["tauTau"]["VSe"]).value})
        && (Tau_idDeepTau2017v2p1VSmu >= {getattr(WorkingPointsTauVSmu, config["deepTauWPs"]["tauTau"]["VSmu"]).value})
        && (Tau_idDeepTau2017v2p1VSjet >= {getattr(WorkingPointsTauVSjet, config["deepTauWPs"]["tauTau"]["VSjet"]).value})
    """)

    df = df.Define("Tau_B2_tauTau_2", f"""
        Tau_B0 && v_ops::pt(Tau_p4) > 20
        && (Tau_idDeepTau2017v2p1VSe >= {getattr(WorkingPointsTauVSe, config["deepTauWPs"]["tauTau"]["VSe"]).value})
        && (Tau_idDeepTau2017v2p1VSmu >= {getattr(WorkingPointsTauVSmu, config["deepTauWPs"]["tauTau"]["VSmu"]).value})
    """)

    df = df.Define("Muon_B2_muMu_1", f"""
        Muon_B0 && v_ops::pt(Muon_p4) > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.15)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)
    df = df.Define("Muon_B2_muMu_2", f"""
        Muon_B0 && v_ops::pt(Muon_p4) > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.3)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.3) )
    """)

    df = df.Define("Electron_B2_eMu_1", f"""
        Electron_B0 && v_ops::pt(Electron_p4) > 20 && Electron_mvaIso_WP80
    """)
    df = df.Define("Muon_B2_eMu_2", f"""
        Muon_B0 && v_ops::pt(Muon_p4) > 20 && (   (Muon_tightId && Muon_pfRelIso04_all < 0.15)
                                    || (Muon_highPtId && Muon_tkRelIso < 0.15) )
    """)

    df = df.Define("Electron_B2_eE_1", f"""
        Electron_B0 && v_ops::pt(Electron_p4) > 20 &&  Electron_mvaIso_WP80
    """)
    df = df.Define("Electron_B2_eE_2", f"""
                    Electron_B2_eE_1
    """)

    cand_columns = []
    for ch in channels:
        leg1, leg2 = channelLegs[ch]
        cand_column = f"httCands_{ch}"
        df = df.Define(cand_column, f"""
            GetHTTCandidates(Channel::{ch}, 0.5, {leg1}_B2_{ch}_1, {leg1}_p4, {leg1}_iso, {leg1}_charge, {leg1}_genMatchIdx,
                                                 {leg2}_B2_{ch}_2, {leg2}_p4, {leg2}_iso, {leg2}_charge, {leg2}_genMatchIdx)
        """)
        cand_columns.append(cand_column)
    cand_filters = [ f'{c}.size() > 0' for c in cand_columns ]
    df = df.Filter(" || ".join(cand_filters), "Reco Baseline 2")
    cand_list_str = ', '.join([ '&' + c for c in cand_columns])
    return df.Define('httCand', f'GetBestHTTCandidate({{ {cand_list_str} }}, event)')

def ThirdLeptonVeto(df):
    df = df.Define("Electron_vetoSel",
                   f"""v_ops::pt(Electron_p4) > 10 && abs(v_ops::eta(Electron_p4)) < 2.5 && abs(Electron_dz) < 0.2 && abs(Electron_dxy) < 0.045
                      && ( Electron_mvaIso_WP90 == true || ( Electron_mvaNoIso_WP90 && Electron_pfRelIso03_all<0.3) )
                     && (httCand.isLeg(Electron_idx, Leg::e)== false)""")
    df = df.Filter("Electron_idx[Electron_vetoSel].size() == 0", "No extra electrons")
    df = df.Define("Muon_vetoSel",
                   f"""v_ops::pt(Muon_p4) > 10 && abs(v_ops::eta(Muon_p4)) < 2.5 && abs(Muon_dz) < 0.2 && abs(Muon_dxy) < 0.045
                      && ( Muon_mediumId || Muon_tightId ) && Muon_pfRelIso04_all<0.3
                      && (httCand.isLeg(Muon_idx, Leg::mu) == false)""")
    df = df.Filter("Muon_idx[Muon_vetoSel].size() == 0", "No extra muons")
    return df

def RecoJetSelection(df):
    df = df.Define("Jet_bIncl", f"v_ops::pt(Jet_p4)>20 && abs(v_ops::eta(Jet_p4)) < 2.5 && ( Jet_jetId & 2 ) && (Jet_puId>0 || v_ops::pt(Jet_p4)>50)")
    df = df.Define("FatJet_bbIncl", "FatJet_msoftdrop > 30 && abs(v_ops::eta(FatJet_p4)) < 2.5")
    df = df.Define("Jet_bCand", "RemoveOverlaps(Jet_p4, Jet_bIncl,{{httCand.leg_p4[0], httCand.leg_p4[1]},}, 2, 0.5)")
    df = df.Define("FatJet_bbCand", "RemoveOverlaps(FatJet_p4, FatJet_bbIncl, {{httCand.leg_p4[0], httCand.leg_p4[1]},}, 2, 0.5)")
    return df.Filter("Jet_idx[Jet_bCand].size()>=2 || FatJet_idx[FatJet_bbCand].size()>=1", "Reco bjet candidates")


def RequestOnlyResolvedRecoJets(df):
    return df.Filter("Jet_idx[Jet_bCand].size()>=2", "Reco resolved bjet candidates")


def GenRecoJetMatching(df):
    df = df.Define("Jet_genJetIdx_matched", "GenRecoJetMatching(event,Jet_idx, GenJet_idx, Jet_bCand, GenJet_B2, GenJet_p4, Jet_p4 , 0.3)")
    df = df.Define("Jet_genMatched", "Jet_genJetIdx_matched>=0")
    return df.Filter("Jet_genJetIdx_matched[Jet_genMatched].size()>=2", "Two different gen-reco jet matches at least")

def DefineHbbCand(df):
    df = df.Define("Jet_HHBtagScore", "GetHHBtagScore(Jet_bCand, Jet_idx, Jet_p4,Jet_btagDeepFlavB, MET_pt,  MET_phi, httCand, period, event)")
    df = df.Define("HbbCandidate", "GetHbbCandidate(Jet_HHBtagScore, Jet_bCand, Jet_p4, Jet_idx)")
    return df