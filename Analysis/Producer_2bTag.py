import ROOT
import Common.Utilities as Utilities
import sys
import os
import math
from TauIDSFs_modifier import *
from Producer import *

deepJet_weights = ["weight_bTagSF_Loose_Central", "weight_bTagSF_Medium_Central","weight_bTagSF_Tight_Central"]


deepJet_WPs = {
    "Loose": 0.049,
    "Medium":0.2783,
    "Tight":0.71
}
pNet_WPs = {
    "Loose": 0.049,
    "Medium":0.2783,
    "Tight":0.71
}

if __name__ == "__main__":
    import argparse
    import PlotKit.Plotter as Plotter
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--period', required=False, type=str, default = 'Run2_2018')
    parser.add_argument('--version', required=False, type=str, default = 'v2_deepTau_v2p1')
    parser.add_argument('--vars', required=False, type=str, default = 'tau1_pt')
    parser.add_argument('--mass', required=False, type=int, default=500)
    parser.add_argument('--signalScale', required=False, type=int, default=1)
    parser.add_argument('--new-weights', required=False, type=bool, default=True)
    parser.add_argument('--onlyWeightCentral', required=False, type=bool, default=False)
    parser.add_argument('--useTauIDWeight', required=False, type=bool, default=False)
    parser.add_argument('--customName', required=False, type=str, default = '')
    parser.add_argument('--workingPoint', required=False, type=str, default = 'Medium')
    parser.add_argument('--useParticleNet', required=False, type=bool, default = False)
    args = parser.parse_args()
    if args.new_weights:
        print(f"using new weights")
    abs_path = os.environ['CENTRAL_STORAGE']
    anaTuplePath= f"/eos/home-k/kandroso/cms-hh-bbtautau/anaTuples/Run2_2018/{args.version}/"
    page_cfg = "config/plot/cms_stacked.yaml"
    page_cfg_custom = "config/plot/2018.yaml"
    hist_cfg = "config/plot/histograms.yaml"
    inputs_cfg = "config/plot/inputs.yaml"
    with open(inputs_cfg, 'r') as f:
        inputs_cfg_dict = yaml.safe_load(f)

    plotter = Plotter.Plotter(page_cfg=page_cfg, page_cfg_custom=page_cfg_custom, hist_cfg=hist_cfg, inputs_cfg=inputs_cfg_dict)

    dataframes = {}
    for sample in files.keys():
        rootFiles = [anaTuplePath+f + ".root" for f in files[sample]]
        df = ROOT.RDataFrame("Events", rootFiles)
        anaskimmer = AnaSkimmer(df, args.version.split('_')[-1])
        if(sample in signals):
            for input in inputs_cfg_dict:
                name = input['name']
                if(name == sample):
                    input['title']+= f"mass {args.mass} GeV (\sigma = {args.signalScale} pb)"
                    input['scale'] = args.signalScale
            anaskimmer.df = anaskimmer.df.Filter(f"X_mass=={args.mass}")
        anaskimmer.skimAnatuple()
        if args.useParticleNet:
            dataframes[sample] = anaskimmer.df.Filter(f"""b1_particleNetAK4_B >= {pNet_WPs[args.workingPoint]} && b2_particleNetAK4_B >= {pNet_WPs[args.workingPoint]}""")
        else:
            dataframes[sample] = anaskimmer.df.Filter(f"""b1_btagDeepFlavB >= {deepJet_WPs[args.workingPoint]} && b2_btagDeepFlavB >= {deepJet_WPs[args.workingPoint]}""")
            weights_to_apply.append(f"weight_bTagSF_{args.workingPoint}_Central")
    defineWeights(dataframes, args.new_weights,args.version.split('_')[-1], args.onlyWeightCentral, args.useTauIDWeight)

    all_histograms = {}
    vars = args.vars.split(',')
    all_sums = createSums(dataframes)

    for var in vars:
        hists = createHistograms(dataframes, var, plotter)
        all_histograms[var] = hists

    hists_to_plot = {}
    all_histograms=GetValues(all_histograms)
    all_sums=GetValues(all_sums)

    for var in vars:
        hists_to_plot[var] = {}
        for sample in all_histograms[var].keys():
            for region in regions:
                hists_to_plot[var][sample] = all_histograms[var][sample]['region_A']
        hists_to_plot[var]['QCD'] = Estimate_QCD(all_histograms[var], all_sums)
        custom1= {'cat_text':'2bTag'}
        fullName = f"2bTag_{var}_XMass{args.mass}_deepTau{args.version.split('_')[-1]}"
        if args.customName != "":
            fullName += f"_{args.customName}"
        if args.new_weights:
            fullName += f"_newWeights"
        plotter.plot(var, hists_to_plot[var], f"output/plots/{fullName}.pdf", custom=custom1)
        for sample in  hists_to_plot[var].keys():
            print(f"{sample}, {hists_to_plot[var][sample].Integral()}")