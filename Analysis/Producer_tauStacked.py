import ROOT
import Common.Utilities as Utilities
import sys
import os
import math
from TauIDSFs_modifier import *
from Producer import *

stack_categories = ["inclusive", "2tau", "1tau1notau", "2notau","data", "GluGluToRadion"]
#stack_categories = ["2tau", "1tau1notau", "2notau","data"]

def createSingleHistogram(var, plotter):
    x_bins = plotter.hist_cfg[var]['x_bins']
    if type(plotter.hist_cfg[var]['x_bins'])==list:
        x_bins_vec = Utilities.ListToVector(x_bins, "double")
        return ROOT.TH1D("", "", x_bins_vec.size()-1, x_bins_vec.data())
    else:
        n_bins, bin_range = x_bins.split('|')
        start,stop = bin_range.split(':')
        return ROOT.TH1D("", "",int(n_bins), float(start), float(stop))



if __name__ == "__main__":
    import argparse
    import PlotKit.Plotter as Plotter
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument('--period', required=False, type=str, default = 'Run2_2018')
    parser.add_argument('--version', required=False, type=str, default = 'v2_deepTau_v2p1')
    parser.add_argument('--vars', required=False, type=str, default = 'tau1_pt')
    parser.add_argument('--mass', required=False, type=int, default=500)
    parser.add_argument('--signalScale', required=False, type=int, default=5)
    parser.add_argument('--new-weights', required=False, type=bool, default=False)
    parser.add_argument('--onlyWeightCentral', required=False, type=bool, default=False)
    parser.add_argument('--useTauIDWeight', required=False, type=bool, default=False)
    parser.add_argument('--customName', required=False, type=str, default = '')
    args = parser.parse_args()
    if args.new_weights:
        print(f"using new weights")
    abs_path = os.environ['CENTRAL_STORAGE']
    anaTuplePath= f"/eos/home-k/kandroso/cms-hh-bbtautau/anaTuples/Run2_2018/{args.version}/"
    page_cfg = "config/plot/cms_stacked.yaml"
    page_cfg_custom = "config/plot/2018.yaml"
    hist_cfg = "config/plot/histograms.yaml"
    inputs_cfg = "config/plot/inputs_tauStacked.yaml"
    with open(inputs_cfg, 'r') as f:
        inputs_cfg_dict = yaml.safe_load(f)

    plotter = Plotter.Plotter(page_cfg=page_cfg, page_cfg_custom=page_cfg_custom, hist_cfg=hist_cfg, inputs_cfg=inputs_cfg_dict)

    dataframes = {}
    for cat in stack_categories:
        dataframes[cat] = {}
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
        if sample=='data':
            dataframes["inclusive"][sample] = anaskimmer.df
            dataframes["data"][sample] = anaskimmer.df
        elif sample in signals:
            dataframes["GluGluToRadion"][sample] = anaskimmer.df
        else:
            dataframes["inclusive"][sample] = anaskimmer.df
            dataframes["2tau"][sample] = anaskimmer.df.Filter('tau1_gen_kind == 5 && tau2_gen_kind == 5')
            dataframes["1tau1notau"][sample] = anaskimmer.df.Filter('((tau1_gen_kind == 5 && tau2_gen_kind != 5) || (tau1_gen_kind != 5 && tau2_gen_kind == 5))')
            dataframes["2notau"][sample] = anaskimmer.df.Filter('tau1_gen_kind != 5 && tau2_gen_kind != 5')

    all_histograms = {}
    all_sums = {}
    vars = args.vars.split(',')
    for var in vars:
        all_histograms[var] = {}


    for cat in dataframes.keys():
        defineWeights(dataframes[cat], args.new_weights,args.version.split('_')[-1], args.onlyWeightCentral, args.useTauIDWeight)
        all_sums[cat] = createSums(dataframes[cat])
        for var in vars:
            hists = createHistograms(dataframes[cat], var, plotter)
            all_histograms[var][cat] = hists

    hists_to_plot = {}
    all_histograms=GetValues(all_histograms)
    all_sums=GetValues(all_sums)
    print(all_sums)
    for cat  in all_sums.keys():
        total_yield = 0
        for sample in all_sums[cat].keys():
            total_yield+=all_sums[cat][sample]["region_A"]
        print(cat, total_yield)


    # 1. estimate QCD

    for var in vars:
        hists_to_plot[var] = {}
        for cat in all_histograms[var].keys():
            if cat=="data" or cat == 'inclusive' or cat == 'GluGluToRadion': continue
            hists_to_plot[var][cat] = createSingleHistogram(var, plotter)
            if cat == '2notau':
                QCD_dict = all_histograms[var]["inclusive"]
                QCD_sum_dict = all_sums["inclusive"]
                hists_to_plot[var][cat].Add(Estimate_QCD(QCD_dict, QCD_sum_dict))
            for sample in all_histograms[var][cat].keys():
                if sample=='data' : continue
                hists_to_plot[var][cat].Add(all_histograms[var][cat][sample]['region_A'])
            custom1= {'cat_text':'inclusive'}
            fullName = f"tauStacked_{var}_XMass{args.mass}_deepTau{args.version.split('_')[-1]}"
            if args.customName != "":
                fullName += f"_{args.customName}"
            if args.new_weights:
                fullName += f"_newWeights"
        hists_to_plot[var]["data"] = all_histograms[var]["data"]["data"]['region_A']
        hists_to_plot[var]["GluGluToRadion"] = all_histograms[var]["GluGluToRadion"]["GluGluToRadion"]['region_A']
        #print(hists_to_plot[var].keys())
        plotter.plot(var, hists_to_plot[var], f"output/plots/{fullName}.pdf", custom=custom1)
        for cat in hists_to_plot[var].keys():
            print(f"{var}, {cat}, {hists_to_plot[var][cat].Integral()}")
