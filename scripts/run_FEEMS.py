#!/usr/bin/env python

# Python script to run FEEMS/FEEMSmix from a terminal versus interactively
# Please be mindful of default options when fitting and plotting 
# If a full fit has already been run, then you can comment out specific 
# chunks and feed in appropriate values where needed
# April 2025

# importing libraries
# base
import numpy as np
import pandas as pd
from importlib import resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
# import statsmodels.api as sm
import pickle

# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as clr

# feems
from feems.utils import prepare_graph_inputs, cov_to_dist
from feems.objective import comp_mats
from feems.viz import draw_FEEMSmix_surface, plot_FEEMSmix_summary
from feems import SpatialGraph, Viz, Objective
from feems.cross_validation import run_cv, run_cv_joint

import argparse

parser = argparse.ArgumentParser('feems_script')
parser.add_argument('-d', '--dsample', help = 'downsample number, or 0 if not downsampling', type = int, default = 0)
parser.add_argument('-c', '--cross_validate', help = 'do cross validation?', type = int, default = 0)
parser.add_argument('-l', '--lamb', help = 'lambda smoothing value', type = int, default = 10)
parser.add_argument('-n', '--node_specific_variance', help = 'do FEEMS 2.0?', type = int, default = 0)
parser.add_argument('-b', '--bootstrap', help = 'bootstrap?', type = int, default = 0)
parser.add_argument('-r', '--rep', help = 'rep_number', type = int, default = 1)
parser.add_argument('-s', '--shuffle', help = 'shuffle coordinates?', type = int, default = 0)

args = parser.parse_args()

print('Downsampling to {}'.format(args.dsample))
print('Cross validation? {}'.format(args.cross_validate))
print('Lambda = {}'.format(args.lamb))
print('FEEMS 2.0? {}'.format(args.node_specific_variance))
print('bootstrap? {}'.format(args.bootstrap))
print('Rep {}'.format(args.rep))
print('Shuffling? {}'.format(args.shuffle))

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

#---------- CHANGE VARIABLES BELOW ----------
head_dir = "/group/jrigrp10/maize-linguistics"
dataset = 'SeeD'
path_to_plink = "{}/data/feems/all_linguistics_samples_{}".format(head_dir, dataset)
path_to_sample_coords = "{}/data/feems/all_samples_coordinates_{}.txt".format(head_dir, dataset)
K = 3 # number of long-range edges to fit
path_to_output_dir = "{}/results/feems".format(head_dir)

# optional, if chosen change outer=outer in prepare_graph_output
# outer = np.loadtxt("/path/to/outer/boundary/file")

#---------- INPUT FILES ----------
print('\nReading in input data...')
(bim, fam, G) = read_plink(path_to_plink)
coord = np.loadtxt(path_to_sample_coords)
outer = np.loadtxt("{}/data/feems/outer_coords_mesoamerica_manual.txt".format(head_dir)) 

print('Finish loading plink files')

## shuffle coordinates
if args.shuffle:
        print('shuffling coordinates')
        np.random.shuffle(coord)

# ## import pruned SNPs
pruned = pd.read_table("{}/data/feems/all_linguistics_samples_{}.prune.in".format(head_dir, dataset), names = ['keep'])
bim_filtered = bim[bim['snp'].isin(pruned['keep'])]
G_filtered = np.array(G)[bim_filtered['i'].values]
print('Finished pruning')


# imputing any missing genotypes
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
# genotypes = imp.fit_transform((np.array(G)).T)
genotypes = imp.fit_transform((np.array(G_filtered)).T)
print("n_samples={}, n_snps={}\n".format(genotypes.shape[0], genotypes.shape[1]))

# discrete global grid (DGG), could supply custom triangular grid
# data_path = str(resources.files('feems') / 'data')
# grid_path = "{}/grid_250.shp".format(data_path)
grid_path = "{}/data/feems/global_grid_res9.shp".format(head_dir)

outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=True, 
                                             buffer=0,
                                             outer=outer)

#---------- GRAPH SETUP ----------
print('\nSetting up graph...')
dsample = args.dsample
bootstrap = args.bootstrap
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True, downsample = dsample, bootstrap = bootstrap)

print('Finished setting up graph')

# change projection here
projection = ccrs.PlateCarree()

#---------- CROSS VALIDATION ----------
if args.cross_validate:

        print('\nRunning cross-validation scheme...')
        ## this chunk only needs to be run ONCE per data set
        # choosing a discrete log-space grid between 0.01 and 100 
        lamb_grid = np.geomspace(0.01, 100., 10, endpoint=True)[::-1]
        lamb_q_grid = np.geomspace(0.01, 100., 5, endpoint=True)[::-1]

        # using only 5-fold here for faster runtime 
        # but recommended is leave-one-out (default: n_folds = None)
        cv_err = run_cv_joint(sp_graph, lamb_grid, lamb_q_grid, n_folds=5, factr=1e10)
        mean_cv_err = np.mean(cv_err, axis=0)

        lamb_q_cv = lamb_q_grid[np.where(mean_cv_err == np.min(mean_cv_err))[0][0]]
        lamb_cv = lamb_grid[np.where(mean_cv_err == np.min(mean_cv_err))[1][0]]
        print(r"\nlambda_CV values: ({}, {})".format(lamb_cv, lamb_q_cv))

#---------- BASELINE FEEMS FIT ----------
print('\nFitting baseline FEEMS...')
lamb_cv = args.lamb
if lamb_cv == 12:
        lamb_q_cv = 0.01
elif lamb_cv == 2:
        lamb_q_cv = 100
else:
        lamb_q_cv = 10

if args.node_specific_variance:
        sp_graph.fit(lamb = lamb_cv, lamb_q = lamb_q_cv, optimize_q='n-dim')
else:
        sp_graph.fit(lamb = lamb_cv, optimize_q = None)

# visualizing the baseline migration surface
projection = ccrs.EquidistantConic(central_longitude = -118.3, central_latitude = 12.74)
new_colors = [
  "#994000",
  "#CC5800",
  "#FF8F33",
  "#FFAD66",
  "#FFCA99",
  "#FFE6CC",
  "#BEBEBE",
  "#CCFDFF",
  "#99F8FF",
  "#66F0FF",
  "#33E4FF",
  "#00AACC",
  "#007A99",
]

fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
ax.add_feature(cfeature.BORDERS, edgecolor='gray', linewidth=0.3)
ax.add_feature(cfeature.LAND, facecolor='#f7f7f7')

v = Viz(ax, sp_graph, 
        projection=projection, 
        edge_width=.6, edge_alpha=1, edge_zorder=100, 
        sample_pt_size=20, sample_pt_color="black", sample_pt_alpha = 0.10, sample_pt_zorder = 1,
        obs_node_size=7.5, obs_node_alpha = 0.75, obs_node_zorder = 1,  
        cbar_font_size=8,
        cbar_ticklabelsize=8,
        cbar_width="15%",
        cbar_height="4%",)

setattr(v, 'eems_colors', new_colors)
setattr(v, 'edge_cmap', clr.LinearSegmentedColormap.from_list(
    "eems_colors", v.eems_colors, N=256
))
getattr(v, 'edge_cmap')

v.draw_map()
# v.draw_obs_nodes(use_ids=False) 
# v.draw_samples()

v.draw_edges(use_weights=True)
v.draw_edge_colorbar()

fig.savefig(path_to_output_dir + "/baselineFEEMS-lambda_{}-dsample_{}-nsv_{}-cv_{}-boot_{}-rep_{}-shuffle_{}.png".format(args.lamb, args.dsample, 
        args.node_specific_variance, args.cross_validate, args.bootstrap, args.rep, args.shuffle), bbox_inches = 'tight')

print('Finish printing baseline graph')

#---------- FINDING OUTLIERS ---------- 
# outliers_df = sp_graph.extract_outliers(fraction_of_pairs=0.01)

# #---------- FEEMSMIX FIT ----------
# print('\nFitting FEEMSmix with K = {} LREs...'.format(K))
# # see docs/notebooks/further-exploration.ipynb for other modes in which to run this fit
# seq_results = sp_graph.sequential_fit(
#     outliers_df=outliers_df, 
#     lamb=lamb_cv, lamb_q=lamb_q_cv, optimize_q='n-dim', 
#     nedges=K, nedges_to_same_deme=2, top=10,
#     search_area='all',
#     fraction_of_pairs=0.01
# )

# # storing the FEEMSmix output in a pickle for easy access
# filehandler = open(path_to_output_dir+"/seq_results.pkl", 'wb')
# pickle.dump(seq_results, filehandler)

# # visualizing the LREs over the baseline fit
# fig = plt.figure(dpi=300)
# ax = fig.add_subplot(1, 1, 1, projection=color_projection)  
# v = Viz(ax, sp_graph, projection=color_projection, edge_width=.5, 
#         edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
#         obs_node_size=7.5, sample_pt_color="black", 
#         cbar_font_size=10)
# v.draw_map(); v.draw_edges(use_weights=True); v.draw_edge_colorbar()
# v.draw_LREs(seq_results)
# v.draw_c_colorbar()
# fig.savefig(path_to_output_dir+"/LREs.pdf")

# # final summary of FEEMSmix fit
# plot_FEEMSmix_summary(seq_results, sequential=True)
# fig = plt.gcf()
# fig.savefig(path_to_output_dir+"/FEEMSmix_summary.pdf")

