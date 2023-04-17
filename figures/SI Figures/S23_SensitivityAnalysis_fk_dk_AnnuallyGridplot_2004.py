from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from Libs.Standard import StandardLPJG
from A02_MCWD_Dataset_Analysis.Pylibs.MCWDevaluation import MCWDFile as MF
from A03_Hydraulics_Implementation.AutoAnalysis.Source.SetupFiles import SetupFile
from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis
from A03_Hydraulics_Implementation.AutoAnalysis.PaperFigures_N2.Revisions2022.BiomassAmazonConverter import BiomassAmazonAnomalyConverter
from sys import platform

import matplotlib.pyplot as plt
import matplotlib.colors

import pandas as pd
import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from scipy.stats import gaussian_kde as kde
from sklearn.metrics import mean_squared_error
import scipy.stats
import warnings
from matplotlib.colors import Normalize
from matplotlib import cm

import string
lowerletters = string.ascii_lowercase[:26]


filename = "/Users/pp/data/Simulations/A03_Hydraulics/AB/2023_revision/AnnuallyOutSens.nc"
config ="/Users/pp/data/Simulations/A03_Hydraulics/AB/2023_revision/Values.tsv"

file = AmazonBasinStratified(filePath=filename, configPath=config, configheader=0)

vars = file.stratifiedOutput.PftNames()


#vars = ['aetplant']


fks = [0.72, 0.8, 0.88]
numcols = len(fks)
dks = [7.2, 8.0, 8.8]
numrows = len(dks)

vars.remove('cmass_litter')
vars.remove('lambda')
vars.remove('omega')
vars.remove('litter_root_inc')
vars.remove('litter_leaf_inc')

data2004_base = file.GetDataSlicesRange('cmasstotal', 2004, 2004, 0)
dataFrameCond = file.inputDataFrame
xnorm = (dataFrameCond['fk'] == 0.8) & (dataFrameCond['dk'] == 8.0)
slice = data2004_base[xnorm]
biomass_indexes = slice > 10.0


data2004_base = file.GetDataSlicesRange('cav_xylem', 2004, 2004, 0)
dataFrameCond = file.inputDataFrame
xnorm =(dataFrameCond['fk'] == 0.8) & (dataFrameCond['dk'] == 8.0)
slice = data2004_base[xnorm]
biomass_indexes = slice < 0.05


for var in vars:
        print(var)
        file = AmazonBasinStratified(filePath=filename, configPath=config, configheader=0)
        data2004 = file.GetDataSlicesRange(var, 2004, 2004, -1)


        data_cmass = data2004

        dataFrameCond = file.inputDataFrame

        img = np.zeros((6, 3, 3))


        mins = 10000000;
        maxs = -10000000;
        for gc in range(0, 6):
            for fki in range(0, len(fks)):
                for dki in range(0, len(dks)):
                    xnorm = (dataFrameCond['fk'] == 0.8) & (dataFrameCond['dk'] == 8.0)
                    xd = (dataFrameCond['fk'] == fks[fki]) & (dataFrameCond['dk'] == dks[dki])

                    indices = dataFrameCond[xd].index
                    indices_norm = dataFrameCond[xnorm].index

                    data_cmass_slice = data_cmass[indices]
                    data_cmass_slice_norm = data_cmass[indices_norm]

                    mean_diff = np.mean(data_cmass_slice[biomass_indexes[:, gc], gc] - data_cmass_slice_norm[biomass_indexes[:, gc], gc])

                    mins = min(mins, mean_diff)
                    maxs = max(maxs, mean_diff)



        max_abs = max(-mins, maxs)

        norm = mpl.colors.Normalize(-max_abs, max_abs)
        colors = [[norm(-max_abs), "tab:red"],
                  [norm(0.0), "white"],
                  [norm(max_abs), "tab:blue"]]
        cmap = mpl.colors.LinearSegmentedColormap.from_list("", colors)

        fig = plt.figure(figsize=(15, 9))
        index = 1
        for gc in range(0, 6):


            ax = fig.add_subplot(2, 3, index)

            img = np.zeros((3, 3))



            for fki in range(0, len(fks)):
                for dki in range(0, len(dks)):

                    xnorm = (dataFrameCond['fk'] == 0.8) & (dataFrameCond['dk'] == 8.0)
                    xd = (dataFrameCond['fk'] == fks[fki]) & (dataFrameCond['dk'] == dks[dki])

                    indices = dataFrameCond[xd].index
                    indices_norm = dataFrameCond[xnorm].index

                    data_cmass_slice = data_cmass[indices]
                    data_cmass_slice_norm = data_cmass[indices_norm]

                    if (fks[fki] == 0.8) & (dks[dki] == 8.0):
                        mean = np.mean(data_cmass_slice[biomass_indexes[:, gc], gc])
                        ax.text(dki / (numrows) + 0.5 / numrows,
                                fki / (numcols) + 0.5 / numcols,
                                np.round(mean, 2), horizontalalignment='center',
                                verticalalignment='center',
                                transform=ax.transAxes, size=10)

                    else:
                        mean = np.mean(data_cmass_slice[biomass_indexes[:, gc], gc])
                        mean_diff = np.mean(data_cmass_slice[biomass_indexes[:, gc], gc] -
                                            data_cmass_slice_norm[biomass_indexes[:, gc], gc] )
                        img[fki, dki] = mean_diff
                        ax.text(dki / (numrows) + 0.5 / numrows,
                                1 - (fki / (numcols) + 0.5 / numcols),
                                np.round(mean, 1), horizontalalignment='center',
                                verticalalignment='center',
                                transform=ax.transAxes, size=10)

            imsh = ax.imshow(img, cmap=cmap, norm= norm,  interpolation='nearest')
            ax.text(0.03, 0.94, lowerletters[index-1] + ")", horizontalalignment='left',
                       verticalalignment='center',
                       transform=ax.transAxes, size=12,
                       bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.set_xticks(np.arange(3), dks)
            ax.set_yticks(np.arange(3), fks)
            ax.set_ylabel('dk')
            ax.set_xlabel('fk')

            index += 1
        plt.subplots_adjust(top=0.8, left=0.15, wspace=0.25, hspace=0.2)
        cax = plt.axes([0.92, 0.11, 0.02, 0.69])
        bar = plt.colorbar(imsh, cax=cax, orientation="vertical")
        bar.set_label(r'$\Delta$' + var + ' [' + file.stratifiedOutput.GetUnit('Pft-Out', var) + ']', fontsize=12)
        plt.savefig("Sens_all/2004/" + var + ".png", dpi=300, bbox_inches='tight',pad_inches=0)
        plt.close()









