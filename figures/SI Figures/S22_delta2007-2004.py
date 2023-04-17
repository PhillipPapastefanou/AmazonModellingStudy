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


filename = "/Users/pp/data/Simulations/A03_Hydraulics/sens_1000/Review_Sensitivity/AnnuallyOutSens.nc"
config ="/Users/pp/data/Simulations/A03_Hydraulics/sens_1000/Review_Sensitivity/Values.tsv"

file = AmazonBasinStratified(filePath=filename, configPath=config, configheader=0)

vars = file.stratifiedOutput.PftNames()


#vars = ['aetplant']


b_s = [0.3, 0.5, 0.7]
numcols = len(b_s)
m_leaf_s = [1, 1.25, 1.5]
numrows = len(m_leaf_s)

vars.remove('cmass_litter')
vars.remove('lambda')
vars.remove('omega')

for var in vars:
    try:
        print(var)
        file = AmazonBasinStratified(filePath=filename, configPath=config, configheader=0)
        data2004 = file.GetDataSlicesRange(var, 2004, 2004, 0)
        data2007 = file.GetDataSlicesRange(var, 2008, 2008, 0)

        data_cmass = data2007 - data2004

        dataFrameCond = file.inputDataFrame

        img = np.zeros((6, 3, 3))


        mins = 10000000;
        maxs = -10000000;
        for gc in range(0, 6):
            for bi in range(0, len(b_s)):
                for m_leafi in range(0, len(m_leaf_s)):
                    xnorm = (dataFrameCond['b'] == 0.5) & (dataFrameCond['m_leaf'] == m_leaf_s[m_leafi]) & (dataFrameCond['m_root'] == 0.75)
                    xd = (dataFrameCond['b'] == b_s[bi]) & (dataFrameCond['m_leaf'] == m_leaf_s[m_leafi]) & (
                                dataFrameCond['m_root'] == 0.75)

                    indices = dataFrameCond[xd].index
                    indices_norm = dataFrameCond[xnorm].index

                    data_cmass_slice = data_cmass[indices]
                    data_cmass_slice_norm = data_cmass[indices_norm]

                    mean_diff = np.mean(data_cmass_slice[:, gc] - data_cmass_slice_norm[:,gc])

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



            for bi in range(0, len(b_s)):
                for m_leafi in range(0, len(m_leaf_s)):
                    xnorm = (dataFrameCond['b'] == 0.5) & (dataFrameCond['m_leaf'] == 1.5) & (dataFrameCond['m_root'] == 0.75)
                    xd = (dataFrameCond['b'] == b_s[bi]) & (dataFrameCond['m_leaf'] == m_leaf_s[m_leafi]) & (
                                dataFrameCond['m_root'] == 0.75)

                    indices = dataFrameCond[xd].index
                    indices_norm = dataFrameCond[xnorm].index

                    data_cmass_slice = data_cmass[indices]
                    data_cmass_slice_norm = data_cmass[indices_norm]

                    if (b_s[bi] == 0.5) & (m_leaf_s[m_leafi]== 1.5):
                        mean = np.mean(data_cmass_slice_norm[:, gc])
                        ax.text(m_leafi / (numrows) + 0.5 / numrows,
                                bi / (numcols) + 0.5 / numcols,
                                np.round(mean, 2), horizontalalignment='center',
                                verticalalignment='center',
                                transform=ax.transAxes, size=10)

                    else:
                        mean = np.mean(data_cmass_slice[:, gc])
                        mean_diff = np.mean(data_cmass_slice[:, gc] - data_cmass_slice_norm[:,gc])
                        img[bi, m_leafi] = mean_diff
                        ax.text(m_leafi / (numrows) + 0.5 / numrows,
                                1-(bi / (numcols) + 0.5 / numcols),
                                np.round(mean, 1) , horizontalalignment='center',
                                verticalalignment='center',
                                transform=ax.transAxes, size=10)

            imsh = ax.imshow(img, cmap=cmap, norm= norm,  interpolation='nearest')
            ax.text(0.03, 0.94, lowerletters[index-1] + ")", horizontalalignment='left',
                       verticalalignment='center',
                       transform=ax.transAxes, size=12,
                       bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.set_xticks(np.arange(3), m_leaf_s)
            ax.set_yticks(np.arange(3), b_s)
            ax.set_ylabel('b')
            ax.set_xlabel('m_leaf')

            index += 1
        plt.subplots_adjust(top=0.8, left=0.15, wspace=0.25, hspace=0.2)
        cax = plt.axes([0.92, 0.11, 0.02, 0.69])
        bar = plt.colorbar(imsh, cax=cax, orientation="vertical")
        bar.set_label(r'$\Delta$' + var + ' [' + file.stratifiedOutput.GetUnit('Pft-Out', var) + ']', fontsize=12)
        plt.savefig("Sens_all_leaf/diff2008-2004/diff_" + var + ".png", dpi=300, bbox_inches='tight',pad_inches=0)
        plt.close()


    except(RuntimeError):
        print(var + " failed!")
    except:
        print(var + " failed!")









