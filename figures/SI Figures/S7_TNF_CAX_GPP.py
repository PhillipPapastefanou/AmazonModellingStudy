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


class FigureS_TNF_CAX_GPP_MONTHLY(BaseAnalysis):


    def __init__(self, setup):

        cols = ['tab:blue', 'tab:green']

        locations = ['TNF', 'CAX', 'K34']

        standardLPJs = []
        for file in setup.stdFilePaths:
            standardLPJs.append(StandardLPJG(file, True))

        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))


        obs_CAX =pd.read_csv('/Users/pp/Documents/GPP_CAX.csv', sep=',',
            header=None).values


        for var in ['gpp']:

        #for var in listAnalysis[0].stratifiedOutput.PftNames():
            fig = plt.figure(figsize=(13, 13))

            index = 1

            data_first = []
            data_range = []

            for d in range(0, 1):
                # To convert from kgC/month gC/day
                dataf = listAnalysis[d].GetDataSlicesYear(var, 1996) * 1000/30
                dataR = listAnalysis[d].GetDataSlicesRange(var, 1996, 2010) * 1000/30
                data_first.append(dataf)
                data_range.append(dataR)



            for gc in range(0, 2):

                ax = fig.add_subplot(3,2, 2 *index -1)
                x = np.arange(1996, 2011, 1/12)

                for d in range(0, 1):
                    data_first_gc = data_first[d][:, gc]
                    data_first_gc = np.reshape(data_first_gc, (2, -1))

                    data_range_gc = data_range[d][:, gc, :]
                    data_range_gc = np.reshape(data_range_gc, (2, 37, -1))

                    data_year_c = data_first_gc[0]

                    data_range_c = data_range_gc[0]
                    data_range_d = data_range_gc[1]



                    for i in range(0, data_range_c.shape[0]):
                        ax.plot(x, data_range_c[i], c='black', lw=1, alpha=0.1)
                    ax.plot(x, np.mean(data_range_c, axis=0), c=cols[d])
                    ax.plot(x, np.quantile(data_range_c, axis=0, q=0.25), c=cols[d], dashes=[2, 2])
                    ax.plot(x, np.quantile(data_range_c, axis=0, q=0.75), c=cols[d], dashes=[2, 2])

                # ax.set_ylim((5, 25))
                ax.set_xlabel("Time")
                if gc == 0:
                    ax.set_xlim((1998, 2005))
                else:
                    ax.set_xlim((2001, 2010))

                ax.set_ylabel("GPP [g C $\mathrm{d}^{-1}$]")
                ax.text(0.04, 0.96, lowerletters[2* index - 2 ] + ")", horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes, size=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                ax.text(0.96, 0.96, locations[gc], horizontalalignment='right',
                        verticalalignment='top',
                        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                ax.text(0.96, 0.04, 'control', horizontalalignment='right',
                        verticalalignment='bottom',
                        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                if gc == 1:
                    ax.plot(obs_CAX[:, 0], obs_CAX[:, 1], c='black')

                ax = fig.add_subplot(3,2, 2 *index )

                for d in range(0, 1):
                    data_first_gc = data_first[d][:, gc]
                    data_first_gc = np.reshape(data_first_gc, (2, -1))

                    data_range_gc = data_range[d][:, gc, :]
                    data_range_gc = np.reshape(data_range_gc, (2, 37, -1))

                    data_year_c = data_first_gc[0]

                    data_range_c = data_range_gc[0]
                    data_range_d = data_range_gc[1]

                    for i in range(0, data_range_c.shape[0]):
                        ax.plot(x, data_range_d[i] - data_range_c[i], c='black', lw=1, alpha=0.1)
                    ax.plot(x, np.mean(data_range_d - data_range_c, axis=0), c=cols[d])
                    ax.plot(x, np.quantile(data_range_d - data_range_c, axis=0, q=0.25), c=cols[d], dashes=[2, 2])
                    ax.plot(x, np.quantile(data_range_d - data_range_c, axis=0, q=0.75), c=cols[d], dashes=[2, 2])

                # ax.set_ylim((-12,  0.5))
                ax.set_xlabel("Time")
                ax.set_ylabel("GPP [gC $\mathrm{day}^{-1}$]")
                ax.text(0.04, 0.96, lowerletters[2 *index - 1] + ")", horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes, size=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                ax.text(0.96, 0.96, locations[gc], horizontalalignment='right',
                        verticalalignment='top',
                        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                ax.text(0.96, 0.04, 'drought - control', horizontalalignment='right',
                        verticalalignment='bottom',
                        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                if gc == 0:
                    ax.set_xlim((1998, 2005))
                else:
                    ax.set_xlim((2001, 2010))

                index += 1

            plt.subplots_adjust(top=0.8, left=0.15, wspace=0.25, hspace=0.2)

            import os
            if not os.path.exists(setup.folderPath + "/SelectedVars/"):
                os.mkdir(setup.folderPath + "/SelectedVars/")

            plt.savefig(setup.folderPath + "/SelectedVars/" + "FigureS_TNF_CAX_monthly_" + var + ".png", dpi=600, bbox_inches='tight',
                        pad_inches=0)

            plt.close()








