from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from Libs.Standard import StandardLPJG
from A02_MCWD_Dataset_Analysis.Pylibs.MCWDevaluation import MCWDFile as MF
from A02_MCWD_Dataset_Analysis.Pylibs.MCWD_Analysis21 import MCWDFile as MF_rel
from A02_MCWD_Dataset_Analysis.Pylibs.PrecAnomalyDry21 import PrecAnomaly
from A02_MCWD_Dataset_Analysis.Pylibs.scPDSI2021 import scPDSI

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



class Figure2_SingleAnomalies_Multi(BaseAnalysis):

    def p2009(self, x):
        return 0.1 * (0.3778 - 0.052 * x)

    def p2009_const(self, x):
        return 0.1 * (0.3778 - 0.052 * (-120.0))

    def makeColours(self, vals):
        colours = np.zeros((len(vals), 3))
        norm = Normalize(vmin=vals.min(), vmax=vals.max())

        # Can put any colormap you like here.
        colours = [cm.ScalarMappable(norm=norm, cmap='jet').to_rgba(val) for val in vals]

        return colours

    def f(self, x, a, b, c):
        return -np.exp(a * x + b) + c

    def sumOfSquaredError(self, parameterTuple, px, py):
        warnings.filterwarnings("ignore")  # do not print warnings by genetic algorithm
        val = self.f(px, *parameterTuple)
        return np.sum((py - val) ** 2.0)

    def predband(self, x, xd, yd, p, func, conf=0.95):
        # x = requested points
        # xd = x data
        # yd = y data
        # p = parameters
        # func = function name
        alpha = 1.0 - conf  # significance
        N = xd.size  # data sample size
        var_n = len(p)  # number of parameters
        # Quantile of Student's t distribution for p=(1-alpha/2)
        q = scipy.stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
        # Stdev of an individual measurement
        se = np.sqrt(1. / (N - var_n) * \
                     np.sum((yd - func(xd, *p)) ** 2))
        # Auxiliary definitions
        sx = (x - xd.mean()) ** 2
        sxd = np.sum((xd - xd.mean()) ** 2)
        # Predicted values (best-fit model)
        yp = func(x, *p)
        # Prediction band
        dy = q * se * np.sqrt(1.0 + (1.0 / N) + (sx / sxd))
        # Upper & lower prediction bands.
        lpb, upb = yp - dy, yp + dy
        return lpb, upb

    def latex_float(self, x, prec):
        exp = int(np.floor(np.log10(abs(x))))
        trail = np.round(x / 10 ** exp, prec)
        return r"{0} \times 10^{{{1}}}".format(trail, int(exp))

    def __init__(self, setup):

        from matplotlib import cm
        self.StartSetup()
        self.GetAGBFractionCMassTotal = np.vectorize(self.getAGBFractionCMassTotal)
        self.p2009v = np.vectorize(self.p2009)
        self.p2009_constv = np.vectorize(self.p2009_const)

        biomassConverter = BiomassAmazonAnomalyConverter()

        if platform == "darwin":
            raisg_mask = "/Users/pp/Dropbox/ClimateData/AmazonBasin/AB-SHAPE/Amazon_shape.shp"
            inputData = pd.read_csv('/Users/pp/Dropbox/UNI/Projekte/A03_Hydraulics_Implementation/Analysis/AllCavCurves05-2020.tsv',
                sep='\t',
                header=None).values
            phillips_points = pd.read_csv(
                "/Users/pp/Dropbox/UNI/Projekte/A03_Hydraulics_Implementation/AutoAnalysis/Source/Phillips2009-AGBchange.csv",
                sep=',',
                header=None).values
            posIndexes = np.squeeze(
                pd.read_csv(r'/Users/pp/Dropbox/UNI/Projekte/A03_Hydraulics_Implementation/EllipsisPositionInAB.tsv',
                            sep='\t',
                            header=None).values)

        else:
            raisg_mask = r"F:\Dropbox\ClimateData\AmazonBasin\AB-SHAPE\Amazon_shape.shp"
            inputData = pd.read_csv(
                r'F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\Analysis\AllCavCurves05-2020.tsv',
                sep='\t',
                header=None).values
            phillips_points = pd.read_csv(
                r"F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\AutoAnalysis\Source\Phillips2009-AGBchange.csv",
                sep=',',
                header=None).values
            posIndexes = np.squeeze(
                pd.read_csv(r'F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\EllipsisPositionInAB.tsv',
                            sep='\t',
                            header=None).values)


        # Convert from Mathematica indexing to python indexing
        posIndexes -= 1



        standardLPJs = []
        for file in setup.stdFilePaths:
            standardLPJs.append(StandardLPJG(file, True))

        mcwdFiles = []
        for file in setup.mcwdPaths:
            mcwdFiles.append(StandardLPJG(file, True))


        fileSCPDSI ='/Users/pp/Dropbox/ClimateData/AmazonBasin/SCPDSI/2005-2010/GLDAS20_AB_Monthly_05-scPDSI.txt'
        fileMCWD ='/Users/pp/Dropbox/ClimateData/AmazonBasin/MCWD/GLDAS20_AB_Yearly_05_MCWD.nc'
        filePREC ='/Users/pp/Dropbox/ClimateData/AmazonBasin/MonthlyPrec/GLDAS20_AB_Monthly_05.nc'

        scPDSIfile = scPDSI(fileSCPDSI)
        absSCPSIData = scPDSIfile.ParseRelativeDeviation(1950, 2010)

        fileMcwd2 = MF_rel(fileMCWD)
        absMcwdData = fileMcwd2.ParseRelativeDeviation(2001, 2010)

        precFile = PrecAnomaly(filePREC)
        absPrecData = precFile.ParseRelativeDeviation(2001, 2010)


        yearbeginMCWD = 2003
        yearEndMCWD = 2010

        yearbeginAGB = 1985
        yearEndAGB = 2010

        yois = [2005, 2007, 2009, 2010]

        dataMCWD = np.zeros((len(mcwdFiles), 1946, yearEndMCWD - yearbeginMCWD + 1))

        MCWDFiles = []
        for file in setup.mcwdPaths:
            MCWDFiles.append(MF(file))

        i = 0
        for MCWDFile in MCWDFiles:
            slice =  MCWDFile.ParseRange(yearbeginMCWD, yearEndMCWD, 2000, 2010)
            dataMCWD[i] = slice
            i += 1



        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))


        datasetNames = setup.hydFileNames
        paramData = []


        fig = plt.figure(figsize=(18, 9))
        index = 1



        #impacts = np.zeros((4, 500))

        rsquares = []

        for metric_index in range(0, 4):
            for d in range(0, 1):

                for yoi in yois:

                    datasetImpacts = []
                    rsquares_dat = []

                    biomass = listAnalysis[d].GetDataSlicesCondContinous("cmasstotal", yearbeginAGB, yearEndAGB, 1.0)
                    # standardLPJs[d].GetBiomassBeforeAfter(yearbeginAGB, yearEndAGB, True)
                    biomass_std = standardLPJs[d].GetDataSlicesCondContinous("cmasstotal", yearbeginAGB, yearEndAGB,
                                                                             1.0)

                    diff_anomaly = biomassConverter.ConvertToAbsoluteAnomaly_CVS(biomass, yoi, yearbeginAGB, yearEndAGB, yearbeginAGB)
                    diff_anomaly_std = biomassConverter.ConvertToAbsoluteAnomaly_CVS(biomass_std, yoi, yearbeginAGB, yearEndAGB,
                                                                yearbeginAGB)

                    # mcwdSlice = listAnalysis[d].GetMCWD(yearEnd)
                    mcwdSlice = dataMCWD[d, :, yoi - yearbeginMCWD]
                    rightText = ""

                    # diff[(dB < 2.0) | (diff > 0.5)] = np.nan
                    # diffStd[(dBStd < 2.0) | (diffStd > 0.5) ] = np.nan

                    ax = fig.add_subplot(len(yois), 4, index)

                    ptsLables = []
                    txtLables = []

                    sliceSpec = np.nanmean(diff_anomaly, axis=0)
                    sliceSpecStd = diff_anomaly_std

                    # posIndexes = np.arange(0,1946)

                    if metric_index == 0 :
                        mcwdSliced = mcwdSlice
                        p0 = [-0.01, -0.03, -2]
                        col = 'tab:blue'
                        xtext = 'absolute MCWD anomaly\n[mm $year^{-1}$]'

                    elif metric_index == 1:
                        mcwdSliced = absMcwdData[yoi - 2001]
                        p0 = [-2, -0.03, 0]
                        col = 'tab:red'
                        xtext = 'relative MCWD anomaly'

                    elif metric_index == 2:
                        mcwdSliced = absPrecData[yoi - 2001]
                        col = 'tab:green'
                        p0 = [-2, -0.03, 0]
                        xtext = 'rel dry season prec. anomaly\n[mm $year^{-1}$]'

                    elif metric_index == 3:
                        mcwdSliced = absSCPSIData[yoi - 2000]
                        col = 'black'
                        p0 = [-1, -0.03, 0]
                        xtext = 'rel scPDSI anomaly\n[$year^{-1}$]'

                    posIndexes = (biomass_std[:, 0] > 2.0) & (sliceSpec > -1000) & (np.isfinite(mcwdSliced))

                    sliceSpecd = sliceSpec[posIndexes] * 10.0
                    sliceSpecdStd = sliceSpecStd[posIndexes] * 10.0
                    mcwdSliced = mcwdSliced[posIndexes]

                    # sliceSpecd = sliceSpec*  10.0
                    # sliceSpecdStd = sliceSpecStd * 10.0
                    # mcwdSliced = mcwdSlice

                    dataXM = mcwdSliced
                    dataYM = sliceSpecd
                    dataYMSTD = sliceSpecdStd

                    ax.set_xlabel(xtext ,fontsize = 8)

                    if (index % 4 == 1):
                        ax.set_ylabel(r'AGB anomaly [MgC $ha^{-1}$  $year^{-1}$]', fontsize = 6)
                    # ax.text(0.97, 0.06, '${R}^2=$' + str(np.round(r_squared, 2)), horizontalalignment='right', transform=ax.transAxes,
                    #        size=12,
                    #        bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                    # ax.text(0.97, 0.17, r"$\mathrm{RMSE}_\mathrm{Hyd} = " + str(np.round(hyd_RMSE, 1)) + "$", horizontalalignment='right', verticalalignment='center',
                    #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'), c = "tab:blue")
                    # ax.text(0.97, 0.07, r"$\mathrm{RMSE}_\mathrm{Std} = " + str(np.round(std_RMSE, 1)) + "$", horizontalalignment='right', verticalalignment='center',
                    #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'), c = "tab:purple", alpha = 0.8)
                    ax.text(0.03, 0.94, lowerletters[index - 1] + ") " + str(yoi), horizontalalignment='left',
                            verticalalignment='center',
                            transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                    ax.text(0.15, 0.95, rightText, horizontalalignment='left', verticalalignment='center',
                            transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                    if metric_index == 0:
                        ax.set_xlim([-250, 50])
                        xStd = np.linspace(-120, 40, 10)
                        xExt = np.linspace(-240, -120, 10)

                    elif metric_index == 1:
                        xStd = np.linspace(np.min(mcwdSliced), np.max(mcwdSliced), 10)
                        xExt = np.linspace(np.min(mcwdSliced), np.max(mcwdSliced), 10)

                    ax.set_ylim([-20, 5])



                    pxxmin = np.min(mcwdSliced)
                    pxx = np.linspace(pxxmin, np.max(mcwdSliced), 100)

                    popt, pcov = curve_fit(self.f, dataXM, dataYMSTD, p0=np.asarray([-0, -0.05, -0]), maxfev=100000)
                    a = popt[0]
                    b = popt[1]
                    c = popt[2]

                    pyy_std = self.f(pxx, *popt)
                    std_RMSE = mean_squared_error(self.f(-phillips_points[:, 0], a, b, c), phillips_points[:, 1] * 0.5,
                                                  squared=False)

                    popt, pcov = curve_fit(self.f, dataXM, dataYM, p0=np.asarray(p0), maxfev=1000000)
                    a = popt[0]
                    b = popt[1]
                    c = popt[2]

                    hyd_RMSE = mean_squared_error(self.f(dataXM, *popt), dataYM,
                                                  squared=False)

                    # lpb, upb = self.predband(pxx,dataXM, dataYM, *popt, 0.95)
                    # lpb, upb = self.predband(pxx, dataXM, dataYM, popt, self.f, conf=0.95)

                    pyy = self.f(pxx, *popt)

                    ptsLable = ax.scatter(dataXM, dataYM, color=col, alpha=0.7, s=0.5)
                    ptsLable2 = ax.scatter(dataXM, dataYMSTD, color='tab:purple', alpha=0.3, s=0.5)
                    ptsLables.append(ptsLable)
                    ptsLables.append(ptsLable2)

                    ax.plot(pxx, pyy, linewidth=1.5, alpha=1, color=col)
                    # plt.plot(pxx, lpb, 'k--', linewidth=0.75, alpha=1, color = "tab:blue")
                    # plt.plot(pxx, upb, 'k--', linewidth=0.75, alpha=1, color = "tab:blue")
                    ax.plot(pxx, pyy_std, linewidth=1.0, alpha=0.6, color="tab:purple")

                    ax.text(0.97, 0.17, r"$\mathrm{RMSE}_\mathrm{Hyd} = " + str(np.round(hyd_RMSE, 1)) + "$",
                            horizontalalignment='right', verticalalignment='center',
                            transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'),
                            c=col)

                    if i == 37:
                        lineStd, = ax.plot(pxx, pyy, linewidth=1.5, color='tab:purple', alpha=1.0)
                    rsquares.append(rsquares_dat)

                    if (yoi == 2005) & (metric_index == 0):
                        pobns = ax.scatter(-phillips_points[:, 0], phillips_points[:, 1] * 0.5, color='black', s=10)

                    index += 1






        txtLables.append("Obs $(n=22)$")





        #for i in range(0, len(legend.legendHandles)-1):
       #     legend.legendHandles[i]
       # legend.legendHandles[len(legend.legendHandles)-1].set_sizes([15.0])


        #for handle in legend2.legendHandles:
            #handle.set_sizes([15.0])


        plt.subplots_adjust(right=0.45, wspace=0.3, hspace=0.5)

       #bar.set_label('$\Delta$AGB [MgC/ha]', fontsize=12)



        df = pd.DataFrame(paramData)

        import os
        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)
        plt.savefig(setup.folderPath + "/" + "Figure2_drought_relationships.png", dpi=600, bbox_inches='tight',
                    pad_inches=0)

        plt.close()

        self.Succesfull()


