from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from Libs.Standard import StandardLPJG
from A02_MCWD_Dataset_Analysis.Pylibs.MCWDevaluation import MCWDFile as MF
from A03_Hydraulics_Implementation.AutoAnalysis.Source.SetupFiles import SetupFile
from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis
from sys import platform
from A02_MCWD_Dataset_Analysis.Setup2010 import Setup2010
from A02_MCWD_Dataset_Analysis.Pylibs.MCWD_Analysis21 import MCWDFile as MF2
from matplotlib.legend_handler import HandlerBase

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
from sklearn.metrics import mean_squared_error
import scipy.stats
import warnings
import  pandas as pd

from scipy.stats import gaussian_kde as kde
from matplotlib.colors import Normalize
from matplotlib import cm

from sklearn.metrics import mean_squared_error


class HandlerBoxPlot(HandlerBase):
    def create_artists(self, legend, orig_handle,
                   xdescent, ydescent, width, height, fontsize,
                   trans):

        width = 20
        height = 4
        xdescent = 0
        lw = 0.75
        a_list = []

        a_list.append(matplotlib.patches.Rectangle((0.25 *width-xdescent,0.25 *height-xdescent),
                                              0.5*width-ydescent, 0.5*height-ydescent, lw = lw,
                                                   facecolor='none')) #

        a_list.append(matplotlib.lines.Line2D(np.array([0.25, 0.25, 0.75, 0.75, 0.25])*width-xdescent,
                                              np.array([0.25, 0.75, 0.75, 0.25, 0.25])*height-ydescent, lw = lw)) # box

        a_list.append(matplotlib.lines.Line2D(np.array([0.75,1.0])*width-xdescent,
                                              np.array([0.5,0.5])*height-ydescent, lw = lw)) # top vert line

        a_list.append(matplotlib.lines.Line2D(np.array([0,0.25])*width-xdescent,
                                              np.array([0.5,0.5])*height-ydescent, lw = lw)) # bottom vert line

        a_list.append(matplotlib.lines.Line2D(np.array([1,1])*width-xdescent,
                                               np.array([0.25,0.75])*height-ydescent, lw = lw)) # top vert line

        a_list.append(matplotlib.lines.Line2D(np.array([0,0])*width-xdescent,
                                               np.array([0.25,0.75])*height-ydescent, lw = lw)) # top vert line

        for a in a_list:
            a.set_color(orig_handle.get_color())
        return a_list

class Figure2_SingleAnomalies_Yang_CAX_smth(BaseAnalysis):

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

    # Exponential model to be fitted for each PWS
    def f(self, x, a, b, c):
        return -np.exp(a * x + b) + c

    def sumOfSquaredError(self, parameterTuple, px, py):
        warnings.filterwarnings("ignore")  # do not print warnings by genetic algorithm
        val = self.f(px, *parameterTuple)
        return np.sum((py - val) ** 2.0)

    # Prediction band algorithm
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

    # Forumula to calculate the MCWD
    # This routine is needed to esimate the MCWD of the Caxiuana experiment
    def MCWD(self,monthly_prec):
        cwd = 0
        mcwd = 0
        for i in range(monthly_prec.shape[0]):
            wd = monthly_prec[i] - 100
            if wd < 0:
                cwd += wd
            if mcwd > cwd:
                mcwd = cwd
        return mcwd

    def __init__(self, setup):
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        # Get the data from the experiments
        df_normal_prec = pd.read_csv('../MCWD_CAX/CAX_Monthly_PRec.csv', header=None)
        df_abg_rowload = pd.read_csv('../MCWD_CAX/RowlandAGB_TFE.csv', header=None)

        mcwds_normal = np.zeros(7)
        for k in range(0, 7):
            prec = df_normal_prec.values[9 + 12 * k:21 + 12 * k, 1]
            mcwds_normal[k] = self.MCWD(prec)

        # Reduce precipitation according to the application of the throughfall exclusion pads
        mean_MCWD_normal = np.mean(mcwds_normal)
        reduced_prec = df_normal_prec.values

        indexes = np.floor(reduced_prec[:, 0]) == 2004
        reduced_prec[indexes, 1] *= 0.7

        indexes = (np.floor(reduced_prec[:, 0]) != 2004) & (reduced_prec[:, 0] > 2002)
        reduced_prec[indexes, 1] *= 0.5

        mcwds_reduced = np.zeros(7)
        for k in range(0, 7):
            prec = reduced_prec[9 + 12 * k:21 + 12 * k, 1]
            mcwds_reduced[k] = self.MCWD(prec)

        mcwd_reduced_anomalies = mcwds_reduced - mean_MCWD_normal

        diffs = np.diff(df_abg_rowload.values[0:7, 1])
        diffs_upper = np.diff(df_abg_rowload.values[7:14, 1])
        diffs_lower = np.diff(df_abg_rowload.values[14:21, 1])
        df_abg_rowload.values[7, 1] / 2.0 - df_abg_rowload.values[0, 1] / 2.0


        from matplotlib import cm
        self.StartSetup()
        self.GetAGBFractionCMassTotal = np.vectorize(self.getAGBFractionCMassTotal)
        self.p2009v = np.vectorize(self.p2009)
        self.p2009_constv = np.vectorize(self.p2009_const)

        if platform == "darwin":
            raisg_mask = "../AB-SHAPE/Amazon_shape.shp"
            inputData = pd.read_csv('../AllCavCurves05-2020.tsv',
                sep='\t',
                header=None).values
            phillips_points = pd.read_csv(
                "../Phillips2009-AGBchange.csv",
                sep=',',
                header=None).values
            posIndexes = np.squeeze(
                pd.read_csv(r'/Users/pp/Dropbox/UNI/Projekte/A03_Hydraulics_Implementation/EllipsisPositionInAB.tsv',
                            sep='\t',
                            header=None).values)

        impact2005_yang = pd.read_csv(
            '../Yang2018_ClassesAGB_change.csv',
            sep=',', header=None)

        # Impact of drought stress according to the Yang et al. 2018 study
        vals = impact2005_yang.values
        vals_SD = np.diff(vals[0:6,1])
        vals_MD = np.diff(vals[6:12,1])
        vals_ED = np.diff(vals[12:18,1])
        impact_SD = vals_SD[2] - np.mean(np.delete(vals_SD,2))
        impact_MD = vals_MD[2] - np.mean(np.delete(vals_MD,2))
        impact_ED  = vals_ED[2] - np.mean(np.delete(vals_ED,2))
        impact_SD /= 2.0
        impact_MD /= 2.0
        impact_ED /= 2.0

        years = [2005]
        def conditions(arr, i):
            if i == 1:
                return (arr < -0.5) & (arr > -2)
            if i == 2:
                return (arr <= -2) & (arr > -2.5)
            if i == 3:
                return (arr <= -2.5)
            return -1

        fileNames = ["DOLCE", "GLEAM"]
        data2005 = np.zeros((2, 1946))

        e_index = 0

        MCWDFilesConst = []
        scPDSIFiles = []

        setup2 = Setup2010()

        MCWDFilesConst.append(MF2(setup2.MCWDrootPAth + "/" + setup2.files[6][1]))

        i = 0
        for x in MCWDFilesConst:
            data = x.ParseRelativeDeviation(2000, 2009)
            data_abs = x.ParseAbsoluteDeviation(2000, 2009)
            for y in range(0, len(years)):
                data2005[0] = data[years[y] - 2000]
                data2005[1] = data_abs[years[y] - 2000]

        e_index += 1

        df = pd.DataFrame(data2005.T, columns=['rel', 'abs'])
        df['class'] = np.nan
        indexes = df['rel'] < -2.0
        df['class'].iloc[indexes] = 'ED'
        indexes = (df['rel'] > -2.0) & (df['rel'] < - 1.0)
        df['class'].iloc[indexes] = 'SD'
        indexes = (df['rel'] > -1.0) & (df['rel'] < - 0.5)
        df['class'].iloc[indexes] = 'MD'
        c = 'tab:red'
        #ax.boxplot(df[df['class'] == 'ED']['abs']
        #            , positions=[impact_ED], vert=False, patch_artist=True,
        #            boxprops=dict(facecolor=c, color=c)
        #            )

        # Convert from Mathematica indexing to python indexing
        posIndexes -= 1

        standardLPJs = []
        for file in setup.stdFilePaths:
            standardLPJs.append(StandardLPJG(file, True))

        mcwdFiles = []
        for file in setup.mcwdPaths:
            mcwdFiles.append(StandardLPJG(file, True))

        data2005 = np.zeros((len(mcwdFiles), 1946))
        data2010 = np.zeros((len(mcwdFiles), 1946))

        MCWDFiles = []
        for file in setup.mcwdPaths:
            MCWDFiles.append(MF(file))

        i = 0
        for MCWDFile in MCWDFiles:
            MCWDFile.Parse(2000, 2010, [2005, 2010])
            data2005[i] = MCWDFile.DataSlices[0]
            data2010[i] = MCWDFile.DataSlices[1]
            i += 1



        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))


        normPsi = mpl.colors.Normalize(-2.0, 2.0)
        colorsPsi = [[normPsi(-2.0), "tab:red"],
                     [normPsi(-1.0), "tab:orange"],
                     [normPsi(-0.1), "white"],
                     [normPsi(0.0), "white"],
                     [normPsi(0.1), "white"],
                     [normPsi(2.0), "tab:green"]]
        cmapPsi = mpl.colors.LinearSegmentedColormap.from_list("", colorsPsi)

        norm = mpl.colors.Normalize(-180, 0.0)
        colorsMCWD = [[norm(-180.0), (0.65, 0.16, 0.16)],
                      [norm(-100.0), (0.80, 0.58, 0.047)],
                      [norm(-25.0), "white"],
                      [norm(0.0), "white"], ]
        cmapMCWD = mpl.colors.LinearSegmentedColormap.from_list("", colorsMCWD)


        paramData = []
        yearbegin = 1985
        yearEnd = 2005

        fig = plt.figure(figsize=(7, 6))
        index = 0

        import string
        lowerletters = string.ascii_lowercase[:26]

        pts = 1830
        impacts = np.zeros((4, pts))
        dataIndividual = np.zeros(((37+1)* 2, pts))
        rsquares = []
        indexes = [1, 2]

        for d in range(0, 1):

            datasetImpacts = []
            rsquares_dat =[]

            biomass = listAnalysis[d].GetDataSlicesCondContinous("cmasstotal", yearbegin, yearEnd, 0)
            standardLPJs[d].GetBiomassBeforeAfter(yearbegin, yearEnd, True)

            dBStd = standardLPJs[d].BiomassBefore
            dAStd = standardLPJs[d].BiomassAfter

            fracAGB = self.GetAGBFractionCMassTotal(biomass)
            fracAGBstd = self.GetAGBFractionCMassTotal(dBStd)

            biomassAGB = biomass * fracAGB
            biomassAGBDiff = np.dff(biomassAGB, axis = 2)

            diff = biomassAGBDiff[:,:,2005-yearbegin-1] - np.mean(biomassAGBDiff[:,:,0:2004-yearbegin-2], axis= 2)
            diffStd = dAStd - dBStd
            diffStd *= fracAGBstd

            mcwdSlice = data2005[d]
            rightText = ""

            #
            # TOP PLOT
            #

            ax = fig.add_subplot(2, 1, indexes[index])
            cool = 'tab:red'
            boxplt  = ax.boxplot(df[df['class'] == 'ED']['abs']
                        , positions=[impact_ED], vert=False, patch_artist=True,
                        boxprops=dict(facecolor=cool, color='black'), widths=0.6,labels=['']
                        , zorder =2, manage_ticks=False,)
            c = 'tab:red'
            ax.boxplot(df[df['class'] == 'SD']['abs'], positions=[impact_SD], vert=False, patch_artist=True,
                        boxprops=dict(facecolor=c, color='black'), widths=0.6,labels=[''],
                        flierprops = dict(marker='o',markeredgecolor= 'gray', markersize=2, linestyle='none'), zorder = 2, manage_ticks=False)
            c = 'tab:red'
            ax.boxplot(df[df['class'] == 'MD']['abs'], positions=[impact_MD], vert=False, patch_artist=True,
                        boxprops=dict(facecolor=c, color='black'), widths=0.6,labels=[''], zorder = 2, manage_ticks=False)


            ptsLables = []
            txtLables = []

            w = biomassAGB[:,:,0]

            diff = diff[:,posIndexes]
            w = w[:,posIndexes]

            sliceSpec = np.average(diff, axis = 0, weights=w) *10.0
            sliceSpecStd = diffStd
            sliceSpecdStd = sliceSpecStd[posIndexes] * 10.0


            mcwdSliced = mcwdSlice[posIndexes]
            nonNanIndexes = (np.logical_not(np.isnan(mcwdSliced))) \
                            & (np.logical_not(np.isnan(sliceSpec)))

            dataXM = mcwdSliced[nonNanIndexes]
            dataYM = sliceSpec[nonNanIndexes]

            nonNanIndexes_std = (np.logical_not(np.isnan(mcwdSliced))) \
                            & (np.logical_not(np.isnan(sliceSpecdStd)))

            dataXM_std = mcwdSliced[nonNanIndexes_std]
            dataYM_std = sliceSpecdStd[nonNanIndexes_std]


            # MCWD anonamly linespace
            pxx = np.linspace(-230, 50, 100)

            # Fit the exponential function to the standard model (without hydraulics)
            popt, pcov = curve_fit(self.f, dataXM_std, dataYM_std, p0=np.asarray([-0, -0.05, -0]), maxfev=100000)
            a = popt[0]
            b = popt[1]
            c = popt[2]

            pyy_std = self.f(pxx, *popt)
            std_RMSE = mean_squared_error(self.f(-phillips_points[:, 0], a, b, c),  phillips_points[:, 1] * 0.5, squared=False)

            # Fit the exponential function to each of the PWS
            popt, pcov = curve_fit(self.f, dataXM, dataYM, p0=np.asarray([-0.01, -0.03, -2]), maxfev=100000)
            a = popt[0]
            b = popt[1]
            c = popt[2]
            residuals = dataYM - self.f(dataXM, a, b, c)
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((dataYM - np.mean(dataYM)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            dP = np.sqrt(np.diag(pcov))

            hyd_RMSE = mean_squared_error(self.f(-phillips_points[:, 0], a, b, c),  phillips_points[:, 1] * 0.5, squared=False)


            lpb, upb = self.predband(pxx, dataXM, dataYM, popt, self.f, conf=0.95)

            pyy = self.f(pxx, *popt)

            ptsLable =  ax.scatter(dataXM, dataYM, color='tab:blue', alpha = 0.7, s=0.1)
            ptsLable_std =  ax.scatter(dataXM_std, dataYM_std, color='tab:purple', alpha = 0.3, s=0.1)

            ptsLables.append(ptsLable)
            ptsLables.append(ptsLable_std)


            ax.plot(pxx, pyy, linewidth=1.5, alpha=1, color = "tab:blue", zorder = 2)
            dummy,  = ax.plot(pxx, pyy, linewidth=1.5, alpha=1, color = "black", zorder = 0)
            ax.plot(pxx, lpb, 'k--', linewidth=0.75, alpha=1, color = "tab:blue", zorder = 2)
            ax.plot(pxx, upb, 'k--', linewidth=0.75, alpha=1, color = "tab:blue", zorder = 2)
            ax.plot(pxx, pyy_std, linewidth=1.0, alpha=0.6, color = "tab:purple", zorder = 2)


            if i == 37:
                lineStd, = ax.plot(pxx, pyy, linewidth=1.5, color='tab:purple', alpha=1.0)
            rsquares.append(rsquares_dat)


            ax.set_xlabel(r'MCWD anomaly [mm $year^{-1}$]')
            ax.set_ylabel(r'AGB anomaly [MgC $ha^{-1}$  $year^{-1}$]')
            ax.text(0.97, 0.17, r"$\mathrm{RMSE}_\mathrm{Hyd} = " + str(np.round(hyd_RMSE, 1)) + "$", horizontalalignment='right', verticalalignment='center',
                    transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'), c = "tab:blue")
            ax.text(0.97, 0.07, r"$\mathrm{RMSE}_\mathrm{Std} = " + str(np.round(std_RMSE, 1)) + "$", horizontalalignment='right', verticalalignment='center',
                    transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'), c = "tab:purple", alpha = 0.8)
            ax.text(0.03, 0.94, lowerletters[indexes[index]-1] + ")  mean over all species", horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.text(0.15, 0.95, rightText, horizontalalignment='left', verticalalignment='center',
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.set_xlim([-250, 50])
            ax.set_ylim([-15, 5])
            ax.set_yticks([5, -0, -5, -10, -15])

            xStd = np.linspace(-120, 40, 10)
            xExt = np.linspace(-240, -120, 10)

            # 13/04/2021 fix from Mg/m² tö MgC/m²
            pobns = ax.scatter(-phillips_points[:, 0], phillips_points[:, 1]*0.5 , color='black', s=10, zorder = 3)
            errs= ax.errorbar(mcwd_reduced_anomalies[1:], diffs,
                         yerr=df_abg_rowload.values[7, 1] / 2.0 - df_abg_rowload.values[0, 1] / 2.0, c='tab:brown',
                         fmt='s', label = 'Drought experiment derived',
                         capsize=2, markersize = 2, lw =1 )

            index += 1


            #
            # BOTTOM PLOT
            #
            ax = fig.add_subplot(2, 1, indexes[index])

            #minpsi50 = np.min(inputData[:, 0])
            minpsi50 = -4.5
            maxpsi50 = -1.5

            from matplotlib import cm

            viridis = cm.get_cmap('viridis', 60)

            normPsi = mpl.colors.Normalize(minpsi50, maxpsi50)
            colorsPsi = [[normPsi(minpsi50), "#111D5E"],
                         [normPsi(-3.0), "#C70039"],
                         [normPsi(maxpsi50), "#C0E218"]]
            colorsPsi = [[normPsi(minpsi50), "tab:red"],
                         [normPsi(maxpsi50), "tab:gray"]]
            cmapPsi = mpl.colors.LinearSegmentedColormap.from_list("", colorsPsi)

            dataY = []
            indexesOI = np.arange(0, 37)

            ptsLables = []
            txtLables = []
            for i in indexesOI:
                if i != 37:
                    sliceSpec = diff[i, :]
                else:
                    sliceSpec = diffStd


                sliceSpecd = sliceSpec
                mcwdSliced = mcwdSlice

                mcwdSliced = mcwdSlice[posIndexes]

                nonNanIndexes = (np.logical_not(np.isnan(mcwdSliced))) \
                                & (np.logical_not(np.isnan(sliceSpecd)))


                py = sliceSpecd[nonNanIndexes] *10
                px = mcwdSliced[nonNanIndexes]



                # Sanity check to avoid crashing of the routine
                if py.shape[0] > 470:

                    popt, pcov = curve_fit(self.f, px, py, p0=np.asarray([-0.04, -0.23, +2]), maxfev=100000)
                    a = popt[0]
                    b = popt[1]
                    c = popt[2]
                    residuals = py - self.f(px, a, b, c)
                    ss_res = np.sum(residuals ** 2)
                    ss_tot = np.sum((py - np.mean(py)) ** 2)
                    r_squared = 1 - (ss_res / ss_tot)
                    dP = np.sqrt(np.diag(pcov))
                    paramData.append([datasetNames[d], rightText, a, b, c])

                    pxx = np.linspace(-250, 50, 100)
                    pyy = self.f(pxx, *popt)

                    rsquares_dat.append(r_squared)

                    if ((i != 37)):
                        ax.scatter(px, py, s=1, alpha=0.0, edgecolors=None)
                        ptsLable = ax.plot(pxx, pyy, linewidth=1, alpha=1, c = viridis(normPsi(inputData[i, 0])))
                        ptsLables.append(ptsLable)
                        dataY.append(pyy)

            rsquares.append(rsquares_dat)

            ax.set_xlabel(r'$\Delta$MCWD [mm $year^{-1}$]')
            ax.set_ylabel(r'$\Delta$AGB [MgC $ha^{-1}$  $year^{-1}$]')
            ax.text(0.03, 0.94, lowerletters[indexes[index]-1] + ")  Individual", horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.text(0.15, 0.95, rightText, horizontalalignment='left', verticalalignment='center',
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.set_xlim([-250, 50])
            ax.set_ylim([-15, 5])

            xStd = np.linspace(-120, 40, 10)
            xExt = np.linspace(-240, -120, 10)

            index += 1



        j = 0
        for i in indexesOI:
            lablesTxt = r"$\psi_{50} = " + str(np.round(inputData[i, 0], 1)) + "$ MPa\n" + \
                        r"$\psi_{88} = " + str(np.round(inputData[i, 1], 1)) + "$ MPa"#\
            txtLables.append(lablesTxt)
            j +=1


        ptsLables.append(pobns)
        txtLables.append("Obs $(n=22)$")

        from matplotlib.lines import Line2D

        le_pts = Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Inventories $(n=22)$',
                        markersize=6)

        le_hyd = [Line2D([0], [0],marker='o' ,color='w',markerfacecolor='tab:blue', label='Modelled', markersize=0.75),
                           Line2D([0], [0], color='tab:blue', label='Exp-Fit', markersize=15),
                           Line2D([0], [0], color='tab:blue', label='95%-CI-bands', ls='--',
                                  markersize=15)
                           ]

        le_std = [         Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:purple',label='Modelled', markersize=0.75),
                           Line2D([0], [0], color='tab:purple', label='Exp-Fit', markersize=15)
                           ]

        legend1 = fig.legend(handles= [le_pts, errs, dummy], labels= ['Inventories','Drought experiment derived',
                                                                      'LiDAR derived'],
                  loc='center left', handler_map={dummy:HandlerBoxPlot()},
                  markerscale= 1,
                  bbox_to_anchor=(0.455, 0.83),
                   fontsize=8,
                   ncol=1,
                   fancybox=False, shadow=False)

        legend1.legendHandles[2].set_color('tab:red')

        legend2 = fig.legend(handles=le_hyd,
                  loc='center left',
                  title="$\\bf{LPJ-GUESS-HYD}$",
                  markerscale=5,
                  bbox_to_anchor=(0.455, 0.71),
                   fontsize=8,
                   ncol=1,
                   fancybox=False, shadow=False)
        legend2.get_title().set_fontsize('8')  # legend 'Title' fontsize


        legend3 = fig.legend(handles=le_std,
                  loc='center left',
                  title="$\\bf{LPJ-GUESS}$",
                  markerscale=5,
                  bbox_to_anchor=(0.455, 0.59),
                   fontsize=8,
                   ncol=1,
                   fancybox=False, shadow=False)

        legend3.get_title().set_fontsize('8')  # legend 'Title' fontsize


        plt.subplots_adjust(right=0.45, wspace=0.3, hspace=0.25)
        c_map_ax = fig.add_axes([0.47, 0.111, 0.01, 0.34])
        c_map_ax.axes.get_xaxis().set_visible(True)
        c_map_ax.axes.get_yaxis().set_visible(True)

        bar = mpl.colorbar.ColorbarBase(c_map_ax,cmap= viridis, norm = normPsi, orientation="vertical", label = r'$\psi_{50}$ [MPa]')

        df = pd.DataFrame(paramData)

        import os
        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)
        plt.savefig(setup.folderPath + "/" + "Figure2_single-" + str(yearEnd) + "_CAX_YANG.png", dpi=600, bbox_inches='tight',
                    pad_inches=0)

        plt.close()



        df.to_csv(setup.folderPath + "/ParametersFit" + str(yearEnd) + ".tsv")
        self.Succesfull()


