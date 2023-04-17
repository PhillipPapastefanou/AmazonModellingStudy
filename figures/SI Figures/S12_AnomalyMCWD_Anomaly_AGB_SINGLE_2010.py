from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from Libs.Standard import StandardLPJG
from A02_MCWD_Dataset_Analysis.Pylibs.MCWDevaluation import MCWDFile as MF
from A03_Hydraulics_Implementation.AutoAnalysis.Source.SetupFiles import SetupFile
from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis

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


class Figure2_SingleAnomalies_2010(BaseAnalysis):

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

        raisg_mask = r"F:\Dropbox\ClimateData\AmazonBasin\AB-SHAPE\Amazon_shape.shp"
        mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                              ccrs.PlateCarree())

        inputData = pd.read_csv(
            r'F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\Analysis\AllCavCurves05-2020.tsv',
            sep='\t',
            header=None).values


        phillips_points =  pd.read_csv(r"F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\AutoAnalysis\Source\Phillips2009-AGBchange.csv",
                                  sep=',',
                                  header=None).values


        #posIndexes =   np.squeeze(pd.read_csv(r'F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\EllipsisPositionInAB.tsv',
        #                             sep='\t',
         #                            header=None).values)
       ## # Convert from Mathematica indexing to python indexing
        #posIndexes -= 1

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


        datasetNames = setup.hydFileNames

        paramData = []
        yearbegin = 1985
        yearEnd = 2010

        fig = plt.figure(figsize=(7, 6))
        index = 0

        import string
        lowerletters = string.ascii_lowercase[:26]

        impacts = np.zeros((4, 1946))

        dataIndividual = np.zeros(((37+1)* 2, 1946))
        rsquares = []

        indexes = [1,2]

        for d in range(0, 1):

            datasetImpacts = []

            rsquares_dat =[]

            biomass = listAnalysis[d].GetDataSlicesCondContinous("cmasstotal", yearbegin, yearEnd, 1.0, 0)
            standardLPJs[d].GetBiomassBeforeAfter(yearbegin, yearEnd, True)
            dBStd = standardLPJs[d].BiomassBefore
            dAStd = standardLPJs[d].BiomassAfter

            fracAGB = self.GetAGBFractionCMassTotal(biomass)
            fracAGBstd = self.GetAGBFractionCMassTotal(dBStd)

            biomassAGB = biomass * fracAGB

            biomassAGBDiff = np.diff(biomassAGB, axis = 2)


            diff = biomassAGBDiff[:,:,2010-yearbegin-1] - np.mean(biomassAGBDiff[:,:,0:2009-yearbegin-2], axis= 2)



            diffStd = dAStd - dBStd
            diffStd *= fracAGBstd

            # mcwdSlice = listAnalysis[d].GetMCWD(yearEnd)
            mcwdSlice = data2005[d]
            rightText = ""

            # diff[(dB < 2.0) | (diff > 0.5)] = np.nan
            # diffStd[(dBStd < 2.0) | (diffStd > 0.5) ] = np.nan


            diffStd[(dBStd < 1.0)] = np.nan



            ax = fig.add_subplot(2, 1, indexes[index])


            ptsLables = []
            txtLables = []

            #dataXM = []
            #dataYM = []

            sliceSpec = np.nanmean(diff, axis =0)
            sliceSpecStd = diffStd

            #posIndexes = np.arange(0,1946)
            posIndexes = (dBStd > 1.0 )& (sliceSpec > -1000)









            dataXM = mcwdSlice[posIndexes]
            dataYM = sliceSpec[posIndexes]*10
            dataYMSTD = sliceSpecStd[posIndexes]*10






            np.savetxt(setup.folderPath + "\\" + datasetNames[d] + "-Data-Median_2010.tsv",  np.squeeze(np.stack((dataXM, dataYM), axis=0)), delimiter="\t")



            #r_squared = r2_score(self.f(px, a, b, c), py)
            #r_squared = mean_squared_error(self.f(px, a, b, c), py)

            #print(dP)
            pxx = np.linspace(-230, 50, 100)


            #samples = np.squeeze(np.stack((dataXM, dataYM), axis=0))
            #densObj = kde(samples)
            #colours = self.makeColours(densObj.evaluate(samples))

            popt, pcov = curve_fit(self.f, dataXM, dataYMSTD, p0=np.asarray([-0, -0.05, -0]), maxfev=100000)
            a = popt[0]
            b = popt[1]
            c = popt[2]

            pyy_std = self.f(pxx, *popt)

            std_RMSE= mean_squared_error(self.f(-phillips_points[:, 0], a, b, c),  phillips_points[:, 1] * 0.5, squared=False)


            #impacts[0 + 2*d ]= dataXM
            #impacts[1 + 2*d ]= dataYM

            popt, pcov = curve_fit(self.f, dataXM, dataYM, p0=np.asarray([-0.01, -0.03, -2]), maxfev=100000)
            a = popt[0]
            b = popt[1]
            c = popt[2]
            residuals = dataYM - self.f(dataXM, a, b, c)
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((dataYM - np.mean(dataYM)) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            dP = np.sqrt(np.diag(pcov))

            hyd_RMSE= mean_squared_error(self.f(-phillips_points[:, 0], a, b, c),  phillips_points[:, 1] * 0.5, squared=False)

            #lpb, upb = self.predband(pxx,dataXM, dataYM, *popt, 0.95)
            lpb, upb = self.predband(pxx, dataXM, dataYM, popt, self.f, conf=0.95)





            pyy = self.f(pxx, *popt)






            ptsLable =  ax.scatter(dataXM, dataYM, color='tab:blue', alpha = 0.7, s=0.5)
            ptsLable2 =  ax.scatter(dataXM, dataYMSTD, color='tab:purple', alpha = 0.3, s=0.5)
            ptsLables.append(ptsLable)
            ptsLables.append(ptsLable2)



            #meansAndBandsSlice = meanAndBands[meanAndBands[:,0] == datasetNames[d]]

            ax.plot(pxx, pyy, linewidth=1.5, alpha=1, color = "tab:blue")
            plt.plot(pxx, lpb, 'k--', linewidth=0.75, alpha=1, color = "tab:blue")
            plt.plot(pxx, upb, 'k--', linewidth=0.75, alpha=1, color = "tab:blue")
            ax.plot(pxx, pyy_std, linewidth=1.0, alpha=0.6, color = "tab:purple")
            #ax.plot(meansAndBandsSlice[:,1], meansAndBandsSlice[:,3], linewidth=1.5, alpha=1, color = "gray" ,linestyle='--')
            #ax.plot(meansAndBandsSlice[:,1], meansAndBandsSlice[:,4], linewidth=1.5, alpha=1, color = "gray" ,linestyle='--')

            #dataY.append(pyy)

            if i == 37:
                lineStd, = ax.plot(pxx, pyy, linewidth=1.5, color='tab:purple', alpha=1.0)
            rsquares.append(rsquares_dat)


            # lineQuantileInnerHQ, = ax.plot(pxx, hqInner, linewidth=1.2, color='orangered', linestyle=(0, (2, 2)))
            # lineQuantileInnerLQ, = ax.plot(pxx, lqInner, linewidth=1.2, color='orangered', linestyle=(0, (2, 2)))

            # line0, = ax.plot(pxx, np.zeros(len(pxx)), linewidth=1.2, color='silver')

            ax.set_xlabel(r'MCWD anomaly [mm $year^{-1}$]')
            ax.set_ylabel(r'AGB anomaly [MgC $ha^{-1}$  $year^{-1}$]')
            # ax.text(0.97, 0.06, '${R}^2=$' + str(np.round(r_squared, 2)), horizontalalignment='right', transform=ax.transAxes,
            #        size=12,
            #        bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

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
            ax.set_ylim([-20, 5])

            xStd = np.linspace(-120, 40, 10)
            xExt = np.linspace(-240, -120, 10)

            #line4, = ax.plot(xStd, -self.p2009v(xStd), linewidth=1.5, color='dodgerblue')
            #line5, = ax.plot(xExt, -self.p2009v(xExt), linestyle=(0, (2, 2)), linewidth=1.5, color='dodgerblue')
            #line6, = ax.plot(xExt, -self.p2009_constv(xExt), linestyle=(0, (1, 1)), linewidth=1.5, color='dodgerblue')

            #pobns = ax.scatter(-phillips_points[:, 0], phillips_points[:, 1] , color='black', s=10)
            #13.04.2021 fix from Mg/m² tö MgC/m²
            pobns = ax.scatter(-phillips_points[:, 0], phillips_points[:, 1]*0.5 , color='black', s=10)
            index += 1









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
            indexesOI = [34, 0, 30, 10, 27, 36]
            indexesOI = [34, 30, 10, 36]
            indexesOI = [34, 30, 0, 10,  27,]
            indexesOI = np.arange(0,37)

            ptsLables = []
            txtLables = []
            for i in indexesOI:
                if i != 37:
                    sliceSpec = diff[i, :]
                else:
                    sliceSpec = diffStd


                sliceSpecd = sliceSpec
                mcwdSliced = mcwdSlice

                #mcwdSliced = mcwdSlice
                #sliceSpecd = sliceSpec

                nonNanIndexes = (np.logical_not(np.isnan(mcwdSliced))) \
                                & (np.logical_not(np.isnan(sliceSpecd)))


                py = sliceSpecd[nonNanIndexes] *10
                px = mcwdSliced[nonNanIndexes]



                #if i == indexesOI[0]:
                 #   dataIndividual[0 + 37*d] = px




                if py.shape[0] > 470:

                    #dataIndividual[i+1 + 37*d] = py

                    popt, pcov = curve_fit(self.f, px, py, p0=np.asarray([-0.04, -0.23, +2]), maxfev=100000)
                    a = popt[0]
                    b = popt[1]
                    c = popt[2]
                    residuals = py - self.f(px, a, b, c)
                    ss_res = np.sum(residuals ** 2)
                    ss_tot = np.sum((py - np.mean(py)) ** 2)
                    r_squared = 1 - (ss_res / ss_tot)
                    dP = np.sqrt(np.diag(pcov))

                    # r_squared = r2_score(self.f(px, a, b, c), py)
                    # r_squared = mean_squared_error(self.f(px, a, b, c), py)

                    # print(dP)

                    paramData.append([datasetNames[d], rightText, a, b, c])

                    pxx = np.linspace(-250, 50, 100)
                    pyy = self.f(pxx, *popt)
                    # lpb, upb = predband(pxx, px, py, popt, f, conf=0.95)

                    rsquares_dat.append(r_squared)
                    # if (r_squared > 0.1) & (i != 37):
                    if ((i != 37)):
                        ax.scatter(px, py, s=1, alpha=0.0, edgecolors=None)

                        #ptsLable = ax.plot(pxx, pyy, linewidth=1, alpha=1, c = cmapPsi(normPsi(inputData[i, 0])))
                        ptsLable = ax.plot(pxx, pyy, linewidth=1, alpha=1, c = viridis(normPsi(inputData[i, 0])))
                        ptsLables.append(ptsLable)
                        dataY.append(pyy)

                    #if i == 37:
                    #    lineStd, = ax.plot(pxx, pyy, linewidth=1.5, color='tab:gray', alpha=1.0)


            rsquares.append(rsquares_dat)


            # lineQuantileInnerHQ, = ax.plot(pxx, hqInner, linewidth=1.2, color='orangered', linestyle=(0, (2, 2)))
            # lineQuantileInnerLQ, = ax.plot(pxx, lqInner, linewidth=1.2, color='orangered', linestyle=(0, (2, 2)))

            # line0, = ax.plot(pxx, np.zeros(len(pxx)), linewidth=1.2, color='silver')

            ax.set_xlabel(r'$\Delta$MCWD [mm $year^{-1}$]')
            ax.set_ylabel(r'$\Delta$AGB [MgC $ha^{-1}$  $year^{-1}$]')
            # ax.text(0.97, 0.06, '${R}^2=$' + str(np.round(r_squared, 2)), horizontalalignment='right', transform=ax.transAxes,
            #        size=12,
            #        bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

            #ax.text(0.97, 0.05, datasetNames[d], horizontalalignment='right', verticalalignment='center',
            #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            #ax.text(0.03, 0.93, lowerletters[indexes[index]-1] + ")  $\psi_{50} > -2.4$MPa", horizontalalignment='left',
            ax.text(0.03, 0.94, lowerletters[indexes[index]-1] + ")  Individual", horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.text(0.15, 0.95, rightText, horizontalalignment='left', verticalalignment='center',
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            ax.set_xlim([-250, 50])
            ax.set_ylim([-20, 5])

            xStd = np.linspace(-120, 40, 10)
            xExt = np.linspace(-240, -120, 10)


            #line4, = ax.plot(xStd, -self.p2009v(xStd), linewidth=1.5, color='dodgerblue')
            #line5, = ax.plot(xExt, -self.p2009v(xExt), linestyle=(0, (2, 2)), linewidth=1.5, color='dodgerblue')
            #line6, = ax.plot(xExt, -self.p2009_constv(xExt), linestyle=(0, (1, 1)), linewidth=1.5, color='dodgerblue')

            #13.04.2021 fix from Mg/m² tö MgC/m²
            #pobns = ax.scatter(-phillips_points[:, 0], phillips_points[:, 1]*0.5 , color='black', s=10)
            index += 1



        j = 0
        for i in indexesOI:
            lablesTxt = r"$\psi_{50} = " + str(np.round(inputData[i, 0], 1)) + "$ MPa\n" + \
                        r"$\psi_{88} = " + str(np.round(inputData[i, 1], 1)) + "$ MPa"#\
                        #+ r"; $\tilde{R}_{GLDAS}^2 = " + str(
                #np.round(rsquares[0][j], 2)) + "$" +r"; $\tilde{R}_{WATCH}^2 = " + str(
                #np.round(rsquares[1][j], 2)) + "$"
            txtLables.append(lablesTxt)
            j +=1


        ptsLables.append(pobns)
        txtLables.append("Obs $(n=22)$")

        from matplotlib.lines import Line2D

        le_pts = Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label='Obs $(n=22)$', markersize=1.25)



        le_hyd = [Line2D([0], [0],marker='o' ,color='w',markerfacecolor='tab:blue', label='Modelled', markersize=0.75),
                           Line2D([0], [0], color='tab:blue', label='Exp-Fit', markersize=15),
                           Line2D([0], [0], color='tab:blue', label='95%-CI-bands', ls='--',
                                  markersize=15)
                           ]

        le_std = [         Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:purple',label='Modelled', markersize=0.75),
                           Line2D([0], [0], color='tab:purple', label='Exp-Fit', markersize=15)
                           ]

        legend1 = fig.legend(handles= [le_pts],
                  loc='center left',
                  markerscale=5,
                  bbox_to_anchor=(0.455, 0.85),
                   fontsize=8,
                   ncol=1,
                   fancybox=False, shadow=False)


        legend2 = fig.legend(handles=le_hyd,
                  loc='center left',
                  title="$\\bf{LPJ-GUESS-HYD}$",
                  markerscale=5,
                  bbox_to_anchor=(0.455, 0.75),
                   fontsize=8,
                   ncol=1,
                   fancybox=False, shadow=False)
        legend2.get_title().set_fontsize('8')  # legend 'Title' fontsize


        legend3 = fig.legend(handles=le_std,
                  loc='center left',
                  title="$\\bf{LPJ-GUESS}$",
                  markerscale=5,
                  bbox_to_anchor=(0.455, 0.6),
                   fontsize=8,
                   ncol=1,
                   fancybox=False, shadow=False)

        legend3.get_title().set_fontsize('8')  # legend 'Title' fontsize

        #for i in range(0, len(legend.legendHandles)-1):
       #     legend.legendHandles[i]
       # legend.legendHandles[len(legend.legendHandles)-1].set_sizes([15.0])


        #for handle in legend2.legendHandles:
            #handle.set_sizes([15.0])


        plt.subplots_adjust(right=0.45, wspace=0.3, hspace=0.25)

        c_map_ax = fig.add_axes([0.47, 0.111, 0.01, 0.34])
        c_map_ax.axes.get_xaxis().set_visible(True)
        c_map_ax.axes.get_yaxis().set_visible(True)


        bar = mpl.colorbar.ColorbarBase(c_map_ax,cmap= viridis, norm = normPsi, orientation="vertical", label = r'$\psi_{50}$ [MPa]')
        #bar.set_label('$\Delta$AGB [MgC/ha]', fontsize=12)



        df = pd.DataFrame(paramData)

        import os
        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)
        plt.savefig(setup.folderPath + "\\" + "Figure2_single-" + str(yearEnd) + ".png", dpi=600, bbox_inches='tight',
                    pad_inches=0)

        plt.close()



        df.to_csv(setup.folderPath + "\\ParametersFit" + str(yearEnd) + ".tsv")
        #np.savetxt(setup.folderPath + "\\Fig2-Data.tsv", impacts, delimiter="\t")
        #np.savetxt(setup.folderPath + "\\Individual2005ImpactData.tsv", dataIndividual, delimiter="\t")




        self.Succesfull()


