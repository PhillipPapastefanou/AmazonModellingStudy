import matplotlib.pyplot as plt
import numpy as np
import  pandas as pd
import scipy.stats as stats

from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis
import matplotlib as mpl
from matplotlib import gridspec



class Figure3_RatesWithStd_Single_ArrowBox_ABCOORD_smooth_w(BaseAnalysis):

    def moving_average(self, a, n):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    def __init__(self, setup):

        self.StartSetup()

        hubauFiles = ["..//Source//Hubau_AGB_NETCHANGE.csv", "..//Source//Hubau_AGB_Gains.csv", "..//Source//Hubau_AGB_Morts.csv"]
        modelledFiles = [setup.folderPath+"/wNET_SINK_SAP_HEART_OverTime_individual.tsv",
                         setup.folderPath+"/wNPP_SAP_HEART_OverTime_individual.tsv",
                         setup.folderPath+"/wMORT_SAP_HEART_OverTime_individual.tsv",
                         setup.folderPath+"/wNET_SINK_SAP_HEART_OverTime_individual_BC.tsv",
                         setup.folderPath + "/wNPP_SAP_HEART_OverTime_individual_BC.tsv",
                         setup.folderPath+"/wMORT_SAP_HEART_OverTime_individual_BC.tsv"]

        hubauData = []
        hydraulicData = []

        for file in hubauFiles:
            hubauData.append(pd.read_csv(file, header=0).values)

        for file in modelledFiles:
            hydraulicData.append(pd.read_csv(file, header=0, sep="\t").values)

        fig = plt.figure(figsize=(8, 10))

        spec = gridspec.GridSpec(ncols=2, nrows=3,
                                 width_ratios=[8, 1], wspace=-0.05)

        index = 1
        cols = ["#3f71f5","#228B22"]
        yLabs = ["Net carbon\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon gains\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon lossses\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Net carbon sink\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon gains\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon lossses\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$"]

        pltRanges = [[-1.5, 1.5], [1.0, 3.5], [0.5, 4.5], [-1.5, 1.5], [1.0, 3.5], [0.5,4.5]]
        hubauPoints = [-0.016, 0.014, 0.023]


        lewisOffset = -0.3
        lewisMort = np.asarray([4.2, 3.0])
        lewisMortCILow = np.asarray([2.3, 1.5])
        lewisMortCIHigh =np.asarray([6.4,4.9])

        feldpauOffset = +0.1
        feldpauMortTimes =[2010+feldpauOffset]
        feldpauMort = np.asarray([1.95])
        feldpauMortCILow = feldpauMort - np.asarray([1.18])
        feldpauMortCIHigh =np.asarray([2.77]) -feldpauMort

        feldpauNetOffset = +0.0
        feldpauNetTimes =[2010+feldpauOffset]
        feldpauNet = np.asarray([-0.43])
        feldpauNetCILow = feldpauNet - np.asarray([0.18])
        feldpauNetCIHigh =np.asarray([-1.13]) -feldpauNet

        lines = []

        import string
        lowerletters = string.ascii_lowercase[:26]

        self.Errors = []

        for i in range(0, 3):

            if i < 3:
                rightText = r"Amazon basin $\left(n=1946\right)$"
            else:
                rightText = r"Inventory plots $\left(n=47\right)$"


            ax = fig.add_subplot(spec[2*i])

            # only on the left side
            if i > -1:
                x = hubauData[i % 3][:, 0]
                y = hubauData[i % 3][:, 1]
                yL = hubauData[i % 3][:, 2]
                yH = hubauData[i % 3][:, 3]


                sIndex = np.where((1984 <= x) & (x <= 2010.9))
                yObs= y[sIndex]
                yObs = np.mean(yObs.reshape(-1, 10), axis=1)

                line1 = ax.plot(x, y, '--', color="black", alpha = 0.5)
                ax.fill_between(x, yL, yH, alpha=0.1, color="tab:gray")
                index += 1

            if (i == 0) | (i == 3):
                x = hubauData[i % 3][:, 0]
                ax.plot(x, np.zeros(x.shape[0]), '-', color="black", alpha=0.2, linewidth=1)

            cI = 0
            pos = [0.19, 0.07]

            if i == 0:
                ax.annotate('', xy=(-0.17, 0), xycoords='axes fraction', xytext=(-0.17, 0.47),
                            arrowprops=dict(arrowstyle="->", facecolor='black', lw=1))
                ax.annotate('', xy=(-0.17, 0.53), xycoords='axes fraction', xytext=(-0.17, 1),
                            arrowprops=dict(arrowstyle='<-', facecolor='black', lw=1))

                plt.text(0.015, 0.745, "Sink", fontsize=10, transform=plt.gcf().transFigure, rotation = 90)
                plt.text(0.015, 0.64, "Source", fontsize=10, transform=plt.gcf().transFigure, rotation=90)

            scentext = ["AB"]

            for d in ["GLDAS"]:

                for scen in range(0,1):

                    sliceIndex = hydraulicData[i+3*scen][:, 1] == d
                    sliceData = hydraulicData[i+ 3*scen][sliceIndex]
                    x = sliceData[:, 0].astype(np.int)

                    yBase = sliceData[:, 2:(2 + 37)].astype(np.float)
                    yBaseStd = sliceData[:, 2 + 37].astype(np.float) * 10.0

                    y = np.median(yBase, axis=1) * 10
                    yL = np.quantile(yBase, axis=1, q=0.8) * 10.0
                    yH = np.quantile(yBase, axis=1, q=0.2) * 10.0

                    xReP = np.repeat(x[np.newaxis, :], 37, axis=0)
                    sIndex = np.where((1984 <= x) & (x <= 2010))
                    xReP = np.transpose(np.squeeze(xReP[:, sIndex]))
                    yBase = np.squeeze(yBase[sIndex, :])
                    xReP = np.squeeze(xReP.reshape(-1, 1))
                    yBase = np.squeeze(yBase.reshape(-1, 1) * 10.0)

                    slope, intercept, r_value, p_value, std_err = stats.linregress(xReP, yBase)
                    slope_std, intercept_std, r_value_std, p_value_std, std_err_std = stats.linregress(x, yBaseStd)

                    yLinear = x * slope + intercept
                    yLinearStd = x * slope_std + intercept_std

                    y_roll = self.moving_average(y , 2)
                    yL_roll = self.moving_average(yL , 2)
                    yH_roll = self.moving_average(yH , 2)
                    x_offset = x[1:]

                    line = ax.plot(x, yLinear, color=cols[scen], alpha=1.0, lw = 2)
                    if  scen == 0:
                        line = ax.plot(x_offset, y_roll, color=cols[cI], ls="--", alpha=0.5)
                    lines.append(line)

                    if  scen == 0:
                        line3 = ax.plot(x, yLinearStd, color="tab:purple", alpha=0.8)
                        ax.fill_between(x_offset, yL_roll, yH_roll, alpha=0.15, color=cols[cI])


                    p_valuestr = str(np.round(slope, 3)) + " $"
                    if p_value < 0.001:
                        p_valuestr += "   (p-value < 0.001)"

                    else:
                        p_valuestr += "   (p-value = "
                        p_valuestr += str(np.round(p_value, 3)) + ")"

                    p_valuestr_std = str(np.round(slope_std, 3)) + " $"
                    if p_value2 < 0.001:
                        p_valuestr_std += "   (p-value < 0.001)"

                    else:
                        p_valuestr_std += "   (p-value = "
                        p_valuestr_std += str(np.round(p_value_std, 3)) + ")"

                    ax.text(0.98, pos[scen], r"$\mathrm{slope}_\mathrm{Hyd}\left(\mathrm{"+  scentext[scen] + r"}\right)= " + p_valuestr,
                            horizontalalignment='right', verticalalignment='center', c=cols[scen],
                            transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=1.0, edgecolor='white'))

                cI += 1

                if  i > 2:
                    sIndex = np.where((1984 <= x) & (x <= 2010))
                    yEstimate = y[sIndex]
                    self.Errors.append(np.sqrt(np.nanmean((yEstimate-yObs)**2)))


            ax.set_xlim((1984, 2010.5))
            ax.set_ylim(pltRanges[i])
            ax.set_ylabel(yLabs[i])
            ax.set_xlabel("Year")

            ax.text(0.02, 0.15, "$\mathrm{slope}_{Obs}= " +str(np.round(hubauPoints[i - 3], 4))+ "$", horizontalalignment='left', verticalalignment='bottom',
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

            ax.text(0.02, 0.04, r"$\mathrm{slope}_\mathrm{std}= " + p_valuestr_std, horizontalalignment='left',
                    verticalalignment='bottom', c="tab:purple",
                    transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=1.0, edgecolor='white'))

            if  i < 3:
                print(i)
                ax.text(0.02, 0.92, lowerletters[i] + ")", horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes, size=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))


            ax2 = fig.add_subplot(spec[2*i + 1])

            for d in ["GLDAS"]:
                sliceIndex = hydraulicData[i][:, 1] == d
                sliceData = hydraulicData[i][sliceIndex]
                x = sliceData[:, 0].astype(np.int)
                yBase = sliceData[:, 2:(2+37)].astype(np.float)

                sIndex = np.where((1984 <= x) & (x <= 2010))

                slopes = np.zeros(37)
                inters = np.zeros(37)
                for c in range(0,37):
                    slope, intercept, r_value, p_value, std_err = stats.linregress(x[sIndex], 10*yBase[sIndex,c])
                    slopes[c] = slope
                    inters[c] = intercept


            bp= ax2.boxplot(slopes *2010+ inters,  patch_artist=True)

            for patch in bp['boxes']:
                patch.set(facecolor='tab:blue')
            ax2.set_ylim(pltRanges[i])

            ax2.axis('off')


        from matplotlib.pyplot import Line2D

        le_hyd = [
            Line2D([0], [0], color='black', label='Obs', ls='--', markersize=15),
            Line2D([0], [0], color='tab:blue', label='LPJ-GUESS-HYD mean', ls='--', markersize=15),
            Line2D([0], [0], color='tab:blue', label='LPJ-GUESS-HYD quantiles', linewidth = 15, alpha = 0.15),
            Line2D([0], [0], color='tab:blue', label='LPJ-GUESS-HYD mean trend', markersize=15),
            Line2D([0], [0], color='tab:purple', label='LPJ-GUESS mean trend', markersize=15)
                           ]

        l1 = fig.legend(handles= le_hyd,
                  loc='lower center',
                  markerscale=5,
                   bbox_to_anchor=(0.5, 0.03),
                   fontsize=10,
                   ncol=3,
                   fancybox=False, shadow=False)

        plt.subplots_adjust(bottom=0.15, wspace=0.24, hspace=0.25)

        import os
        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)

        plt.savefig(setup.folderPath+ "//Figure3_STd_single_present_clim_rolltwo_weight.png", dpi=300, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        self.Succesfull()





