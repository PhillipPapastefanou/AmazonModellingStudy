import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
import  pandas as pd
import scipy.stats as stats

from scipy.stats import gaussian_kde as kde
from matplotlib.colors import Normalize
from matplotlib import cm

from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis
from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified

class SI_Figure4_MCWD(BaseAnalysis):

    def makeColours(self, vals):
        colours = np.zeros((len(vals), 3))
        norm = Normalize(vmin=vals.min(), vmax=vals.max())

        # Can put any colormap you like here.
        colours = [cm.ScalarMappable(norm=norm, cmap='jet').to_rgba(val) for val in vals]

        return colours

    def __init__(self, setup):


        self.StartSetup()
        self.GetAGBFractionCMassTotal = np.vectorize(self.getAGBFractionCMassTotal)
        modelledFiles = [setup.folderPath+"//NET_SINK_SAP_HEART_OverTime_individual.tsv",
                         setup.folderPath+"//NET_SINK_SAP_HEART_OverTime_individual_BC.tsv"]
        hydraulicData = []


        for file in modelledFiles:
            hydraulicData.append(pd.read_csv(file, header=0, sep="\t").values)

        datasetNames = setup.hydFileNames

        fig = plt.figure(figsize=(10, 8))
        index = 1
        cols = ["tab:blue", "tab:red"]
        yLabs = ["Net carbon sink\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon gains\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon lossses\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Net carbon sink\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon gains\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon lossses\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$"]


        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))

        pltRanges = [[-2.5, 2.5],  [-2.5,2.5]]

        yearbegin = 1985
        yearEnd = 2010



        co2_data= pd.read_csv(r"F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\co2_1764_2014_extended_observed.dat", delimiter= "  ").values


        indexesYears = (co2_data[:, 0] >= yearbegin) & (co2_data[:, 0] <= yearEnd)

        co2_data_sel = co2_data[indexesYears, 1]

        lines = []

        import string
        lowerletters = string.ascii_lowercase[:26]

        self.Errors = []

        textAdd = ["Modelled output", "Climate input"]

        rightText = r"Amazon basin $\left(n=1946\right)$"

        i = 1
        c = 0


        di = 0
        for d in ["GLDAS"]:

            #listAnalysis[0].GetDataSlicesCond("cmasstotal", yearbegin, yearEnd, 1.0, 0)
            #print(listAnalysis[0].stratifiedOutput.PatchNames())

            sliceIndex = hydraulicData[c][:, 1] == d
            sliceData = hydraulicData[c][sliceIndex]
            x = sliceData[:, 0].astype(np.int)
            yBase = sliceData[:, 2:].astype(np.float)


            # ax.text(0.97, 0.92, d, horizontalalignment='right', verticalalignment='center',
            #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            #ax.text(0.02, 0.95, lowerletters[i - 1] + ")  " + textAdd[0], horizontalalignment='left',
            #        verticalalignment='center',
            #        transform=ax.transAxes, size=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            #ax.text(0.02, 0.07, rightText, horizontalalignment='left', verticalalignment='center',
            #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            #ax.text(0.97, 0.95, datasetNames[c], horizontalalignment='right', verticalalignment='center',
             #       transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))



            dBf = listAnalysis[di].GetDataSlicesRange("cmasstotal", yearbegin-1, yearEnd-1, 0)

            fracsAGB = self.GetAGBFractionCMassTotal(dBf)


            #dB =listAnalysis[di].GetDataSlicesRange("npp_sap", yearbegin , yearEnd, 0)
            #dB2 = listAnalysis[di].GetDataSlicesRange("npp_heart", yearbegin, yearEnd, 0)
            #diff = dB+ dB2



            dB =listAnalysis[di].GetDataSlicesRange("npp_sap", yearbegin , yearEnd, 0)
            #dB2 = listAnalysis[di].GetDataSlicesRange("npp_heart", yearbegin, yearEnd, 0.001)
            diff = dB

            mort = listAnalysis[di].GetDataSlicesRange("cmass_loss_heart_bg", yearbegin, yearEnd, 0)
            mort += listAnalysis[di].GetDataSlicesRange("cmass_loss_heart_cav", yearbegin, yearEnd, 0)
            mort += listAnalysis[di].GetDataSlicesRange("cmass_loss_heart_greff", yearbegin, yearEnd, 0)
            mort += listAnalysis[di].GetDataSlicesRange("cmass_loss_sap_bg", yearbegin, yearEnd, 0)
            mort += listAnalysis[di].GetDataSlicesRange("cmass_loss_sap_cav", yearbegin, yearEnd, 0)
            mort += listAnalysis[di].GetDataSlicesRange("cmass_loss_sap_greff", yearbegin, yearEnd, 0)

            mort[dBf < 1.0] = np.nan
            #dBf[dBf < 1.0] = np.nan
            mort *= fracsAGB


            diff[dBf < 1.0] = np.nan
            #dBf[dBf < 1.0] = np.nan
            diff *= fracsAGB


            diff =  diff -mort

            sLabs = ["Carbon gains", "Carbon losses"]

            for s in range(0,2):

                if s == 0:
                    slice = diff
                    ylim = (0, 10)
                elif (s == 1):
                    slice = mort
                    ylim = (0, 10)
                else:
                    slice = diff-mort
                    ylim = (-4, 4)

                # for driverName in ( "vpd", "temp_air", "mcwd_monthly"):
                for driverName in ("mcwd", "mcwd anomaly"):

                    ax = fig.add_subplot(2, 2, i)

                    if driverName == "mcwd":
                        driverVariable = listAnalysis[c].GetDataPatch("mcwd_monthly", yearbegin, yearEnd, 0)
                    else:
                        driverVariable = listAnalysis[c].GetDataPatch("mcwd_monthly", yearbegin, yearEnd, 0)
                        mean = np.mean(driverVariable, axis=1)
                        for j in range(0, driverVariable.shape[1]):
                            driverVariable[:, j] = driverVariable[:, j] - mean
                            #driverVariable[:, j][driverVariable[:, j] > 100.0] = np.nan

                    #if (i == 4) | (i == 8)| (i == 12):
                    #    mean = np.mean(driverVariable, axis=1)
                    #    for j in range(0, driverVariable.shape[1]):
                    #        driverVariable[:, j] = driverVariable[:, j] - mean
                    #        driverVariable[:, j][driverVariable[:, j] > 25.0] = np.nan

                    driverVariableRep = np.repeat(driverVariable[np.newaxis, :, :], 37, axis=0)

                    # meanDriverVariable = np.quantile(driverVariable, 0.5, axis=1)

                    x1 = driverVariableRep.reshape(-1, 1)
                    y1 = slice.reshape(-1, 1) * 10.0

                    xT = x1[~np.isnan(y1)]
                    yT = y1[~np.isnan(y1)]

                    xTau = xT[~np.isnan(xT)]
                    yT = yT[~np.isnan(xT)]

                    yTf = yT[yT > 0.01]
                    xTauf = xTau[yT > 0.01]

                    from scipy import stats


                    #xTau = stats.zscore(xTau)

                    samples = np.squeeze(np.stack((xTauf, yTf), axis=0))
                    indexes = np.random.choice(xTauf.shape[0], 30000, replace=False)

                    rndSamples = samples[:, indexes]

                    slope, intercept, r_value, p_value, std_err = stats.linregress(xTauf, yTf)

                    print('Slope: ', slope, '\nPValue: ', p_value)

                    ax.text(0.02, 0.95, lowerletters[i - 1] + ") "+ driverName, horizontalalignment='left',
                            verticalalignment='center',
                            transform=ax.transAxes, size=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                    #ax.text(0.97, 0.95, driverName, horizontalalignment='right', verticalalignment='center',
                    #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                    ax.text(0.02, 0.07, "slope = "+ str(np.round(slope, 3)), horizontalalignment='left', verticalalignment='center',
                            transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                    densObj = kde(rndSamples)
                    colours = self.makeColours(densObj.evaluate(rndSamples))
                    ax.scatter(rndSamples[0], rndSamples[1], color=colours, s=0.5)
                    #ax.scatter(rndSamples[0], rndSamples[1], s=0.5)

                    xS = np.arange(np.min(rndSamples[0]), np.max(rndSamples[0]),(np.max(rndSamples[0]) -np.min(rndSamples[0]))/100.0)
                    yS = slope*xS + intercept

                    ax.plot(xS, yS,  linestyle="dashed", c = "black")

                    ax.set_ylim(ylim)

                    if (i == 1) | (i == 3)| (i == 9):
                        ax.set_ylabel(sLabs[s])


                    # ax.plot(x, meanDriverVariable)
                    print(i)

                    i += 1







            #ax.set_ylim((-3,3))
           # ax.set_ylabel(yLabs[i])
            # ax.text(0.97, 0.92, d, horizontalalignment='right', verticalalignment='center',
            #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

            #ax.text(0.02, 0.07, rightText, horizontalalignment='left', verticalalignment='center',
            #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))



            di+= 1



        print("Writing file.. ")


        plt.subplots_adjust(bottom=0.15, wspace=0.24, hspace=0.25)

        import os
        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)

        plt.savefig(setup.folderPath+ "//Gains_Losses_MCWD.png", dpi=150, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        self.Succesfull()





