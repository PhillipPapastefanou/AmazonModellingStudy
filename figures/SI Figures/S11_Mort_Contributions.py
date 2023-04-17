import matplotlib.pyplot as plt
import numpy as np
import  pandas as pd

from sklearn.linear_model import LinearRegression

from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis

class IndivMortSlopes(BaseAnalysis):

    def __init__(self, setup):

        self.StartSetup()

        modelledFiles = [
                    setup.folderPath + "/MORT_SAP_HEART_Greff_OverTime_individual.tsv",
                    setup.folderPath + "/MORT_SAP_HEART_BG_OverTime_individual.tsv",
                    setup.folderPath + "/MORT_SAP_HEART_CAV_OverTime_individual.tsv"]

        hubauData = []
        hydraulicData = []

        for file in modelledFiles:
            hydraulicData.append(pd.read_csv(file, header=0, sep="\t").values)

        datasetNames = ["GLDAS"]

        fig = plt.figure(figsize=(8, 6))
        index = 1
        yLabs = ["Net carbon sink\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon gains\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon lossses\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Net carbon sink\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon gains\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$",
                 "Carbon lossses\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$"]

        cols = ["tab:blue", "tab:green", "tab:red"]

        lines = []

        import string
        lowerletters = string.ascii_lowercase[:26]
        import seaborn as sns

        sns.set(style="whitegrid")

        cI = 1
        for d in datasetNames:

            sliceIndex = hydraulicData[0][:, 1] == d
            sliceData = hydraulicData[0][sliceIndex]
            x = sliceData[:, 0].astype(np.int)
            y0 = np.zeros(x.shape[0])

            ax = fig.add_subplot(1, 1, cI)

            avgs = np.zeros((3))
            for i in range(0, 3):
                sliceIndex = hydraulicData[i][:, 1] == d
                sliceData = hydraulicData[i][sliceIndex]
                x = sliceData[:, 0].astype(np.int)
                yBase = sliceData[:, 2:].astype(np.float)
                plt.subplots_adjust(bottom=0.15, top=0.8, left=0.15)
                y = np.median(yBase, axis=1) * 10

                yL = np.quantile(yBase, axis=1, q=0.8) * 10
                yH = np.quantile(yBase, axis=1, q=0.2) * 10

                line = ax.plot(x, y + y0, '--', color=cols[i])
                lines.append(line)
                ax.fill_between(x, yL + y0, yH + y0, alpha=0.5, color=cols[i])
                ax.fill_between(x, yH + y0, y0, alpha=0.2, color=cols[i])

                y0 += y
                avgs[i]  = np.mean(y)

            print("Average mortalities contributions in %")
            print(np.round(avgs/ np.sum(avgs),2) * 100)

            ax.set_xlim((1985, 2010))
            ax.set_ylim([0, 4])
            ax.set_ylabel("Accumulated carbon lossses\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-1}$")
            ax.set_xlabel("Year")

            ax.text(0.97, 0.92, d, horizontalalignment='right', verticalalignment='center',
                    transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            #ax.text(0.02, 0.92, lowerletters[cI - 1] + ")", horizontalalignment='left',
            #        verticalalignment='center',
            #        transform=ax.transAxes, size=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            # ax.text(0.02, 0.07, rightText, horizontalalignment='left', verticalalignment='center',
            #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

            cI += 1

        for d in datasetNames:

            sliceIndex = hydraulicData[0][:, 1] == d
            sliceData = hydraulicData[0][sliceIndex]
            x = sliceData[:, 0].astype(np.int)


            params = np.zeros((3, 37))

            for i in range(0, 3):
                sliceIndex = (hydraulicData[i][:, 1] == d) & (1984 <= hydraulicData[i][:, 0])
                sliceData = hydraulicData[i][sliceIndex]
                x = sliceData[:, 0].astype(np.int)
                yBase = sliceData[:, 2:].astype(np.float)
                #plt.subplots_adjust(bottom=0.15, top=0.8, left=0.15)

                for ds in range(0, 37):
                    y = yBase[:, ds] * 10
                    #print(f" {i} {yBase.shape}")
                    model = LinearRegression().fit(x.reshape((-1,1)), y)

                    params[i, ds] = model.coef_






        # for d in datasetNames:
        #
        #     sliceIndex = hydraulicData[0][:, 1] == d
        #     sliceData = hydraulicData[0][sliceIndex]
        #     x = sliceData[:, 0].astype(np.int)
        #
        #     ax = fig.add_subplot(3, 2, cI)
        #
        #     params = np.zeros((3, 37))
        #
        #     for i in range(0, 3):
        #         sliceIndex = (hydraulicData[i][:, 1] == d) & (1984 <= hydraulicData[i][:, 0])
        #         sliceData = hydraulicData[i][sliceIndex]
        #         x = sliceData[:, 0].astype(np.int)
        #         yBase = sliceData[:, 2:].astype(np.float)
        #         plt.subplots_adjust(bottom=0.15, top=0.8, left=0.15)
        #
        #         for ds in range(0, 37):
        #             y = yBase[:, ds] * 10
        #             print(f" {i} {yBase.shape}")
        #             model = LinearRegression().fit(x.reshape((-1,1)), y)
        #
        #             params[i, ds] = model.coef_
        #
        #     ax = sns.violinplot(data=np.transpose(params), palette=cols, inner='quartile', scale="width", alpha=0.1)
        #     ax.set(ylim=(-0.01, 0.03))
        #     ax.set_xticklabels(["$\mathrm{Mort}_\mathrm{Growth efficiency}$", "$\mathrm{Mort}_\mathrm{Background}$",
        #                         "$\mathrm{Mort}_\mathrm{Hydraulic failure}$"],
        #                        size=12, fontname="Arial")
        #     # ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], size=11)
        #     ax.set_ylabel("Carbon losses slope\nMg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-2}$ ", size=12)
        #     ax.text(0.97, 0.92, d, horizontalalignment='right', verticalalignment='center',
        #             transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
        #     ax.text(0.02, 0.92, lowerletters[cI - 1] + ")", horizontalalignment='left',
        #             verticalalignment='center',
        #             transform=ax.transAxes, size=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
        #
        #     # ax.text(0.02, 0.07, rightText, horizontalalignment='left', verticalalignment='center',
        #     #        transform=ax.transAxes, size=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
        #
        #     cI += 1

        # fig.legend(handles = [line1[0],lines[0][0],lines[1][0]] ,
        #           labels=['Hubau et al. 2020',
        #                   '$\mathrm{LPJ-GUESS-HYD}_\mathrm{GLDAS}$',
        #                   '$\mathrm{LPJ-GUESS-HYD}_\mathrm{WATCH}$'],
        #           loc='lower center',
        #            markerscale=50,
        #            bbox_to_anchor=(0.45,0.03),
        #            fontsize=12,
        #           ncol=3,
        #           fancybox=False, shadow=False)

        #plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.25)

        plt.savefig(setup.folderPath +"/IndivMortalityContrib.png", dpi=300, bbox_inches='tight',
                    pad_inches=0)

