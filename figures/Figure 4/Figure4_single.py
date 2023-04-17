import matplotlib.pyplot as plt
import numpy as np
import  pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression

from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis

class Figure4_HydraulicFailureMortSingleExclude(BaseAnalysis):

    def __init__(self, setup):

        self.StartSetup()

        inputData = pd.read_csv(
            '../AllCavCurves05-2020.tsv',
            sep='\t',
            header=None).values

        modelledFiles = [
                    setup.folderPath + "/MORT_SAP_HEART_Greff_OverTime_individual.tsv",
                    setup.folderPath + "/MORT_SAP_HEART_BG_OverTime_individual.tsv",
                    setup.folderPath + "/MORT_SAP_HEART_CAV_OverTime_individual.tsv"]


        hydraulicData = []

        for file in modelledFiles:
            hydraulicData.append(pd.read_csv(file, header=0, sep="\t").values)

        texts = ["growth-related\n    mortality","age-related\n    mortality", "drought-related\n    mortality", "drought-related\n    mortality\n    no extreme years"]

        fig = plt.figure(figsize=(12, 10))
        cols = ["tab:blue",  "tab:red"]

        import string
        lowerletters = string.ascii_lowercase[:26]
        import seaborn as sns

        sns.set(style="whitegrid")

        cI = 1
        all_slopes = np.zeros((4,37))
        all_intercepts = np.zeros((4,37))
        all_pvalues = np.zeros((4,37))


        for index in range(0, 4):

            ax = fig.add_subplot(2, 2, cI)
            di = 0
            for d in ["GLDAS"]:

                i = index
                if  index == 3:
                    i  = 2

                sliceIndex = (hydraulicData[i][:, 1] == d) & (1984 <= hydraulicData[i][:, 0])
                sliceData = hydraulicData[i][sliceIndex]
                x = sliceData[:, 0].astype(np.int)
                yBase = sliceData[:, 2:].astype(np.float)


                if  index == 3:
                    indexesNoExtreme =  (x != 2005) & (x != 2007)& (x != 2009) & (x != 2010)
                    x = x[indexesNoExtreme]


                slopes = np.zeros(37)
                intercepts = np.zeros(37)
                pValues = np.zeros(37)



                for ds in range(0, yBase.shape[1]):

                    y = yBase[:, ds] * 10

                    if  index == 3:
                        y = y[indexesNoExtreme]


                    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                    slopes[ds] = slope
                    pValues[ds] = p_value
                    intercepts[ds] = intercept

                for ds in range(0, yBase.shape[1]):
                    if pValues[ds] < 0.1:
                        ax.scatter(inputData[ds, 0], slopes[ds], s = 40 , c=cols[di], edgecolors=cols[di])
                    else:
                        ax.scatter(inputData[ds, 0], slopes[ds], s = 40 , facecolors='none',  edgecolors=cols[di])


                all_slopes[index] = slopes
                all_intercepts[index] = intercepts
                all_pvalues[index] = pValues

                di += 1

            ax.set(ylim=(-0.01, 0.05))
            if  (cI == 1) | (cI == 3):
                ax.set_ylabel("Slope of carbon losses \n[Mg C $\mathrm{ha}^{-1}$ $\mathrm{yr}^{-2}$] ", size=14)
            ax.set_xlabel("$\psi_{50}$ [MPa]", size=14)

            ax.tick_params(axis='both', which='major', labelsize=14)
            ax.text(0.02, 0.97, lowerletters[cI - 1] + ") " + texts[cI - 1], horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes, size=16, bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            cI += 1


        plt.subplots_adjust(bottom=0.15, wspace=0.25, hspace=0.25)
        np.savetxt(setup.folderPath +"/Figure4_Slopes.tsv", all_slopes, delimiter='\t')
        np.savetxt(setup.folderPath +"/Figure4_Intercepts.tsv", all_intercepts, delimiter='\t')
        np.savetxt(setup.folderPath +"/Figure4_PValues.tsv", all_pvalues, delimiter='\t')
        plt.savefig(setup.folderPath +"/Figure4_PaperHydraulicFailure_Exclude.png", dpi=300, bbox_inches='tight',
                    pad_inches=0)

        self.Succesfull()

