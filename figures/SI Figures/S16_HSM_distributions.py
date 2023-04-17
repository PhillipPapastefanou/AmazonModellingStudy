from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from A03_Hydraulics_Implementation.AutoAnalysis.Source.SetupFiles import SetupFile
from A02_MCWD_Dataset_Analysis.Pylibs.MCWDevaluation import MCWDFile as MF
from Libs.Standard import StandardLPJG



import matplotlib.pyplot as plt
import matplotlib.colors
import os
from sys import platform
import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis
import pandas as pd


class HSM_Distributions(BaseAnalysis):
    pass

    def __init__(self, setup):

        self.StartSetup()
        self.GetAGBFractionCMassTotal = np.vectorize(self.getAGBFractionCMassTotal)

        if platform == "darwin":
            raisg_mask = r"/Users/pp/Dropbox/ClimateData/AmazonBasin/AB-SHAPE/Amazon_shape.shp"

            inputData = pd.read_csv('/Users/pp/Dropbox/UNI/Projekte/A03_Hydraulics_Implementation/Analysis/AllCavCurves05-2020.tsv',
                sep='\t',
                header=None).values

        else:
            raisg_mask = r"F:\Dropbox\ClimateData\AmazonBasin\AB-SHAPE\Amazon_shape.shp"

        mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                              ccrs.PlateCarree())

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


        pos = np.squeeze(np.argwhere(data2005[0] < -175))

        lons = MCWDFiles[0].GeoFile.lons[pos]
        lats = MCWDFiles[0].GeoFile.lats[pos]

        arr = np.transpose(np.stack((lons, lats)))
        np.savetxt(setup.folderPath + "/"+"Less175MM_MCWD-" + setup.hydFileNames[0]+ ".csv", arr, delimiter=',')


        lons = MCWDFiles[0].GeoFile.lons
        lats = MCWDFiles[0].GeoFile.lats
        arr = np.transpose(np.stack((lons, lats, data2005[0])))
        np.savetxt(setup.folderPath + "/"+"MCWD-" + setup.hydFileNames[0]+ ".csv", arr, delimiter=',')
        arr = np.transpose(np.stack((lons, lats, data2005[1])))
        np.savetxt(setup.folderPath + "/"+"MCWD-" + setup.hydFileNames[1]+ ".csv", arr, delimiter=',')


        listAnalysis = []
        listAnalysisMonth = []
        for i in range(0,1):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))
            listAnalysisMonth.append(AmazonBasinStratified(setup.hydFilePathsMonth[i], setup.hydFileParamsPaths[i]))


        normAGB = mpl.colors.Normalize(-15, 5)
        colorsAGB = [[normAGB(-15), "#9C0200"],
                    # [normAGB(-1.0), "#ff7400"],
                     [normAGB(-5), "#ff9a00"],
                     [normAGB(-0.5), "white"],
                     [normAGB(0.0), "white"],
                     [normAGB(0.5), "white"],
                     [normAGB(5), "tab:blue"]]
        cmapAGB = mpl.colors.LinearSegmentedColormap.from_list("", colorsAGB)

        norm = mpl.colors.Normalize(0.2, 0.9)
        colorsMCWD = [[norm(0.2), 'tab:red'],
                      [norm(0.9), 'tab:blue'] ]
        cmapMCWD = mpl.colors.LinearSegmentedColormap.from_list("", colorsMCWD)

        yearbegin = 1985
        yearEnd = 2010

        datasetNames = setup.hydFileNames

        fig = plt.figure(figsize=(9, 10))
        index = 1

        indexes = [1, 2, 3]

        import string
        lowerletters = string.ascii_lowercase[:26]

        for d in range(0, 1):
            biomass = listAnalysis[d].GetDataSlicesCondContinous("cmasstotal", yearbegin, yearEnd, 1.0, 0)
            fracAGB = self.GetAGBFractionCMassTotal(biomass)
            avg_biomassAGB = np.mean(biomass * fracAGB, axis = 2)


            psi_leaf = listAnalysisMonth[d].GetDataSlicesRange("psileaf", yearbegin, yearEnd)

            psi_leaf = psi_leaf.reshape(37,1946, -1, 12 )

            psi_leaf_min = np.min(psi_leaf, axis = 3)

            psi_leaf_min_mean = np.mean(psi_leaf_min, axis = 2)


            psi50 = inputData[:,0]



            arr = np.zeros(1946)
            for s in range(1946):
                hsm =  psi_leaf_min_mean[:, s] - psi50
                arr[s] =  np.sum(hsm*avg_biomassAGB[:,s])/np.sum(avg_biomassAGB[:,s])


            for i in range(0, 1):

                img = listAnalysis[d].CreateImage(arr)
                img_extent = listAnalysis[d].stratifiedOutput.IMG_extent
                offset = [-3, 3, -3, 3]

                axGeo = fig.add_subplot(1, 1, indexes[index - 1], projection=ccrs.PlateCarree())
                lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='g')
                lat_formatter = LatitudeFormatter()
                axGeo.xaxis.set_major_formatter(lon_formatter)
                axGeo.yaxis.set_major_formatter(lat_formatter)
                axGeo.add_feature(cfeature.BORDERS, edgecolor='tab:grey')
                axGeo.coastlines(resolution='110m', linewidth=1, color='tab:grey')
                # axGeo.set_title("Precipitation")
                axGeo.set_extent(list(np.array(img_extent) + np.array(offset)), crs=ccrs.PlateCarree())
                axGeo.add_feature(mask, edgecolor='black', linewidth=1.3, facecolor="None")
                # titleTxt = axGeo.set_title(titles[i], size=16, loc = 'left')
                #axGeo.text(0.97, 0.92, datasetNames[d], horizontalalignment='right', verticalalignment='center',
                #           transform=axGeo.transAxes, size=10,
                #           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                axGeo.text(0.03, 0.92, lowerletters[indexes[index - 1] - 1] + ")", horizontalalignment='left',
                           verticalalignment='center',
                           transform=axGeo.transAxes, size=24,
                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                #axGeo.text(0.15, 0.92, rightText, horizontalalignment='left', verticalalignment='center',
                #           transform=axGeo.transAxes, size=10,
                #           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                axGeo.set_xticks([-80, -70, -60, -50], crs=ccrs.PlateCarree())
                axGeo.set_yticks([-20, -15, -10, -5, 0, 5], crs=ccrs.PlateCarree())
                axGeo.set_xlabel(r'Longitude')
                axGeo.set_ylabel(r'Latitude')


                imsh2 = axGeo.imshow(img, transform=ccrs.PlateCarree(), extent=img_extent, cmap=cmapMCWD, norm=norm)
                index += 1

        plt.subplots_adjust(right=0.90, wspace=0.1, hspace=0.3)

        cax = plt.axes([0.95, 0.19, 0.03, 0.62])

        bar = plt.colorbar(imsh2, cax=cax, orientation="vertical")
        bar.set_label('Hydraulic safety margins $\psi_{50}$ [MPA]', fontsize=12)


        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)

        plt.savefig(setup.folderPath + "/HSM_distribution" + ".png", dpi=300, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        self.Succesfull()



