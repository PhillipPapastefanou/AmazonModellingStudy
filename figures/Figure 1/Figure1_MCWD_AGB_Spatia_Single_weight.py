from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from A03_Hydraulics_Implementation.AutoAnalysis.Source.SetupFiles import SetupFile
from A02_MCWD_Dataset_Analysis.Pylibs.MCWDevaluation import MCWDFile as MF
from Libs.Standard import StandardLPJG
from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis

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



class Figure1_Single_weight(BaseAnalysis):
    pass

    def __init__(self, setup):

        self.StartSetup()
        self.GetAGBFractionCMassTotal = np.vectorize(self.getAGBFractionCMassTotal)

        if platform == "darwin":
            raisg_mask = r"../AB-SHAPE/Amazon_shape.shp"

        mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                              ccrs.PlateCarree())
        standardLPJs = []
        for file in setup.stdFilePaths:
            standardLPJs.append(StandardLPJG(file, True))

        mcwdFiles = []
        for file in setup.mcwdPaths:
            mcwdFiles.append(StandardLPJG(file, True))

        data_mcwd_2005 = np.zeros((len(mcwdFiles), 1946))
        data_mcwd_2010 = np.zeros((len(mcwdFiles), 1946))

        MCWDFiles = []
        for file in setup.mcwdPaths:
            MCWDFiles.append(MF(file))

        i = 0
        for MCWDFile in MCWDFiles:
            MCWDFile.Parse(2000, 2010, [2005, 2010])
            data_mcwd_2005[i] = MCWDFile.DataSlices[0]
            data_mcwd_2010[i] = MCWDFile.DataSlices[1]
            i += 1

        # Save MCWD anomaly for figure 2
        lons = MCWDFiles[0].GeoFile.lons
        lats = MCWDFiles[0].GeoFile.lats
        arr = np.transpose(np.stack((lons, lats, data_mcwd_2005[0])))
        np.savetxt(setup.folderPath + "/"+"MCWD-" + setup.hydFileNames[0]+ ".csv", arr, delimiter=',')
        arr = np.transpose(np.stack((lons, lats, data_mcwd_2005[1])))
        np.savetxt(setup.folderPath + "/"+"MCWD-" + setup.hydFileNames[1]+ ".csv", arr, delimiter=',')


        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))


        normAGB = mpl.colors.Normalize(-15, 5)
        colorsAGB = [[normAGB(-15), "#9C0200"],
                     [normAGB(-5), "#ff9a00"],
                     [normAGB(-0.5), "white"],
                     [normAGB(0.0), "white"],
                     [normAGB(0.5), "white"],
                     [normAGB(5), "tab:blue"]]
        cmapAGB = mpl.colors.LinearSegmentedColormap.from_list("", colorsAGB)

        norm = mpl.colors.Normalize(-180, 0.0)
        colorsMCWD = [[norm(-180.0), (0.65, 0.16, 0.16)],
                      [norm(-100.0), (0.80, 0.58, 0.047)],
                      [norm(-25.0), "white"],
                      [norm(0.0), "white"], ]
        cmapMCWD = mpl.colors.LinearSegmentedColormap.from_list("", colorsMCWD)

        yearbegin = 1985
        yearEnd = 2005

        datasetNames = setup.hydFileNames

        fig = plt.figure(figsize=(9, 10))
        index = 1

        indexes = [1, 2, 3]

        import string
        lowerletters = string.ascii_lowercase[:26]

        for d in range(0, 1):
            biomass = listAnalysis[d].GetDataSlicesCondContinous("cmasstotal", yearbegin, yearEnd, 0)
            standardLPJs[d].GetBiomassBeforeAfter(yearbegin, yearEnd, True)
            dBStd = standardLPJs[d].BiomassBefore
            dAStd = standardLPJs[d].BiomassAfter

            fracAGB = self.GetAGBFractionCMassTotal(biomass)
            fracAGBstd = self.GetAGBFractionCMassTotal(dBStd)

            biomassAGB = biomass * fracAGB
            biomassAGBDiff = np.diff(biomassAGB, axis = 2)

            # Calculate AGB anomaly
            diff = biomassAGBDiff[:, :, 2005 - yearbegin - 1] - np.mean(biomassAGBDiff[:, :, 0:2004 - yearbegin - 2],
                                                                        axis=2)

            standardLPJs[d].GetBiomassBeforeAfter(yearbegin, yearEnd, True)
            dBStd = standardLPJs[d].BiomassBefore
            dAStd = standardLPJs[d].BiomassAfter

            fracAGBStd = self.GetAGBFractionCMassTotal(dBStd)
            diffStd = dAStd - dBStd
            diffStd *= fracAGBStd

            #Conversion from kgC/m² tp MgC/ha
            diffStd *= 10.0

            #Add a very tiny amount to the biomass to avoid division by zero
            w = biomassAGB[:, :, 0]  + 0.0001

            #Conversion from kgC/m² tp MgC/ha
            sliceSpec = np.average(diff, axis = 0, weights=w) * 0.0

            for i in range(0, 3):
                if i == 1:
                    #diffMean = np.nanmean(diff, axis=0)
                    diffMean = sliceSpec
                    rightText = "LPJ-GUESS-HYD"
                elif i == 2:
                    diffMean = diffStd
                    rightText = "LPJ-GUESS"
                else:
                    # diffMean = listAnalysis[d].GetMCWD(yearEnd)
                    if yearEnd == 2010:
                        diffMean = data_mcwd_2010[d]
                        rightText = ""
                    elif yearEnd == 2005:
                        diffMean = data_mcwd_2005[d]
                        rightText = ""
                    else:
                        diffMean = np.zeros(1946)
                        rightText = ""

                img = listAnalysis[d].CreateImage(diffMean)
                img_extent = listAnalysis[d].stratifiedOutput.IMG_extent
                offset = [-3, 3, -3, 3]

                axGeo = fig.add_subplot(3, 1, indexes[index - 1], projection=ccrs.PlateCarree())
                lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='g')
                lat_formatter = LatitudeFormatter()
                axGeo.xaxis.set_major_formatter(lon_formatter)
                axGeo.yaxis.set_major_formatter(lat_formatter)
                axGeo.add_feature(cfeature.BORDERS, edgecolor='tab:grey')
                axGeo.coastlines(resolution='110m', linewidth=1, color='tab:grey')
                axGeo.set_extent(list(np.array(img_extent) + np.array(offset)), crs=ccrs.PlateCarree())
                axGeo.add_feature(mask, edgecolor='black', linewidth=1.3, facecolor="None")
                axGeo.text(0.03, 0.92, lowerletters[indexes[index - 1] - 1] + ")", horizontalalignment='left',
                           verticalalignment='center',
                           transform=axGeo.transAxes, size=14,
                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                axGeo.text(0.15, 0.92, rightText, horizontalalignment='left', verticalalignment='center',
                           transform=axGeo.transAxes, size=10,
                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

                axGeo.set_xticks([-80, -70, -60, -50], crs=ccrs.PlateCarree())
                axGeo.set_yticks([-20, -15, -10, -5, 0, 5], crs=ccrs.PlateCarree())
                axGeo.set_xlabel(r'Longitude')
                axGeo.set_ylabel(r'Latitude')

                if (i == 1) | (i == 2):
                    imsh1 = axGeo.imshow(img, transform=ccrs.PlateCarree(), extent=img_extent, cmap=cmapAGB,
                                         norm=normAGB)
                else:
                    imsh2 = axGeo.imshow(img, transform=ccrs.PlateCarree(), extent=img_extent, cmap=cmapMCWD, norm=norm)
                index += 1

        plt.subplots_adjust(right=0.90, wspace=0.1, hspace=0.3)

        cax = plt.axes([0.67, 0.666, 0.02, 0.213])

        bar = plt.colorbar(imsh2, cax=cax, orientation="vertical")
        bar.set_label('$\Delta$MCWD [mm]', fontsize=12)

        cax = plt.axes([0.67, 0.111, 0.02, 0.491])
        bar = plt.colorbar(imsh1, cax=cax, orientation="vertical")
        bar.set_label('$\Delta$AGB [MgC/ha]', fontsize=12)

        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)

        plt.savefig(setup.folderPath + "/" + "Figure_1_SINGLE_Geospatial-MCWD-AGB-" + str(yearEnd) + "_weight.png", dpi=300, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        self.Succesfull()



