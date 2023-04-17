import pandas as pd

from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from A03_Hydraulics_Implementation.AutoAnalysis.Source.SetupFiles import SetupFile
from A02_MCWD_Dataset_Analysis.Pylibs.MCWDevaluation import MCWDFile as MF
from Libs.Standard import StandardLPJG
from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis

import matplotlib.pyplot as plt
import matplotlib.colors
import os

import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from sys import platform


class Figure_AGB_CHANGE_OVERALL(BaseAnalysis):
    pass

    def __init__(self, setup):

        self.StartSetup()
        self.GetAGBFractionCMassTotal = np.vectorize(self.getAGBFractionCMassTotal)

        if platform == "darwin":
            raisg_mask = r"/Users/pp/Dropbox/ClimateData/AmazonBasin/AB-SHAPE/Amazon_shape.shp"

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

        data_hubau_brienen = pd.read_csv("/Users/pp/Dropbox/UNI/Projekte/A03_Hydraulics_Implementation/PyFigures/FigureHubauBrienen/BrienenAmazonpositionswithcoords.tsv", sep='\t', header = None)
        data_hubau_brienen.columns =  ['Mathematica_Index', 'lon', 'lat']
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
        np.savetxt(setup.folderPath + "\\"+"Less175MM_MCWD-" + setup.hydFileNames[0]+ ".csv", arr, delimiter=',')


        lons = MCWDFiles[0].GeoFile.lons
        lats = MCWDFiles[0].GeoFile.lats
        arr = np.transpose(np.stack((lons, lats, data2005[0])))
        np.savetxt(setup.folderPath + "\\"+"MCWD-" + setup.hydFileNames[0]+ ".csv", arr, delimiter=',')
        arr = np.transpose(np.stack((lons, lats, data2005[1])))
        np.savetxt(setup.folderPath + "\\"+"MCWD-" + setup.hydFileNames[1]+ ".csv", arr, delimiter=',')


        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))


        normAGB = mpl.colors.Normalize(-10, 10)
        colorsAGB = [[normAGB(-10), "#9C0200"],
                    # [normAGB(-1.0), "#ff7400"],
                     [normAGB(-5), "#ff9a00"],
                     [normAGB(-0.5), "white"],
                     [normAGB(0.0), "white"],
                     [normAGB(0.5), "white"],
                     [normAGB(10), "tab:blue"]]
        cmapAGB = mpl.colors.LinearSegmentedColormap.from_list("", colorsAGB)

        norm = mpl.colors.Normalize(-180, 0.0)
        colorsMCWD = [[norm(-180.0), (0.65, 0.16, 0.16)],
                      [norm(-100.0), (0.80, 0.58, 0.047)],
                      [norm(-25.0), "white"],
                      [norm(0.0), "white"], ]
        cmapMCWD = mpl.colors.LinearSegmentedColormap.from_list("", colorsMCWD)

        yearbegin = 1985
        yearEnd = 2010

        datasetNames = setup.hydFileNames

        fig = plt.figure(figsize=(11, 8))
        index = 1

        indexes = [1, 3, 2, 4]

        import string
        lowerletters = string.ascii_lowercase[:26]

        for d in range(0, 2):
            dB, dA = listAnalysis[d].GetDataSlicesCond("cmasstotal", yearbegin, yearEnd, 3)
            fracAGB = self.GetAGBFractionCMassTotal(dB)
            diff = dA - dB
            diff*= fracAGB

            #nppsB, nppsA = listAnalysis[d].GetDataSlices("npp_sap", yearbegin, yearEnd, 0)
            #npphB, npphA= listAnalysis[d].GetDataSlices("npp_heart", yearbegin, yearEnd, 0)

            #npp_wood = nppsA + npphA
            #npp_wood[dB < 1.0] = np.nan
            #npp_wood*=fracAGB
            #mort = listAnalysis[d].GetDataSlicesYear("cmass_loss_heart_bg", yearEnd, 0)
            #mort += listAnalysis[d].GetDataSlicesYear("cmass_loss_heart_cav", yearEnd, 0)
            #mort += listAnalysis[d].GetDataSlicesYear("cmass_loss_heart_greff", yearEnd, 0)
            #mort += listAnalysis[d].GetDataSlicesYear("cmass_loss_sap_bg", yearEnd, 0)
            #mort += listAnalysis[d].GetDataSlicesYear("cmass_loss_sap_cav", yearEnd, 0)
            #mort += listAnalysis[d].GetDataSlicesYear("cmass_loss_sap_greff", yearEnd, 0)

            #mort[dB < 1.0] = np.nan
            #mort *= fracAGB
            #diff =  (npp_wood - mort)

            standardLPJs[d].GetBiomassBeforeAfter(yearbegin, yearEnd, True)
            dBStd = standardLPJs[d].BiomassBefore
            dAStd = standardLPJs[d].BiomassAfter

            fracAGBStd = self.GetAGBFractionCMassTotal(dBStd)
            diffStd = dAStd - dBStd
            diffStd *= fracAGBStd

            #Conversion from kgC/m² tp MgC/ha
            diffStd *= 10.0
            diff *= 10.0

            for i in range(1, 3):
                if i == 1:
                    diffMean = np.nanmean(diff, axis=0)
                    rightText = "LPJ-GUESS-HYD"
                elif i == 2:
                    diffMean = diffStd
                    rightText = "LPJ-GUESS"
                else:
                    # diffMean = listAnalysis[d].GetMCWD(yearEnd)
                    if yearEnd == 2010:
                        diffMean = data2010[d]
                        rightText = ""
                    elif yearEnd == 2005:
                        diffMean = data2005[d]
                        rightText = ""
                    else:
                        diffMean = np.zeros(1946)
                        rightText = ""

                img = listAnalysis[d].CreateImage(diffMean)
                img_extent = listAnalysis[d].stratifiedOutput.IMG_extent
                offset = [-3, 3, -3, 3]

                axGeo = fig.add_subplot(2, 2, indexes[index - 1], projection=ccrs.PlateCarree())
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
                axGeo.text(0.97, 0.92, datasetNames[d], horizontalalignment='right', verticalalignment='center',
                           transform=axGeo.transAxes, size=10,
                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
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

                axGeo.scatter(data_hubau_brienen['lon'], data_hubau_brienen['lat'], c = '#83f28f',s =14, ec = 'black')

                if (i == 1) | (i == 2):
                    imsh1 = axGeo.imshow(img, transform=ccrs.PlateCarree(), extent=img_extent, cmap=cmapAGB,
                                         norm=normAGB)
                index += 1

        plt.subplots_adjust(right=0.80, wspace=0.1, hspace=0.15)

        cax = plt.axes([0.82, 0.111, 0.02, 0.77])
        bar = plt.colorbar(imsh1, cax=cax, orientation="vertical")
        bar.set_label('$\Delta\mathrm{AGB}_{2010-1985}$ [MgC $\mathrm{ha}^{-1}$]', fontsize=12)

        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)

        plt.savefig(setup.folderPath + "/" + "Figure_Geospatial-AGB-"+str(yearbegin)+ "-" + str(yearEnd) + ".png", dpi=300, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        self.Succesfull()



