from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from Libs.Standard import StandardLPJG
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import matplotlib.colors as colors
from A03_Hydraulics_Implementation.AutoAnalysis.Source.SetupFiles import SetupFile
from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis


class AvitabileBiomassComparisonAbsolute(BaseAnalysis):
    pass

    def __init__(self, setup):

        self.StartSetup()
        self.GetAGBFractionCMassTotal = np.vectorize(self.getAGBFractionCMassTotal)

        raisg_mask = r"F:\Dropbox\ClimateData\AmazonBasin\AB-SHAPE\Amazon_shape.shp"
        mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                              ccrs.PlateCarree())

        standardLPJs = []
        for file in setup.stdFilePaths:
            standardLPJs.append(StandardLPJG(file, True))

        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))

        biomassDataRaw = pd.read_csv(r'F:\Dropbox\UNI\Projekte\A03_Hydraulics_Implementation\Avitabile\AvitabilieStationAGB.tsv',
                                     sep='\t',
                                     header=None).values

        aCoordX = biomassDataRaw[:, 0]
        aCoordY = biomassDataRaw[:, 1]

        from Libs.Scaling import ScalerListToArray

        biomassesAVI = biomassDataRaw[:, 2]
        # -1 because of mathemaitca indexing
        indexes = biomassDataRaw[:, 3].astype(int) - 1

        yearbegin = 2004
        yearEnd = 2005


        biomassesHyd = []
        biomassesStd = []

        datasetNames = setup.hydFileNames

        for d in range(0, 2):
            dB, dA = listAnalysis[d].GetDataSlicesCond("cmasstotal", yearbegin, yearEnd, 0)
            agb = self.GetAGBFractionCMassTotal(dB)*dB
            dataMedian = np.nanmedian(agb, axis=0)
            dataMean = np.nanmean(agb, axis=0)


            biomassesHyd.append(dataMean)

        for d in range(0, 2):
            standardLPJs[d].GetBiomassBeforeAfter(yearbegin, yearEnd)
            dBStd = standardLPJs[d].BiomassBefore

            agbStd = self.GetAGBFractionCMassTotal(dBStd)*dBStd

            biomassesHyd.append(agbStd)


        normAGB = mpl.colors.Normalize(2, 20)
        colorsAGB = [[normAGB(2), "tab:red"],
                     [normAGB(7), "yellow"],
                     [normAGB(12.5), "tab:green"],
                     [normAGB(20), "tab:blue"]]
        cmapAGB = mpl.colors.LinearSegmentedColormap.from_list("", colorsAGB)
        cmapAGB.set_under(color = "tab:red")

        import string
        lowerletters = string.ascii_lowercase[:26]
        fig = plt.figure(figsize=(10, 8))

        for i in range(0, len(biomassesHyd)):

            if  i < 2:
                rightText = "LPJ-GUESS-HYD"
            elif i == 2:
                rightText = "LPJ-GUESS"

            if np.mod(i,2)== 0:
                d = 0
            else :
                d = 1

            offset = [-3, 3, -3, 3]

            axGeo = fig.add_subplot(2, 2, i + 1, projection=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='g')
            lat_formatter = LatitudeFormatter()
            axGeo.xaxis.set_major_formatter(lon_formatter)
            axGeo.yaxis.set_major_formatter(lat_formatter)
            axGeo.add_feature(cfeature.BORDERS, edgecolor='tab:grey')
            axGeo.coastlines(resolution='110m', linewidth=1, color='tab:grey')
            # axGeo.set_title("Precipitation")




            axGeo.add_feature(mask, edgecolor='black', linewidth=1.3, facecolor="None")

            axGeo.set_xticks([-80, -70, -60, -50], crs=ccrs.PlateCarree())
            axGeo.set_yticks([-20, -15, -10, -5, 0, 5], crs=ccrs.PlateCarree())

            lonScaler = ScalerListToArray(aCoordX, 0.5, False)
            latScaler = ScalerListToArray(aCoordY, 0.5, False)
            lonIndexes = lonScaler.Indexes
            latIndexes = latScaler.Indexes
            xlen = lonScaler.len
            ylen = latScaler.len

            img =listAnalysis[0].CreateImage(biomassesHyd[i])
            IMG_extent = listAnalysis[0].stratifiedOutput.IMG_extent
            offset = [-3, 3, -3, 3]

            image = np.empty((ylen + 1, xlen + 1,)) * np.nan
            for j in range(len(aCoordX)):
                image[latIndexes[j], lonIndexes[j]] = biomassesAVI[j]

            axGeo.set_extent(list(np.array(IMG_extent) + np.array(offset)), crs=ccrs.PlateCarree())

            axGeo.scatter(aCoordX, aCoordY, c=biomassesAVI, cmap=cmapAGB, edgecolor='k', norm=normAGB)

            imsh = axGeo.imshow(img , transform=ccrs.PlateCarree(), extent=IMG_extent, cmap=cmapAGB, norm=normAGB)

            axGeo.text(0.97, 0.94, datasetNames[d], horizontalalignment='right', verticalalignment='center',
                       transform=axGeo.transAxes, size=10,
                       bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            axGeo.text(0.03, 0.94, lowerletters[i] + ")", horizontalalignment='left',
                       verticalalignment='center',
                       transform=axGeo.transAxes, size=14,
                       bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
            axGeo.text(0.12, 0.94, rightText, horizontalalignment='left', verticalalignment='center',
                       transform=axGeo.transAxes, size=10,
                       bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

        plt.subplots_adjust(right=0.80, wspace=0.1, hspace=0.3)
        cax = plt.axes([0.82, 0.11, 0.02, 0.77])
        bar = plt.colorbar(imsh, cax=cax, orientation="vertical")
        bar.set_label('Aboveground biomass (AGB) [kgC/m²]', fontsize=12)
        import os
        if not os.path.exists(setup.folderPath):
            os.mkdir(setup.folderPath)

        plt.savefig(setup.folderPath + "\\" + "AvitabileAbovegroundBiomass" + str(yearEnd) + ".png", dpi=300,
                    bbox_inches='tight',
                    pad_inches=0)
        plt.close()
        self.Succesfull()
