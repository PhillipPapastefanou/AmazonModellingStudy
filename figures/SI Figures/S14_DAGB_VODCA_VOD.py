import numpy as np
from netCDF4 import Dataset
from enum import  Enum

from A03_Hydraulics_Implementation.AutoAnalysis.PaperFigures_N2.Revisions2022.BiomassAmazonConverter import \
    BiomassAmazonAnomalyConverter
from Libs.Scaling import ScalerListToArray
from datetime import datetime
from datetime import timedelta
import Libs.SmartOutput.julian
import matplotlib.pyplot as plt
import matplotlib.colors
import os
from sys import platform
import string
from sklearn.metrics import mean_squared_error
import scipy.stats as stats




import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from A03_Hydraulics_Implementation.TRENDY.NetCDF_CVS import NetCDF_CVS_Reader
from A03_Hydraulics_Implementation.TRENDY.NetCDF_CVS import TimeDomain
from A03_Hydraulics_Implementation.TRENDY.extraction import HalfDegreeExtractor
from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from Libs.Standard import StandardLPJG


from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis

class Single_Model_DAGB_VODCA_Deviation(BaseAnalysis):
    pass

    def rmse(self,y1,y2):
        return np.sqrt(np.nanmean((y1 - y2) ** 2))

    def __init__(self, setup):

        self.StartSetup()

        biomassConverter = BiomassAmazonAnomalyConverter()

        lowerletters = string.ascii_lowercase[:26]

        ab_extractor = HalfDegreeExtractor('/Users/pp/Dropbox/ClimateData/Coords/SA_Amazon_basin/Amazon_basin_05.tsv',
                                           '\t')

        raisg_mask = "/Users/pp/Dropbox/ClimateData/AmazonBasin/AB-SHAPE/Amazon_shape.shp"
        mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                              ccrs.PlateCarree())

        filename = '/Users/pp/data/ClimateData/Vod/Vodca/Ku-band/Complete_Annually.nc'
        reader = NetCDF_CVS_Reader(filename, incomplete=False)



        standardLPJs = []
        for file in setup.stdFilePaths:
            standardLPJs.append(StandardLPJG(file, True))

        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))




        index = 1
        fig = plt.figure(figsize=(10, 5))

        d = 0

        yearbegin = 2000
        yearend = 2016
        vod_data = reader.GetCompleteSpaceOverTimeRange('vod', yearbegin, yearend)

        yois = [2010]
        yearbeginAGB = 1985
        yearEndAGB = 2010
        index = 1

        for yoi in yois:

            vod_data_slice = vod_data[yoi - yearbegin] - vod_data[yoi -yearbegin- 1]
            ab_extractor.GetCVSCoords(reader.lons, reader.lats)
            # ab_extractor.GetCVSCoords(reader.lats,reader.lons)
            diff_vod = ab_extractor.Extract_OMA_data(vod_data_slice, lon_first=False)

            biomass = listAnalysis[d].GetDataSlicesCondContinous("cmasstotal", yearbeginAGB, yearEndAGB, 1.0)
            # standardLPJs[d].GetBiomassBeforeAfter(yearbeginAGB, yearEndAGB, True)
            biomass_std = standardLPJs[d].GetDataSlicesCondContinous("cmasstotal", yearbeginAGB, yearEndAGB,
                                                                     1.0)* 0.8 * 10.0

            diff_anomaly = biomassConverter.ConvertToAbsoluteAnomaly_OMA(biomass, yoi, yearbeginAGB, yearEndAGB,
                                                                         yearbeginAGB)
            diff_anomaly_std = biomassConverter.ConvertToAbsoluteAnomaly_OMA(biomass_std, yoi, yearbeginAGB, yearEndAGB,
                                                                             yearbeginAGB)

            ab_extractor.GetCVSCoords(reader.lons, reader.lats)


            w_greater = biomass[:,:, 0] > 0.1
            biomass_weight = biomass[:,w_greater[0],0]
            diff_anomaly_hyd = diff_anomaly[:,w_greater[0]]
            diff_anomaly_std = diff_anomaly_std[w_greater[0]]
            diff_vod = diff_vod[w_greater[0]]




            #diff_anomaly_hyd_mean = np.average(diff_anomaly_hyd, axis=0, weights=biomass_weight) * 0.8 * 10.0
            diff_anomaly_hyd_mean = np.nanmean(diff_anomaly_hyd, axis=0) * 0.8 * 10.0


            #agb_m_std *= 0.8 * 10

            for j in range(0, 2):
                axGeo = fig.add_subplot(1, 2, index)

                if j == 0:
                    m_text = 'LPJ-GUESS-HYD'
                    data_y = diff_anomaly_hyd_mean
                    axGeo.scatter(diff_vod, data_y)
                    #data = np.zeros((2, 37*diff_anomaly_hyd.shape[1]))
                    #for i in range(0,37):
                    #    data[0, i*diff_anomaly_hyd.shape[1]:(i+1)*diff_anomaly_hyd.shape[1]] = diff_vod
                    #    data[1, i*diff_anomaly_hyd.shape[1]:(i+1)*diff_anomaly_hyd.shape[1]] =   diff_anomaly_hyd[i] *0.8 *10
                    #axGeo.scatter(data[0], data[1], c = 'tab:blue', alpha =0.1, s = 2)
                    ind_npt = (np.isnan(diff_vod)) |  (np.isnan(diff_anomaly_hyd_mean) )

                    slope, intercept, r_value, p_value, std_err = stats.linregress(diff_vod[~ind_npt], diff_anomaly_hyd_mean[~ind_npt])
                    axGeo.set_ylabel("AGB anomaly [MgC $\mathrm{ha}^{-1}$]")
                    #rmses = self.rmse(agb_m_hyd_mean, diff_vod)


                else:
                    m_text = 'LPJ-GUESS-STD'
                    data_y = diff_anomaly_std
                    axGeo.scatter(diff_vod, data_y)
                    ind_npt = (np.isnan(diff_vod)) |  (np.isnan(data_y) )
                    slope, intercept, r_value, p_value, std_err = stats.linregress(diff_vod[~ind_npt], data_y[~ind_npt])
                    axGeo.set_ylabel("AGB anomaly [MgC $\mathrm{ha}^{-1}$]")
                    #rmses = self.rmse(data_y, diff_vod)


                #axGeo.plot(np.arange(-0.15,0.15,0.01), np.arange(-0.15,0.15,0.01)*slope + intercept, c= 'black', ls = '--')
                axGeo.set_xlim((-0.11, 0.11))
                #axGeo.set_ylim((-12, 4))


                axGeo.set_xlabel("VOD anomaly [-]")

                axGeo.text(0.04, 0.92, lowerletters[index - 1] + ")", horizontalalignment='left',
                           verticalalignment='center',
                           transform=axGeo.transAxes, size=18,
                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                p_valuestr = 'p = '
                p_valuestr += str(np.round(p_value, 3))

                axGeo.text(0.98, 0.97,
                   m_text,
                     horizontalalignment='right', verticalalignment='top', c='black',
                transform=axGeo.transAxes, size=12, bbox=dict(facecolor='white', alpha=1.0, edgecolor='white'))
                index += 1


        plt.subplots_adjust(wspace=0.2, hspace=0.1, bottom=0.2)
        plt.savefig(setup.folderPath + "/" + "Scatter_SI_DAGB_VOD_Deviation.png", dpi=300, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        self.Succesfull()
