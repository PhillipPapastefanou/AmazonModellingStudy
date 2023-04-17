import numpy as np
from netCDF4 import Dataset
from enum import  Enum
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
from A03_Hydraulics_Implementation.TRENDY.NetCDF_CVS import NetCDF_CVS_Reader, CalendarType
from A03_Hydraulics_Implementation.TRENDY.NetCDF_CVS import TimeDomain
from A03_Hydraulics_Implementation.TRENDY.extraction import HalfDegreeExtractor
from A03_Hydraulics_Implementation.AutoAnalysis.Source.AmazonBasinStratified import AmazonBasinStratified
from Libs.Standard import StandardLPJG


from A03_Hydraulics_Implementation.AutoAnalysis.Source.BaseAnalysis import BaseAnalysis

class Single_Model_LAI_Deviation(BaseAnalysis):
    pass

    def rmse(self,y1,y2):
        return np.sqrt(np.nanmean((y1 - y2) ** 2))

    def __init__(self, setup):

        self.StartSetup()


        lowerletters = string.ascii_lowercase[:26]

        ab_extractor = HalfDegreeExtractor('/Users/pp/Dropbox/ClimateData/Coords/SA_Amazon_basin/Amazon_basin_05.tsv',
                                           '\t')

        raisg_mask = "/Users/pp/Dropbox/ClimateData/AmazonBasin/AB-SHAPE/Amazon_shape.shp"
        mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                              ccrs.PlateCarree())

        path = '/Users/pp/data/LAI/GeoTiff/1991_2020_last_three/Yearly_avg_1991_2020.nc'
        # path = '/Users/pp/data/LAI/GeoTiff/1991-2020_october2/Yearly_avg_1991_2020.nc'

        reader = NetCDF_CVS_Reader(path, incomplete=True)

        reader.yearBegin = 1991
        reader.yearEnd = 2020
        reader.multiplier = 1
        reader.calendar_type = CalendarType.standard
        reader.time_dim_name = 'time'



        standardLPJs = []
        for file in setup.stdFilePaths:
            standardLPJs.append(StandardLPJG(file, True))

        listAnalysis = []
        for i in range(0, len(setup.hydFilePaths)):
            listAnalysis.append(AmazonBasinStratified(setup.hydFilePaths[i], setup.hydFileParamsPaths[i]))




        index = 1
        fig = plt.figure(figsize=(10, 5))

        d = 0




        yois = [2010]

        index = 1

        for yoi in yois:
            agb_m_hyd = listAnalysis[d].GetDataSlicesCondContinous("lai", yoi, yoi, 0.0, -1)
            agb_m_std = standardLPJs[d].GetDataSlicesCondContinous("lai", yoi, yoi, 0.0)

            lai = reader.GetCompleteSpaceOverTimeRange('Band1', yoi, yoi)
            lai = lai.astype(float)
            lai[lai < 0.1] = np.nan
            lai /= 100.0

            ab_extractor.GetCVSCoords(reader.lons, reader.lats)

            agb_vod_ab = ab_extractor.Extract_OMA_data(lai, lon_first=False)

            w_greater = agb_m_hyd[:,:] > 0.1
            agb_m_hyd = agb_m_hyd[:,w_greater[0]]
            agb_m_std = agb_m_std[w_greater[0]]
            agb_vod_ab = agb_vod_ab[w_greater[0]]

            #agb_m_hyd_mean = np.average(agb_m_hyd, axis=0, weights=agb_m_hyd) * 0.8 * 10.0
            agb_m_hyd_mean = np.average(agb_m_hyd, axis=0, weights=agb_m_hyd)


            agb_m_std *= 1

            for j in range(0, 2):
                axGeo = fig.add_subplot(1, 2, index)

                if j == 0:
                    m_text = 'LPJ-GUESS-HYD'
                    data_y = agb_m_hyd_mean
                    #axGeo.scatter(agb_vod_ab, data_y)
                    data = np.zeros((2, 37*agb_m_hyd.shape[1]))
                    for i in range(0,37):
                        data[0, i*agb_m_hyd.shape[1]:(i+1)*agb_m_hyd.shape[1]] = agb_vod_ab
                        data[1, i*agb_m_hyd.shape[1]:(i+1)*agb_m_hyd.shape[1]]  = agb_m_hyd[i]
                    axGeo.scatter(data[0], data[1], c = 'tab:blue', alpha =0.1, s = 2)
                    #slope, intercept, r_value, p_value, std_err = stats.linregress(data[0], data[1])
                    axGeo.set_ylabel("LPJ-GUESS-HYD LAI []")
                    rmses = self.rmse(agb_m_hyd_mean, agb_vod_ab)


                else:
                    m_text = 'LPJ-GUESS-STD'
                    data_y = agb_m_std
                    axGeo.scatter(agb_vod_ab, data_y)
                    slope, intercept, r_value, p_value, std_err = stats.linregress(agb_vod_ab, data_y)
                    axGeo.set_ylabel("LPJ-GUESS-STD LAI []")
                    rmses = self.rmse(agb_m_std, agb_vod_ab)


                axGeo.plot(np.arange(4, 10,1), np.arange(4,10,1), c= 'black', ls = '--')
                axGeo.set_xlim((4, 9))
                axGeo.set_ylim((4, 9))


                axGeo.set_xlabel("Modis based LAI [MgC $\mathrm{ha}^{-1}$]")

                axGeo.text(0.04, 0.90, lowerletters[index - 1] + ")", horizontalalignment='left',
                           verticalalignment='center',
                           transform=axGeo.transAxes, size=18,
                           bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
                p_valuestr = '$RMSE = $ '
                p_valuestr += str(np.round(rmses, 1))

                axGeo.text(0.98, 0.07,
                   p_valuestr,
                    horizontalalignment='right', verticalalignment='center', c='black',
                   transform=axGeo.transAxes, size=14, bbox=dict(facecolor='white', alpha=1.0, edgecolor='white'))
                index += 1


        plt.subplots_adjust(wspace=0.2, hspace=0.1, bottom=0.2)
        plt.savefig(setup.folderPath + "/" + "Scatter_SI_LAI_MODEL_Deviation.png", dpi=300, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        self.Succesfull()
