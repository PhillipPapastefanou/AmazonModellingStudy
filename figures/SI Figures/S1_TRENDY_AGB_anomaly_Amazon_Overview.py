
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


import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from A03_Hydraulics_Implementation.TRENDY.NetCDF_CVS import NetCDF_CVS_Reader
from A03_Hydraulics_Implementation.TRENDY.extraction import HalfDegreeExtractor
from A03_Hydraulics_Implementation.AutoAnalysis.PaperFigures_N2.Revisions2022.BiomassAmazonConverter \
    import BiomassAmazonAnomalyConverter


import os
rtpath = '../TRENDY/cVeg/conv'
files = [f for f in os.listdir(rtpath) if os.path.isfile(os.path.join(rtpath, f))]

file_str = rtpath + '/' + files[0]

ab_extractor = HalfDegreeExtractor('../Amazon_basin_05.tsv', '\t')
raisg_mask = r"../AB-SHAPE/Amazon_shape.shp"
mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                      ccrs.PlateCarree())

normAGB = mpl.colors.Normalize(-5, 5)
colorsAGB = [[normAGB(-5), "#9C0200"],
             [normAGB(-2.5), "#ff9a00"],
             [normAGB(-0.5), "white"],
             [normAGB(0.0), "white"],
             [normAGB(0.5), "white"],
             [normAGB(5), "tab:blue"]]
cmapAGB = mpl.colors.LinearSegmentedColormap.from_list("", colorsAGB)

yois = [2005, 2007, 2009, 2010]

for yoi in yois:
    index = 1
    fig = plt.figure(figsize=(16, 10))

    yearbeginAGB = 1985
    yearEndAGB = 2010


    for file in files:
        filename = os.path.join(rtpath, file)
        reader = NetCDF_CVS_Reader(filename, incomplete=False)

        cVeg = reader.GetCompleteSpaceOverTimeRange('cVeg', yearbeginAGB, yearEndAGB)
        unit_str = reader.GetVarUnit('cVeg')

        converter = BiomassAmazonAnomalyConverter()
        diff_anomaly = converter.Convert_CVS(cVeg, yoi, yearbeginAGB, yearEndAGB, yearbeginAGB)

        ab_extractor.GetCVSCoords(reader.lons, reader.lats)
        img_ab = ab_extractor.CreateCVSImage(diff_anomaly)
        img_extent = ab_extractor.IMG_extent

        if (unit_str == 'kg m-2') | (unit_str == 'kg m^-2') | (unit_str == 'kg/m$^2$'):
            multiplier = 1.0
        elif (unit_str == 'kgC/m2') | (unit_str == 'kgC m-2') | (unit_str == 'kg C m-2'):
            multiplier = 1.0
        else:
            print("Invalid cVeg unit")
            exit(99)

        # To convert from kg AGB C m-2 to Mg C AGB Ha-1
        img_ab *= 10.0
        # Flip according to the latitude axis
        img_ab = img_ab[::-1, :]

        axGeo = fig.add_subplot(4, 4, index, projection=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='g')
        lat_formatter = LatitudeFormatter()
        axGeo.xaxis.set_major_formatter(lon_formatter)
        axGeo.yaxis.set_major_formatter(lat_formatter)
        axGeo.add_feature(cfeature.BORDERS, edgecolor='tab:grey')
        axGeo.coastlines(resolution='110m', linewidth=1, color='tab:grey')
        axGeo.set_extent(np.array(img_extent), crs=ccrs.PlateCarree())
        axGeo.add_feature(mask, edgecolor='black', linewidth=1, facecolor="None")
        axGeo.text(0.97, 0.08, file[:len(file) - 11], horizontalalignment='right', verticalalignment='center',
                   transform=axGeo.transAxes, size=10,
                   bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
        imsh = axGeo.imshow(img_ab, transform=ccrs.PlateCarree(),
                            extent=img_extent, cmap=cmapAGB, norm=normAGB)

        index += 1

    cax = plt.axes([0.12, 0.15, 0.78, 0.02])

    bar = plt.colorbar(imsh, cax=cax, orientation="horizontal")
    bar.set_label('Absolute AGB anomaly [Mg/ha]', fontsize=12)
    plt.subplots_adjust(wspace=0.1, hspace=0.1, bottom=0.2)
    plt.savefig('cVeg_AB_Anomaly_' + str(yoi) + '.png', dpi=300, bbox_inches='tight',
                    pad_inches=0)
    plt.close()










