
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
import scipy.stats as stats
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
rtpath = '/Users/pp/data/TRENDY/cVeg/conv'
files = [f for f in os.listdir(rtpath) if os.path.isfile(os.path.join(rtpath, f))]

file_str = rtpath + '/' + files[0]

ab_extractor = HalfDegreeExtractor('../Amazon_basin_05.tsv', '\t')

raisg_mask = r"../AB-SHAPE/Amazon_shape.shp"
mask = ShapelyFeature(Reader(raisg_mask).geometries(),
                      ccrs.PlateCarree())

normAGB = mpl.colors.Normalize(-5, 5)
colorsAGB = [[normAGB(-5), "#9C0200"],
             # [normAGB(-1.0), "#ff7400"],
             [normAGB(-2.5), "#ff9a00"],
             [normAGB(-0.5), "white"],
             [normAGB(0.0), "white"],
             [normAGB(0.5), "white"],
             [normAGB(5), "tab:blue"]]
cmapAGB = mpl.colors.LinearSegmentedColormap.from_list("", colorsAGB)


fig = plt.figure(figsize=(16, 10))



yearbeginAGB = 1984
yearEndAGB = 2010
index = 1

for file in files:

    ax = fig.add_subplot(4, 4, index)
    filename = os.path.join(rtpath, file)
    reader = NetCDF_CVS_Reader(filename)

    cVeg = reader.GetCompleteSpaceOverTimeRange('cVeg', yearbeginAGB, yearEndAGB)
    unit_str = reader.GetVarUnit('cVeg')

    converter = BiomassAmazonAnomalyConverter()
    cVeg_diff = converter.ConvertToAbsoluteDeltaAll_OMA(cVeg, axis=0)

    ab_extractor.GetCVSCoords(reader.lons, reader.lats)
    cVeg_ab_diff = ab_extractor.Extract_OMA_data(cVeg_diff)

    # Deal with all different units of the trendy outputs
    if (unit_str == 'kg m-2') | (unit_str == 'kg m^-2') | (unit_str == 'kg/m$^2$'):
        multiplier = 1.0
    elif (unit_str == 'kgC/m2') | (unit_str == 'kgC m-2') | (unit_str == 'kg C m-2'):
        multiplier = 1.0
    else:
        print("Invalid cVeg unit")
        exit(99)

    cVeg_ab_diff *= multiplier

    x = np.arange(1985, 2011)

    # To convert from kg AGB C m-2 to Mg C AGB Ha-1
    cVeg_ab_diff *= 10.0

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, np.nanmean(cVeg_ab_diff, axis=0))
    y = x * slope + intercept

    p_valuestr = str(np.round(slope, 3)) + " $"
    if p_value < 0.001:
        p_valuestr += "   (p-value < 0.001)"

    else:
        p_valuestr += "   (p-value = "
        p_valuestr += str(np.round(p_value, 3)) + ")"

    # titleTxt = axGeo.set_title(titles[i], size=16, loc = 'left')
    ax.text(0.97, 0.93, file[:len(file) - 11], horizontalalignment='right', verticalalignment='center',
            transform=ax.transAxes, size=10,
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

    
    ax.plot(x, np.nanmean(cVeg_ab_diff, axis=0), c='tab:blue', lw = 2)
    ax.plot(x, np.nanquantile(cVeg_ab_diff, axis=0, q = 0.25),c='tab:blue', lw = 1)
    ax.plot(x, np.nanquantile(cVeg_ab_diff, axis=0, q = 0.75),c='tab:blue', lw = 1)
    line = ax.plot(x, y, color='black', ls="--", alpha=0.5)

    ax.set_ylim((-0.75, 0.75))
    ax.text(0.03, 0.07, r"$\mathrm{slope} =  " + p_valuestr, horizontalalignment='left',
            verticalalignment='center', c="black",
            transform=ax.transAxes, size=8,)

    index += 1
plt.subplots_adjust(wspace=0.3, hspace=0.3, bottom=0.2)
plt.savefig('AB_NetSink_Trends.png', dpi=300, bbox_inches='tight',
            pad_inches=0)
plt.close()









