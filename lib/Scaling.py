
import numpy as np

class ScalerListToArray:

    def __init__(self, listOfGeo, res, reverse):
        self.listOfGeo = listOfGeo
        self.min, self.max = min(listOfGeo), max(listOfGeo)
        self.len = int((self.max - self.min) / res)
        self.a = (self.len) / (self.max - self.min)
        self.b = -self.min * (self.len) / (self.max - self.min)

        self.Indexes = listOfGeo.copy()
        self.Indexes *= self.a
        self.Indexes += self.b
        self.Indexes = self.Indexes.astype(int)


        self.Range = np.arange(self.min, self.max + res, res)

        if (reverse):
            self.Indexes = self.len - self.Indexes


class ScalerArrayToList:
    def __init__(self, coordList, dataArray, res, lons, lats, lonFirst):

        if(lats[0] < lats[1]):
            latIndexes = self.Rescale(coordList[:, 1], lats.shape[0] - 1, np.min(lats), np.max(lats), reverse = False)
        else:
            latIndexes = self.Rescale(coordList[:, 1], lats.shape[0] - 1, np.min(lats), np.max(lats), reverse = True)

        lonIndexes = self.Rescale(coordList[:, 0], lons.shape[0] - 1, np.min(lons), np.max(lons), reverse = False)




        self.Data = np.zeros(coordList.shape[0])

        if lonFirst:
            for i in range(0, coordList.shape[0]):
                self.Data[i] = dataArray[lonIndexes[i].astype(int), latIndexes[i].astype(int)]
        else:
            for i in range(0, coordList.shape[0]):
                self.Data[i] = dataArray[latIndexes[i].astype(int), lonIndexes[i].astype(int)]




    def Rescale(self, list, size,  min, max, reverse):
        indexes = list.copy()
        if reverse:
            indexes *= -(size) / (max - min)
        else:
            indexes *= (size) / (max - min)

        indexes -= min * size / (max - min)
        return indexes




