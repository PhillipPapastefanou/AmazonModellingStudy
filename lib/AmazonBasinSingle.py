from Libs.OutputSampling import OutputFile
from Libs.SmartOutput.outputfile import TimeDomain
from Libs.Standard import MCWDCalc
import pandas
import numpy as np

class AmazonBasinSingle:
    def __init__(self, filePath, configheader = None):
        self.outputfile = OutputFile(filePath)
        self.filePath = filePath
        self.npft = self.outputfile.nPfts


    def CreateImage(self, data):
        return  self.outputfile.CreateImage(0.5, data)

    def GetDataSlices(self, varName, yearBefore, yearAfter, pftId = -1):
        if pftId == -1:
            pftId = self.npft - 1

        #dataBefore = self.outputfile.GetAllPFTTotal(varName, yearBefore, yearBefore, 0, pftId)
        #dataAfter = self.outputfile.GetAllPFTTotal(varName, yearAfter, yearAfter, 0, pftId)

        #return np.squeeze(dataBefore),  np.squeeze(dataAfter)

    def GetDataSlicesRange(self, varName, yearBefore, yearAfter, pftId = -1):
        if pftId == -1:
            pftId = self.npft - 1

        if pftId == self.npft-1:
            data = self.outputfile.GetAllPFTTotal(varName, yearBefore, yearAfter)
        return np.squeeze(data)

    def GetDataSlicesPFTTotal(self, varName, yearBefore, yearAfter):

        dataBefore = self.outputfile.GetAllPFTTotal(varName, yearBefore, yearBefore)
        dataAfter = self.outputfile.GetAllPFTTotal(varName, yearAfter, yearAfter)

        return np.squeeze(dataBefore),  np.squeeze(dataAfter)


    def GetDataSlicesYearPFTTotal(self, varName, year):

        data =  self.outputfile.GetAllPFTTotal(varName, year, year)

        return np.squeeze(data)

    def GetDataSlicesPatchYear(self, varName, year, pftId = -1):
        if pftId == -1:
            pftId = self.npft - 1

        data =  self.stratifiedOutput.GetAllGridcellsPatch(varName, year, year, 0, self.nIns)


        return np.squeeze(data)

    def GetDataSlicesYearPft(self, varName, year, pftID_begin, pftID_end):

        data = self.stratifiedOutput.GetAllGridcellsPft(varName, year, year, 0, self.nIns, pftID_begin, pftID_end)

        return np.squeeze(data)

    def GetDataSlicesCond(self, varName, yearBefore, yearAfter, less, pftId = -1):

        if pftId == -1:
            pftId = self.npft - 1

        dataBefore = self.stratifiedOutput.GetAllGridcells(varName, yearBefore, yearBefore, 0, self.nIns, pftId)
        dataAfter = self.stratifiedOutput.GetAllGridcells(varName, yearAfter, yearAfter, 0, self.nIns, pftId)

        ind = dataBefore < less

        dataBefore[ind] = np.NAN
        dataAfter[ind] = np.NAN

        return np.squeeze(dataBefore),  np.squeeze(dataAfter)

    def GetDataSlicesPatch(self, varName, yearBefore, yearAfter):

        dataBefore = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearBefore, yearBefore, 0, self.nIns)
        dataAfter = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearAfter, yearAfter, 0, self.nIns)

        return np.squeeze(dataBefore),  np.squeeze(dataAfter)

    def GetDataPatch(self, varName, yearBegin, yearEnd, insfile):

        dataBefore = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearBegin, yearEnd, insfile, insfile)

        return np.squeeze(dataBefore)

    def GetDataPatchNoIns(self, varName, yearBegin, yearEnd):

        dataBefore = self.outputfile.GetPatchSingleAllGridcells(varName, yearBegin, yearEnd)

        return np.squeeze(dataBefore)

    def GetMCWD(self, year):

        return MCWDCalc(self.filePath, year).DeltaMCWD