
import time
import numpy as np

class BaseAnalysis(object):

    #def __init__(self):
    #    self.GetAGBFractionCMassTotal = np.vectorize(self.__getAGBFractionCMassTotal)

    def StartSetup(self):
        self.start = time.time()
        print("Starting:  " + self.getName() + "...")
    def Succesfull(self):
        self.end =  time.time()
        print("Done! (" + str(int(np.round(self.end - self.start))) + " seconds.)");


    def getName(self):
        return self.__class__.__name__

    #cmass in kgC/m²
    def getAGBFractionCMassTotal(self, x):
        if x < 0.05:
            return 0.0
        else:
            return 1.11997 - 0.40231* x ** (-0.0570775)



