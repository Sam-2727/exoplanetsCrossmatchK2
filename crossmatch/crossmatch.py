import pandas as pd
import numpy as np
from astropy.table import Table
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
import math
class Crossmatch():
    def __init__(self,exoplanetsFilename,dr2Filename,binariesFilename):
        self.eFile = exoplanetsFilename
        self.dr2File = dr2Filename
        self.bFile = binariesFilename
    def xmatchDr2(self):
        exoplanets = Table.read(self.eFile,format="fits")
        dr2 = Table.read(self.dr2File,format="fits")
        xmatchGaia = pd.DataFrame()
        eSub  = exoplanets[exoplanets["Label"]=="C"]
        for i in range(0,len(eSub["EPIC"])):
            toAppend = dr2[dr2["epic_number"]==eSub[eSub["Label"]=="C"]["EPIC"][i]]
            if(len(toAppend)>1):
                distances = np.sqrt(((dr2["ra"]-dr2["ra_epic"])*np.cos(dr2["dec"]*180/np.pi))**2+(dr2["dec"]-dr2["dec_epic"])**2)
                xmatchGaia = xmatchGaia.append(toAppend[distances==min(distances)])
            elif(len(toAppend)==0):
                xmatchGaia = xmatchGaia.append([[np.nan]*len(xmatchGaia.columns)],columns=xmatchGaia.columns)
            else:
                xmatchGaia  = xmatchGaia.append(toAppend)
        self.dr2xmatch = xmatchGaia.reset_index(inplace=True)
        #returning comparison for inspection
        return [eSub,self.dr2xmatch]
    def convertEDR3(self):
        designations = []
        for i in range(0,len(self.dr2xmatch)):
            designations.append(self.drxmatch["designation"].iloc[i][9:])
        job1 = Gaia.launch_job("SELECT * FROM gaiaedr3.dr2_neighbourhood WHERE dr2_source_id IN (" + \
            str(list(designations)).replace("\'","").replace(", nan","").replace("nan,","")[1:-1] + ")")
        resultsTotal = job1.get_results().to_pandas()
        edr3IDs = []
        for i in range(0,len(self.dr2xmatch)):
            toAppend = resultsTotal[resultsTotal["dr2_source_id"]==designations[i]]
            if(len(toAppend)==1):
                edr3IDs.append(toAppend["dr3_source_id"].iloc[0])
            else:
                edr3IDs.append(toAppend[toAppend["angular_distance"]==min(toAppend["angular_distance"])]["dr3_source_id"].iloc[0])
        self.edr3 = edr3IDs
        return edr3IDs
    def xmatch(self,criteria=None):
        binaries = Table.read("/Users/sam/Documents/Jupyter Notebooks/EDR3/all_columns_catalog.fits",format="fits").to_pandas()
        if(not criteria):
            criteria=pd.Series([True]*len(self.eSub))
        import math
        crossmatchesBinaries = pd.DataFrame()
        crossmatchesExoplanets = pd.DataFrame()
        crossmatchTFOP = pd.DataFrame()
        num = 0
        num2 = 0
        listHosts = []
        for i in range(0,len(self.edr3)):
            check1 = binaries[binaries["source_id1"]==self.edr3[i]]
            check2 = binaries[binaries["source_id2"]==self.edr3[i]]
            if(len(check1)==1):
                #print(i)
                num+=1
                crossmatchesBinaries = crossmatchesBinaries.append(check1)
                #crossmatchTFOP = crossmatchTFOP.append(reorg.iloc[i])
                crossmatchesExoplanets = crossmatchesExoplanets.append(self.eSub.iloc[i])
                listHosts.append(0)
            elif(len(check2)==1):
                #print(i)
                num2+=1
                crossmatchesBinaries = crossmatchesBinaries.append(check2)
                #crossmatchTFOP = crossmatchTFOP.append(reorg.iloc[i])
                crossmatchesExoplanets = crossmatchesExoplanets.append(self.eSub.iloc[i])
                listHosts.append(1)
        self.binariesC = crossmatchesBinaries[criteria].reset_index(inplace=True)
        self.exoplanetsC = crossmatchesExoplanets[criteria].reset_index(inplace=True)
        return [crossmatchesBinaries,crossmatchesExoplanets]
    def getTIC(self):
        temp1 = Catalogs.query_criteria(catalog="TIC",GAIA=self.binariesC["source_id1"])
        temp2 = Catalogs.query_criteria(catalog="TIC",GAIA=self.binariesC["source_id1"])
        TICc1 = pd.DataFrame()
        TICc2 = pd.DataFrame()
        for i in range(0,len(self.eSub)):
            toAppend1 = self.binariesC[self.binariesC["source_id1"]==i]
            if(len(toAppend1)>1):
                if(len(toAppend1[~np.isnan(toAppend1["Kmag"])])==1):
                    TICc1 = TICc1.append(toAppend1[~np.isnan(toAppend1["Kmag"])])
                elif(len(toAppend1[~np.isnan(toAppend1["Kmag"])])>1):
                    print("weird")
                    TICc1 = TICc1.append(toAppend1[~np.isnan(toAppend1["Kmag"])].iloc[0])
                else:
                    TICc1 = TICc1.append(toAppend1.iloc[0])
            else:
                TICc1.append(pd.DataFrame([[np.nan]*len(TICc1)]))
                print("weird2")
            toAppend2 = self.binariesC[self.binariesC["source_id2"]==i]
            if(len(toAppend2)>1):
                if(len(toAppend2[~np.isnan(toAppend2["Kmag"])])==1):
                    TICc2 = TICc2.append(toAppend2[~np.isnan(toAppend2["Kmag"])])
                elif(len(toAppend2[~np.isnan(toAppend2["Kmag"])])>1):
                    print("weird")
                    TICc2 = TICc2.append(toAppend2[~np.isnan(toAppend2["Kmag"])].iloc[0])
                else:
                    TICc2 = TICc2.append(toAppend2.iloc[0])
            else:
                TICc2.append(pd.DataFrame([[np.nan]*len(TICc2)]))
                print("weird2")
        self.TICc1 = TICc1
        self.TICc2  = TICc2
        return [TICc1,TICc2]

    
            
            