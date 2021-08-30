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
        exoplanets = Table.read(self.eFile,format="fits").to_pandas()
        dr2 = Table.read(self.dr2File,format="fits").to_pandas()
        xmatchGaia = pd.DataFrame()
        eSub  = exoplanets[exoplanets["Label"]==b"C"]
        for i in range(0,len(eSub["EPIC"])):
            toAppend = dr2[dr2["epic_number"]==eSub["EPIC"].iloc[i]]
            if(len(toAppend)>1):
                distances = np.sqrt(((toAppend["ra"]-toAppend["ra_epic"])*np.cos(toAppend["dec"]*180/np.pi))**2+(toAppend["dec"]-toAppend["dec_epic"])**2)
                xmatchGaia = xmatchGaia.append(toAppend[distances==min(distances)].iloc[0])
            elif(len(toAppend)==0):
                xmatchGaia = xmatchGaia.append(pd.DataFrame([[np.nan]*len(xmatchGaia.columns)],columns=xmatchGaia.columns))
            else:
                xmatchGaia  = xmatchGaia.append(toAppend)
            lenXmatch = len(xmatchGaia)
        self.dr2xmatch = xmatchGaia.reset_index(drop=True)
        indexesNan = []
        for i in range(0,len(self.dr2xmatch)):
            try:
                if(np.isnan(self.dr2xmatch["designation"].iloc[i])):
                    indexesNan.append(i)
            except:
                pass
        self.dr2xmatch = self.dr2xmatch.drop(indexesNan).reset_index(drop=True)
        self.eSub = eSub.reset_index(drop=True).drop(indexesNan).reset_index(drop=True)
        #returning comparison for inspection
        return [eSub,self.dr2xmatch]
    def convertEDR3(self):
        designations = []
        for i in range(0,len(self.dr2xmatch)):
            designations.append(self.dr2xmatch["designation"].iloc[i][9:].decode('ascii'))
        job1 = Gaia.launch_job("SELECT * FROM gaiaedr3.dr2_neighbourhood WHERE dr2_source_id IN (" + \
            str(list(designations)).replace("\'","").replace(", nan","").replace("nan,","")[1:-1] + ")")
        resultsTotal = job1.get_results().to_pandas()
        edr3IDs = []
        for i in range(0,len(self.dr2xmatch)):
            toAppend = resultsTotal[resultsTotal["dr2_source_id"]==int(designations[i])]
            if(len(toAppend)==1):
                edr3IDs.append(toAppend["dr3_source_id"].iloc[0])
            else:
                edr3IDs.append(toAppend[toAppend["angular_distance"]==min(toAppend["angular_distance"])]["dr3_source_id"].iloc[0])
        self.edr3 = edr3IDs
        self.dr2xmatch["EDR3"] = self.edr3
        return edr3IDs
    def xmatch(self,criteria=None):
        binaries = Table.read("/Users/sam/Documents/Jupyter Notebooks/EDR3/all_columns_catalog.fits",format="fits").to_pandas()
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
        crossmatchesBinaries = crossmatchesBinaries.reset_index(drop=True)
        crossmatchesExoplanets = crossmatchesExoplanets.reset_index(drop=True)
        if(not criteria):
            criteria=(crossmatchesBinaries["ruwe1"]<1.4)&(crossmatchesBinaries["ruwe2"]<1.4)\
                & (crossmatchesBinaries["R_chance_align"]<0.1)
        self.binariesC = crossmatchesBinaries[criteria].reset_index(drop=True)
        self.exoplanetsC = crossmatchesExoplanets[criteria].reset_index(drop=True)
        return [crossmatchesBinaries,crossmatchesExoplanets]
    def getEDR3(self):
        job1 = Gaia.launch_job("SELECT * FROM gaiaedr3.dr2_neighbourhood WHERE dr3_source_id IN (" + \
            str(list(self.binariesC["source_id1"])).replace("\'","").replace(", nan","").replace("nan,","")[1:-1] + ")")
        resultsTotal = job1.get_results().to_pandas()
        dr2IDs1 = []
        edr3Reorg1 = pd.DataFrame()
        for i in range(0,len(self.binariesC)):
            toAppend = resultsTotal[resultsTotal["dr3_source_id"]==int(self.binariesC["source_id1"].iloc[i])]
            #print(i,toAppend)
            if(len(toAppend)==1):
                edr3Reorg1 = edr3Reorg1.append(toAppend.iloc[0])
                dr2IDs1.append(toAppend["dr2_source_id"].iloc[0])
            else:
                edr3Reorg1 = edr3Reorg1.append(toAppend[toAppend["angular_distance"]==min(toAppend["angular_distance"])].iloc[0])
                dr2IDs1.append(toAppend[toAppend["angular_distance"]==min(toAppend["angular_distance"])]["dr3_source_id"].iloc[0])

        job1 = Gaia.launch_job("SELECT * FROM gaiaedr3.dr2_neighbourhood WHERE dr3_source_id IN (" + \
        str(list(self.binariesC["source_id2"])).replace("\'","").replace(", nan","").replace("nan,","")[1:-1] + ")")
        resultsTotal = job1.get_results().to_pandas()
        dr2IDs2 = []
        edr3Reorg2 = pd.DataFrame()
        for i in range(0,len(self.binariesC)):
            toAppend = resultsTotal[resultsTotal["dr3_source_id"]==int(self.binariesC["source_id2"].iloc[i])]
            #print(i,toAppend)
            if(len(toAppend)==1):
                edr3Reorg2 = edr3Reorg2.append(toAppend.iloc[0])
                dr2IDs2.append(toAppend["dr2_source_id"].iloc[0])
            else:
                edr3Reorg2 = edr3Reorg2.append(toAppend[toAppend["angular_distance"]==min(toAppend["angular_distance"])].iloc[0])
                dr2IDs2.append(toAppend[toAppend["angular_distance"]==min(toAppend["angular_distance"])]["dr3_source_id"].iloc[0])
        self.binariesC["DR2_source_id1"] = dr2IDs1
        self.binariesC["DR2_source_id2"] = dr2IDs2
        edr3Reorg1 = edr3Reorg1.reset_index(drop=True)
        edr3Reorg2 = edr3Reorg2.reset_index(drop=True)
        criteria = (edr3Reorg1["magnitude_difference"]<1)&(edr3Reorg2["magnitude_difference"]<1)
        self.binariesC = self.binariesC[criteria].reset_index(drop=True)
        self.exoplanetsC= self.exoplanetsC[criteria].reset_index(drop=True)
        return edr3
    def getTIC(self):
        temp1 = Catalogs.query_criteria(catalog="TIC",GAIA=self.binariesC["DR2_source_id1"]).to_pandas()
        temp2 = Catalogs.query_criteria(catalog="TIC",GAIA=self.binariesC["DR2_source_id2"]).to_pandas()
        TICc1 = pd.DataFrame()
        TICc2 = pd.DataFrame()
        for i in range(0,len(self.binariesC)):
            toAppend1 = temp1[temp1["GAIA"]==str(self.binariesC["DR2_source_id1"].iloc[i])]
            #print(toAppend1)
            if(len(toAppend1)>1):
                if(len(toAppend1[~np.isnan(toAppend1["Kmag"])])==1):
                    TICc1 = TICc1.append(toAppend1[~np.isnan(toAppend1["Kmag"])])
                elif(len(toAppend1[~np.isnan(toAppend1["Kmag"])])>1):
                    print("weird")
                    TICc1 = TICc1.append(toAppend1[~np.isnan(toAppend1["Kmag"])].iloc[0])
                else:
                    TICc1 = TICc1.append(toAppend1.iloc[0])
            elif(len(toAppend1)==1):
                TICc1 = TICc1.append(toAppend1)
            else:
                TICc1 = TICc1.append(pd.DataFrame([[np.nan]*len(TICc1)]))
                print("weird2")
            toAppend2 = temp2[temp2["GAIA"]==str(self.binariesC["DR2_source_id2"].iloc[i])]
            if(len(toAppend2)>1):
                if(len(toAppend2[~np.isnan(toAppend2["Kmag"])])==1):
                    TICc2 = TICc2.append(toAppend2[~np.isnan(toAppend2["Kmag"])])
                elif(len(toAppend2[~np.isnan(toAppend2["Kmag"])])>1):
                    print("weird")
                    TICc2 = TICc2.append(toAppend2[~np.isnan(toAppend2["Kmag"])].iloc[0])
                else:
                    TICc2 = TICc2.append(toAppend2.iloc[0])
            elif(len(toAppend2)==1):
                TICc2 = TICc2.append(toAppend2)    
            else:
                TICc2.append(pd.DataFrame([[np.nan]*len(TICc2)]))
                print("weird2")
        self.TICc1 = TICc1
        self.TICc2  = TICc2
        return [TICc1,TICc2]
    
            
            