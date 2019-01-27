'''*********************************************
Authors: Alex Crawford
Date Created: 6/9/15
Date Modified: 5/11/18
Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell of a particular sector each year.

Inputs: 
    concentration threshold (concThresh) -- value between 0 and 1
    years of interest (years) -- two integers in a range function
    startmonth -- the month on which to start the annual cycle -- March (3) is
        a good idea because it is the closest to the maximum
    
Outputs: A csv file with the concentration threshold and setor noted in the 
    file name. Retreat and advance days are recorded as a "DOY", with 1 being
    the 1st day of January. If retreat below the
    concentration threshold never occurs, the minimum day is recorded instead.
*********************************************'''
# Import clock:
from time import clock
# Start script stopwatch. The clock starts running when time is imported
start = clock()

'''*******************************************
Set up Modules
*******************************************'''
print "Importing Modules"

import os
from osgeo import gdal, gdalconst, gdalnumeric
import numpy as np
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
print "Declaring Variables"

### Input Variables ###
cts = [0.15] # A number between 0 and 1 for the concentration threshold
n = 2 # Moving Average Size (n = # observation on either side of current day;
        ## so 1 = 3-point, 2 = 5-point, 4 = 9-point)

### Time Variables ###
ymin, ymax = 1986, 2017 # to 1998 to 2002
maxmo = [1,4]
minmo = [8,10]

### Path Variables ###
path = "/Volumes/MIRANDA/"
inpath = path+"/SeaIce/SSMI/Daily"
inpathnrt = path+"/SeaIce/SSMI/NearRealTime"
mapath = path+"/SeaIce/SSMI/SmoothedMA"+str(n*2+1)
outpath = path+"/SeaIce/AdvanceRetreat"

#### File Variables ###
LatN = "SSMI_Lats_Full.tif"
LonN = "SSMI_Lons_Full.tif"
maskN = "SSMI_Regions2.tif"
outputtype = gdal.GDT_Float64

ext = ".tif"
maskMin = 2
maskMax = 19

'''*******************************************
Main Analysis
*******************************************'''
print "Main Analysis"
# Load Lat and Long Rasters
lats = gdalnumeric.LoadFile(path+"/Projections/"+LatN)
lons = gdalnumeric.LoadFile(path+"/Projections/"+LonN)
mask = gdalnumeric.LoadFile(path+"/Projections/"+maskN)

ref = gdal.Open(path+"/Projections/"+LatN,gdalconst.GA_ReadOnly)
refG = ref.GetGeoTransform()

# Time Set Up
years = range(ymin,ymax+1)
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13",\
    "14","15","16","17","18","19","20","21","22","23","24","25","26","27",\
    "28","29","30","31"]
    
noretreat = np.zeros(len(years))
files = os.listdir(inpath) + os.listdir(inpathnrt)
files.sort()

for y in years:
    print " Year: "+str(y)
    dateref = [y,1,0]
    
    ### Prep Inputs ###
    # Load Current Year
    try: # If a smoothed SIC array already exists, open it
        arr1 = md.unpickle(mapath+"/SIC_SmoothedMA"+str(n*2+1)+"_"+str(y)+".pkl")
    
    except: # If not, you must create it anew first
        # Load Annual Files
        fileyList = [f for f in files if ( f.endswith(ext) and (int(f[3:7]) == y) )]
        fileyList_2 = [f for f in files if ( f.endswith(ext) and (int(f[3:7]) == y+1 and int(f[7:9]) == 1 and int(f[9:11]) <= n) )]
        fileyList_0 = [f for f in files if ( f.endswith(ext) and (int(f[3:7]) == y-1 and int(f[7:9]) == 12 and int(f[9:11]) > 31-n) )]
        # Load Arrays into a single 3-D array
        arrList = []
        for f in fileyList_0 + fileyList + fileyList_2: # Skip the last file because it's only referenced if no advance occurs
            # Reclassify array and mask out anything that's not in the regions of the mask
            try:
                arrList.append(md.reclassIceArray( gdalnumeric.LoadFile(inpath+"/"+f) )*np.where(np.where((mask <= maskMax),mask,np.nan) >= maskMin, 1, np.nan))
            except:
                arrList.append(md.reclassIceArray( gdalnumeric.LoadFile(inpathnrt+"/"+f) )*np.where(np.where((mask <= maskMax),mask,np.nan) >= maskMin, 1, np.nan))
            
        # Smooth array
        print "--- Applying moving average... (this may take an hour or so)"
        arr1 = np.apply_along_axis(md.movingAverage,0,np.array(arrList),n)
        if len(fileyList_2) == 0:
            ei = len(fileyList_0)+len(fileyList)
        else:
            ei = -len(fileyList_2)
        md.pickle(arr1[len(fileyList_0):ei,:,:],path+"/SeaIce/SSMI/SmoothedMA"+str(n*2+1)+"/SIC_SmoothedMA"+str(n*2+1)+"_"+str(y)+".pkl")
    
    # Load Following Year
    try:
        arr2 = md.unpickle(mapath+"/SIC_SmoothedMA"+str(n*2+1)+"_"+str(y+1)+".pkl")
 
    except: # If not, you must create it anew first
        # Load Annual Files
        fileyList = [f for f in files if ( f.endswith(ext) and (int(f[3:7]) == y+1) )]
        fileyList_2 = [f for f in files if ( f.endswith(ext) and (int(f[3:7]) == y+2 and int(f[7:9]) == 1 and int(f[9:11]) <= n) )]
        fileyList_0 = [f for f in files if ( f.endswith(ext) and (int(f[3:7]) == y and int(f[7:9]) == 12 and int(f[9:11]) > 31-n) )]
        # Load Arrays into a single 3-D array
        arrList = []
        for f in fileyList_0 + fileyList + fileyList_2: # Skip the last file because it's only referenced if no advance occurs
            # Reclassify array and mask out anything that's not in the regions of the mask
            try:
                arrList.append(md.reclassIceArray( gdalnumeric.LoadFile(inpath+"/"+f) )*np.where(np.where((mask <= maskMax),mask,np.nan) >= maskMin, 1, np.nan))
            except:
                arrList.append(md.reclassIceArray( gdalnumeric.LoadFile(inpathnrt+"/"+f) )*np.where(np.where((mask <= maskMax),mask,np.nan) >= maskMin, 1, np.nan))
            
        # Smooth array
        print "--- Applying moving average... (this may take an hour or so)"
        arr2 = np.apply_along_axis(md.movingAverage,0,np.array(arrList),n)
        if len(fileyList_2) == 0:
            ei = len(fileyList_0)+len(fileyList)
        else:
            ei = -len(fileyList_2)
        md.pickle(arr2[len(fileyList_0):ei,:,:],path+"/SeaIce/SSMI/SmoothedMA"+str(n*2+1)+"/SIC_SmoothedMA"+str(n*2+1)+"_"+str(y+1)+".pkl")
    
    # Load file list 
    fileList = [f for f in files if ( f.endswith(ext) and ( (int(f[3:7]) == y) or (int(f[3:7]) == y+1 ) ))]
    doys = np.array( [md.daysBetweenDates([y,1,1],[int(f[3:7]),int(f[7:9]),int(f[9:11])]) for f in fileList] )
    
    arr = np.concatenate( (arr1, arr2), axis=0 )
    
    # Calculate minimum & maximum value for year
    maxi0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y,maxmo[0]+1,1]) )[0][0]
    maxi1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y,maxmo[1]+1,1]) )[0][-1]
    
    mini0 = np.where( doys >= md.daysBetweenDates([y,1,1],[y,minmo[0],1]) )[0][0]
    mini1 = np.where( doys <= md.daysBetweenDates([y,1,1],[y,minmo[1]+1,1]) )[0][-1]
    
    maxi2 = np.where( doys >= md.daysBetweenDates([y,1,1],[y+1,maxmo[0],1]) )[0][0]
    maxi3 = np.where( doys <= md.daysBetweenDates([y,1,1],[y+1,maxmo[1]+1,1]) )[0][-1]
    
    Maxes1 = np.amax(arr[maxi0:maxi1],0)
    Mins = np.amin(arr[mini0:mini1],0)
    Maxes2 = np.amax(arr[maxi2:maxi3],0)
    
    print "--- Calculating Retreat/Advance..."
    for concThresh in cts:
        CT = str(int(concThresh*100))
        ### Prep Outputs ###
        mndArr, mxdArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
        frdArr, lrdArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
        fadArr, ladArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
        opArr, opcArr = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
        
        ### Calculate Retreat and Advance Events ###
        validcells = np.where( (np.isnan(Mins) == 0) & (np.isnan(Maxes1) == 0) & (np.isnan(Maxes2) == 0) )
        for i in range(len(validcells[0])):
            r,c = validcells[0][i], validcells[1][i] # Assign row and column
            
            # Calculate index for minimum & maximum
            MaxesI1 = int(np.median(np.where(arr[maxi0:maxi1,r,c] == Maxes1[r,c]))) # Gives first occurrence of minimum if multiples present
            MinsI = mini0 + int(np.median(np.where(arr[mini0:mini1,r,c] == Mins[r,c]))) # Gives first occurrence of minimum if multiples present
            MaxesI2 = maxi2 + int(np.median(np.where(arr[maxi2:maxi3,r,c] == Maxes2[r,c]))) # Gives first occurrence of minimum if multiples present
            
            # Store Minimum Day
            mndArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[MinsI][3:7]),int(fileList[MinsI][7:9]),int(fileList[MinsI][9:11]),0,0,0])
            mxdArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[MaxesI1][3:7]),int(fileList[MaxesI1][7:9]),int(fileList[MaxesI1][9:11]),0,0,0])
            mxd2rr = md.daysBetweenDates(dateref,[int(fileList[MaxesI2][3:7]),int(fileList[MaxesI2][7:9]),int(fileList[MaxesI2][9:11]),0,0,0])
            
            # If it's always above the concentration threshold...
            if Mins[r,c] >= concThresh: 
                opArr[r,c], opcArr[r,c] = 0, 0
                
                noretreat[y-min(years)] = noretreat[y-min(years)]+1
            
            # If it's never above the concentration threshold... 
            elif (Maxes1[r,c] < concThresh) & (Maxes2[r,c] < concThresh):
                opArr[r,c], opcArr[r,c] = 365, 365
                
            # Otherwise...
            else:
                above = np.where(arr[:,r,c] >= concThresh)[0] # Indices above concentration
                below = np.where(arr[:,r,c] < concThresh)[0] # Indices below concentration
                
                # First Retreat Day
                # First index after Maxes1 and before/on Mins for which concentration is below threshold
                try:
                    fri = below[np.where((below <= MinsI) & (below > MaxesI1))][0]
                    frdArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[fri][3:7]),int(fileList[fri][7:9]),int(fileList[fri][9:11]),0,0,0])
                except:
                    frdArr[r,c] = np.nan
                
                # Last Retreat Day
                # Last index after Maxes1 and before/on Mins for which concentration is below threshold
                try:
                    lri = above[np.where((above < MinsI) & (above >= MaxesI1))][-1]
                    lrdArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[lri][3:7]),int(fileList[lri][7:9]),int(fileList[lri][9:11])+1,0,0,0])
                except:
                    lrdArr[r,c] = np.nan
    
                # First Advance Day
                # First index after Mins and before/on Maxes2 for which concentration is above threshold
                try:
                    fai = above[np.where((above > MinsI) & (above <= MaxesI2))][0]  # Last index after minimum for which concentration is below threshold
                    fadArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[fai][3:7]),int(fileList[fai][7:9]),int(fileList[fai][9:11]),0,0,0])
                except:
                    fadArr[r,c] = np.nan
                
                # Last Advance Day
                # Last index after Mins anbd before/on Maxes2 for which concentration is below threshold
                try:
                    lai = below[np.where((below >= MinsI) & (below < MaxesI2))][-1]  # First index after minimum for which concentration is above threshold
                    ladArr[r,c] = md.daysBetweenDates(dateref,[int(fileList[lai][3:7]),int(fileList[lai][7:9]),int(fileList[lai][9:11])+1,0,0,0])
                except:
                    ladArr[r,c] = np.nan
                
                # Open Water Periods
                if (Maxes1[r,c] < concThresh)  & (Maxes2[r,c] >= concThresh): # When it starts above threshold but ends below
                    opArr[r,c] =  ladArr[r,c] - mxdArr[r,c] #np.min([365, ladArr[r,c] - MaxesI1])
                    opcArr[r,c] = fadArr[r,c] - mxdArr[r,c] #np.min([365, fadArr[r,c] - MaxesI1])
                    
                    lrdArr[r,c] = np.nan
                    frdArr[r,c] = np.nan
                    
                elif (Maxes1[r,c] >= concThresh)  & (Maxes2[r,c] < concThresh): # When it starts above threshold but ends below
                    opArr[r,c] =  mxd2rr - frdArr[r,c] #np.min([365, MaxesI2 - frdArr[r,c]])
                    opcArr[r,c] = mxd2rr - lrdArr[r,c] #np.min([365, MaxesI2 - lrdArr[r,c]])
                    
                    fadArr[r,c] = np.nan
                    ladArr[r,c] = np.nan
                    
                else: # Simple Case
                    opArr[r,c] =  ladArr[r,c] - frdArr[r,c] #np.min([365, ladArr[r,c] - frdArr[r,c]])
                    opcArr[r,c] =  fadArr[r,c] - lrdArr[r,c] #np.min([365, fadArr[r,c] - lrdArr[r,c]])
        
        ### Write Outputs ###
        md.writeNumpy_gdalObj(Mins,outpath+"/C"+CT+"/MIN/Value/Minimum_"+str(y)+ext,ref,outputtype)
        md.writeNumpy_gdalObj(mndArr,outpath+"/C"+CT+"/MND/Value/MinimumDay_"+str(y)+ext,ref,outputtype)
        md.writeNumpy_gdalObj(frdArr,outpath+"/C"+CT+"/FRD/Value/FirstRetreatDay_"+str(y)+ext,ref,outputtype)
        md.writeNumpy_gdalObj(lrdArr,outpath+"/C"+CT+"/LRD/Value/LastRetreatDay_"+str(y)+ext,ref,outputtype)
        md.writeNumpy_gdalObj(fadArr,outpath+"/C"+CT+"/FAD/Value/FirstAdvanceDay_"+str(y)+ext,ref,outputtype)
        md.writeNumpy_gdalObj(ladArr,outpath+"/C"+CT+"/LAD/Value/LastAdvanceDay_"+str(y)+ext,ref,outputtype)
        md.writeNumpy_gdalObj(opcArr,outpath+"/C"+CT+"/OPC/Value/OpenPeriodContinous_"+str(y)+ext,ref,outputtype)
        md.writeNumpy_gdalObj(opArr,outpath+"/C"+CT+"/OP/Value/OpenPeriod_"+str(y)+ext,ref,outputtype)
        
        # Remove objects from memory
        del fai, lai, fri, lri, frdArr, lrdArr, fadArr, ladArr, opArr, opcArr
        del MinsI, validcells, above, below, MaxesI1, MaxesI2, mxdArr, mxd2rr, mndArr
    
    del Maxes1, Maxes2, Mins, doys
    del maxi0, maxi1, maxi2, maxi3, mini0, mini1, arr, arr1, arr2

    # Print elapsed time
    print 'No Retreat Cells:'+ str(noretreat[y-min(years)]) +'; Elapsed time:',round(clock()-start,2),'seconds'

print "Complete."