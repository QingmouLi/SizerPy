# -*- coding: utf-8 -*-
"""
SizerPy.py is created on Thu Feb 26 16:41:36 2017. It uses
the isPCA method to enlarge rasters to arbitrary size of the
original one.
            
We follow the GPLv2.

   
Disclaim: Though authors can use it in their research and 
applications, there is no any gurantee or implicit claiming 
for its rightness. It is on users' challenge.

refs :
    Li, Q. and Delher S. A. 2018. SizerPy: a Python tool for resizing raster 
        data by interpolation in the eigenspace, Computers & geosciences, xxx, 
        p. xx-xx
    li, Q. and Dehler S. A. 2015. Inverse Spatial principal component analysis 
        for geophysical data interpolation. Journal of Applied Geophysics, 115, 
        p.79-91.



authors:  Qingmou Li and Sonya A. Dehler
          Natrual Resourcese of Canada

                    Qingmou.Li@Canada.ca
                    sonya.dehler@canada.ca
 

"""
import scipy.interpolate
import numpy as np
import gdal
import sys
import argparse

def IsPCA(arr,factor, d1Interpolate):
    '''resize rows with a given factor using the isPCA method, 
       factor can be digital number, such as 
        0.1, 0.25, 0.5, 1, 5 10 
    d1Interpolate: 1D interpolation method
                  =0  linear
                  =1  slinear  #bilinear
                  =2  quadratic
                  =3 cubic spline
    '''
    U,S,V=np.linalg.svd(arr,full_matrices=False) 
    newRows = int(U.shape[0]*factor)     #-factor+1
    U1=np.zeros((newRows,U.shape[1]),dtype=np.float64)
    Ux  =np.linspace(0.,1.0,U.shape[0])                   
    U1x =np.linspace(0.,1.0,newRows)  
    for j in range(U.shape[1]):        
        Uy=U[:,j]        
        if  d1Interpolate  == 0:
            interp=scipy.interpolate.interp1d(Ux,Uy,kind='linear')            
        elif d1Interpolate == 1:
            interp=scipy.interpolate.interp1d(Ux,Uy,kind='slinear')    #first order polynomial        
        elif d1Interpolate == 2:
            interp=scipy.interpolate.interp1d(Ux,Uy,kind='quadratic')  #second order polynomial
        elif d1Interpolate == 3:
            interp=scipy.interpolate.interp1d(Ux,Uy,kind='cubic')      #third order spline polynomial
        else:
            interp=scipy.interpolate.interp1d(Ux,Uy,kind='linear') 
        U1[:,j]=interp(U1x)   
    return np.transpose(np.dot(U1[:,:],np.dot(np.diag(S[:]),V[:,:])))   
    
def IsPCASizer(arr,factor, d1Interpolate=1):
    """
    resize a matrix with a given factor in the eigenspace
       arr[r,c] is the input array,   
       factor is a real number as the times () to resize the 
       array, such as 0.125, 0.5, 2., 4., 5. 10. etc. 
       the return is a new array having shape of [r*factor,c*factor]. 
    """    
    return IsPCA(IsPCA(arr,factor,d1Interpolate),factor,d1Interpolate)   
        
class SizeRaster:
    '''
    A multiple band raster resizing class 
    '''
    def __init__(self,inRasName):
        '''
        class constructor with raster name as the input parameter
        once the class is initialized, it can be used to resize the 
        original raster any times
        inRasName:  raster name for GDAL supported raster formats
        '''
        self.inFileName= inRasName          #keep a copy of file name
                
    def SizeBands(self, outFileName,factor, d1Interpolate=1):
        '''
        outFileName:   resized raster in geotiff format
        factor:        real number for resizing the raster, such as 0.125,0.5,2,8,10 etc
        d1Interpolate: 1D interpolation method using in eigenspace
                       0: linear;  
                       1: bilinear;  
                       2: quatrtic;  
                       3 cubic
        '''
        src_ds = gdal.Open(self.inFileName, gdal.gdalconst.GA_ReadOnly) 
        if src_ds is None:
            print 'Unable to open ' + self.inFileName
            sys.exit(1)
        projectionfrom = src_ds.GetProjection()
        self.geotransform   = src_ds.GetGeoTransform()                  
        self.cols = src_ds.GetRasterBand(1).XSize
        self.rows = src_ds.GetRasterBand(1).YSize      
        bandCount = src_ds.RasterCount
        newRows   = int(self.rows*factor)
        newCols   = int(self.cols*factor)
        #create output dataset
        gtiff = gdal.GetDriverByName('GTiff')  
        outdataset= gtiff.Create(outFileName, newCols, newRows,bandCount, gdal.GDT_Float32) #datatype      
        outGeotransform = self.newGeotrans(factor)
        outdataset.SetGeoTransform(outGeotransform)
        outdataset.SetProjection(projectionfrom)
        
        for b in range(bandCount):              #loop through all bands
            inBand  = src_ds.GetRasterBand(b+1)
            outBand = outdataset.GetRasterBand(b+1)
            noDataValue = src_ds.GetRasterBand(b+1).GetNoDataValue()                  
            rasData = inBand.ReadAsArray()  #0,0,self.cols,self.rows)                
            rasData[(rasData == noDataValue)] = 0.0    
            sizedArr = IsPCASizer(rasData,factor, d1Interpolate)
            outBand.SetNoDataValue(noDataValue)
            outBand.WriteArray(sizedArr,0,0)
            inBand  = None
            outBand = None
            
        src_ds     = None     
        outdataset = None        

    def newGeotrans(self, Factor):  #colFactor,rowFactor):
        '''
        calculate the new geotransform after sizing the input 
        raster
        '''
        A=np.zeros((6,6),np.float64)
        b=np.zeros(6,np.float64)
        newRows=int(self.rows*Factor)
        newCols=int(self.cols*Factor)        
        #1 point     
        b[0] = self.geotransform[0]+self.geotransform[1]+self.geotransform[2]
        b[1] = self.geotransform[3]+self.geotransform[4]+self.geotransform[5]
        A[0] = np.array([1.0,0.0,0.0,0.0,0.0,0.0])
        A[1] = np.array([0.0,0.0,0.0,1.0,0.0,0.0])       
        #2 point
        b[2] = self.geotransform[0]+(self.cols)*self.geotransform[1]+\
            (self.rows)*self.geotransform[2]
        b[3] = self.geotransform[3]+(self.cols)*self.geotransform[4]+\
            (self.rows)*self.geotransform[5]        
        A[2] = np.array([1.0,newCols,newRows,0.0,0.0,0.0])
        A[3] = np.array([0.0,0.0,0.0, 1.0,newCols,newRows])
       #3 point xmax ymin        
        b[4] = self.geotransform[0]+(self.cols)*self.geotransform[1]+ \
            (self.rows)*self.geotransform[2]
        b[5] = self.geotransform[3]+(self.cols)*self.geotransform[4]+\
            self.geotransform[5]        
        A[4] = np.array([1.0,newCols,0.0,0.0,0.0,0.0])
        A[5] = np.array([0.0,0.0,0.0, 1.0, newCols,0.0])    
        x=np.linalg.solve(A,b)
        return x
        
class SizeBitmap:
    '''
    A multiple band bitmap resizing class 
    does not matter with this module
    '''
    pass    

class SeismicTraces:
    '''
    seismic traces revision  
    does not matter with this module
    '''
    pass    

def plotU_U1(arr,factor=5, d1Interpolate=0): #figure 2
    """
    plot spatial loading scores 
    does not matter with this module
    """
    pass
#
def plotV(arr): #
    """
    plot bases of eigenspace  
    does not matter with this module
    """
    pass

def RasArr(tifFile):
    """
    reading a gdal supported raster file and get an numpy array
    """
    src_ds = gdal.Open(tifFile, gdal.gdalconst.GA_ReadOnly) 
    if src_ds is None:
        print 'Unable to open ' + tifFile
        sys.exit(1)
    inBand  = src_ds.GetRasterBand(1)
    noDataValue = inBand.GetNoDataValue()                  
    rasData = inBand.ReadAsArray()
    rasData[(rasData == noDataValue)]=0.0
    return rasData
    
def plotU_U1_difInterpolate(arr,factor=5, d1Interpolate=0): #figure 2
    """    
    does not matter with this module
    """
    pass

    
def plotS(arr): #figure 
    """
    does not matter with this module
    """
    pass

   
def RasterParser():
    '''
    define and parse parameters for raster resizing
    '''    
    parser=argparse.ArgumentParser(add_help=True, version='1.0.0', prog=\
        'SizeRaster.py')
    parser.add_argument('-r','--rasterFileName',action='store',dest='inputRaster'\
       ,help='the input raster file name, such as d:/temp/t.tif')
    parser.add_argument('-o','--outRaster',action='store',dest='outRaster'\
       ,help='the ouput raster file name, such as d:/temp/t.tif')
    parser.add_argument('-f','--zoomFactor',action='store',dest='factor'\
       ,type=float, help='zoom factor of the input raster')  
    parser.add_argument('-i','--d1interpolate',action='store',dest='d1Interpolate' \
       ,type=int, default=1, help='1D interpolation method, \
                                   0: linear,\
                                   1: bilinear\
                                   2: quatradic \
                                   3: cubic spline') 
    args=parser.parse_args() 
    return args    

def bitmapParser():
    
    parser=argparse.ArgumentParser(add_help=True, version='1.0.0', prog=\
        'bitmapSizer.py')
    parser.add_argument('-b','--inputBitmapFileName',action='store',dest='inputBitmap'\
       ,help='the input bitmap file name, such as d:/temp/t.png')
    parser.add_argument('-o','--outBitmap',action='store',dest='outBitmap'\
       ,help='the ouput bitmap file name, such as d:/temp/t.png')
    parser.add_argument('-f','--zoomFactor',action='store',dest='factor'\
       ,type=float, help='zoom factor of the input raster')  
    parser.add_argument('-i','--d1interpolate',action='store',dest='d1Interpolate' \
       ,type=int, default=1, help='1D interpolation method, \
                                   0: linear,\
                                   1: bilinear\
                                   2: quatradic \
                                   3: cubic spline') 
    args=parser.parse_args() 
    return args      

if __name__ == "__main__":
    """
    this is model, ready for call by SizerPy.py
    """     
    print ("This is a model, the calling program name is SizerRaster.py\n")
    