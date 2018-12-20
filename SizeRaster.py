# -*- coding: utf-8 -*-
"""
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
import SizerPy as sp

if __name__ == '__main__': 
    '''   
    demonstrate of resizing (GDAL) raster data set by calling SizerPy  
    '''
    args=sp.RasterParser()                          #get parameters
    inRasterObj = sp.SizeRaster(args.inputRaster)   #initialize object
    inRasterObj.SizeBands(args.outRaster, factor=args.factor, \
                 d1Interpolate = args.d1Interpolate) #resizing and saving rasters
    print 'done'