# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 13:39:10 2019

@author: Z


Day 1 
26 Proserpina Data
"""

#Exporting the data sets

#importing libraries
import urllib.request as url
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import os
from pathlib import Path
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from skimage import exposure
import skimage
import skimage.morphology as morph
from skimage.feature import peak_local_max
#from skimage import data, img_as_float
from scipy import ndimage as ndi
import math
from matplotlib.colors import LogNorm
from skimage.feature import peak_local_max
#import the files in the directory as a list

#FLAT
d1 = Path('Proserpina/26Proserpina/20120119/Day 1/FLAT')
day1 = os.listdir(d1)
headerd1 = fits.open(d1/ day1[0])
#getting data and median for the Proserpina data
for i in range(0,int(len(day1))):
    #i is the file
    zd1 = fits.getdata(d1/ day1[i])
    if i == 0:
        yd1 = zd1
    else:
        yd1 = np.dstack([yd1,zd1])
        #stacking the files and transposing them so that the number of files becomes insignificant
xd1 = yd1.transpose(2,0,1)
#median is a better measure from mean since it doesnt get affected by big offsets
FLATDAY1 = np.median(xd1,axis=0)




#Dark
d1 = Path('Proserpina/26Proserpina/20120119/Day 1/DARK')
day1 = os.listdir(d1)
headerd1 = fits.open(d1/ day1[0])
#getting data and median for the Proserpina data
for i in range(0,int(len(day1))):
    #i is the file
    zd1 = fits.getdata(d1/ day1[i])
    if i == 0:
        yd1 = zd1
    else:
        yd1 = np.dstack([yd1,zd1])
        #stacking the files and transposing them so that the number of files becomes insignificant
xd1 = yd1.transpose(2,0,1)
#median is a better measure from mean since it doesnt get affected by big offsets
DARKDAY1 = np.median(xd1,axis=0)


#Bias
d1 = Path('Proserpina/26Proserpina/20120119/Day 1/BIAS')
day1 = os.listdir(d1)
headerd1 = fits.open(d1/ day1[0])
#getting data and median for the Proserpina data
for i in range(0,int(len(day1))):
    #i is the file
    zd1 = fits.getdata(d1/ day1[i])
    if i == 0:
        yd1 = zd1
    else:
        yd1 = np.dstack([yd1,zd1])
        #stacking the files and transposing them so that the number of files becomes insignificant
xd1 = yd1.transpose(2,0,1)
#median is a better measure from mean since it doesnt get affected by big offsets
BIASDAY1 = np.median(xd1,axis=0)




#DATA
d1 = Path('Proserpina/26Proserpina/20120119/Day 1/DATA')
day1 = os.listdir(d1)
headerd1 = fits.open(d1/ day1[0])
#getting data and median for the Proserpina data
for i in range(0,int(len(day1))):
    #i is the file
    zd1 = fits.getdata(d1/ day1[i])
    if i == 0:
        yd1 = zd1
    else:
        yd1 = np.dstack([yd1,zd1])
        #stacking the files and transposing them so that the number of files becomes insignificant
xd1 = yd1.transpose(2,0,1)
#median is a better measure from mean since it doesnt get affected by big offsets
DATADAY1 = np.median(xd1,axis=0)


#######evaluationg#######

finaldata=((DATADAY1-BIASDAY1)-(DARKDAY1-BIASDAY1))/(FLATDAY1-BIASDAY1)
data=finaldata
#print("finaldata=", finaldata)

#####plot final data######
plt.figure(15,figsize=(15,15))
plt.title('Gray Calibrated Intensity Plot', fontsize=20)
plt.imshow(finaldata, cmap='gray', norm=LogNorm(), vmin=2)
plt.xlabel('pixel', fontsize=16)
plt.ylabel('pixel', fontsize=16)
plt.show()

plt.figure(15,figsize=(15,15))
plt.title('Reverse Gray Calibrated Intensity Plot', fontsize=20)
plt.imshow(finaldata, cmap='gray_r', norm=LogNorm(), vmin=0.005)
plt.xlabel('pixel', fontsize=16)
plt.ylabel('pixel', fontsize=16)
plt.show()

#################program####################

#centroid
def stdev():
    
   sigma=0
   Mean_value = (np.sum(data))/((len(data))**2)
   for x in range(len(data)):
        for y in range(len(data)):
            sigma = sigma+((data[x][y])-(Mean_value))**2
   sigma=np.sqrt(sigma/((len(data))**2))
   print ('Mean_value for NGC7331:', Mean_value)
   print ('sigma for NGC731 :', sigma)
   return data , Mean_value, sigma
#print("STD for each single data are as follow:",  stdev())

 



#using this function to evaluate the Peak of each fence area
def Peak_Fence(fence):
    
    X_Peak =0
    Y_Peak = 0 
    Total_Peak = 0.0
    for x in range(len(fence)-1):
        for y in range(len(fence)-1):
            if (fence[y][x] > Total_Peak):
                Total_Peak = fence[y][x]
                X_Peak= x
                Y_Peak= y
                #print('X Peaks for fence :' ,X_Peak)
                #print('Y Peaks for fence :' ,Y_Peak)
    return X_Peak , Y_Peak
            



#Placing the dear centroeid in our fence
def place_centero(fence, Border):
    fenceborder=len(fence)
    sumxoffset= 0.0
    sumxnoofset=0.0
    sumyoffset=0.0
    sumynoofset=0.0
    for x in range(fenceborder):
        for y in range(fenceborder):
            if (fence[y][x]> Border):            
                sumxoffset= sumxoffset +((x+1)*(fence[y][x]))
                sumxnoofset= sumxnoofset + (fence[y][x])
    for y in range(fenceborder):
        for x in range(fenceborder):
            if (fence[y][x]> Border):
                sumyoffset= sumyoffset +((y+1)*(fence[y][x]))
                sumynoofset= sumynoofset + (fence[y][x])
    finalcenterx=(sumxoffset/sumxnoofset)-1    
    finalcentery=(sumyoffset/sumynoofset)-1 
    #print("center x is =", finalcenterx, "final center y is=", finalcentery)
    return finalcenterx, finalcentery
  
    




def Fence(fence, Mean_value, sigma):
   Border = Mean_value + 3*sigma # cover up to %95 data
   #print("Border=",Border)
   X_Peak , Y_Peak = Peak_Fence(fence) 
   Total_list= []
   for x in range(len(fence)-1):
       list= 0 
       for y in range(len(fence)-1):
           if (fence[y][x] > Border):
               if (x==0) or (x==len(fence)):
                   #print("Broken fence because of x")
                   return []
               if (y==0) or (y==len(fence)):
                   #print("Broken fence becauswe of y")
                   return []
               else:
                   list = list + 1 #updating the list
               
           else:
               Total_list.append(list)
               
               #print("total list=", Total_list)
   width = max(Total_list)
   #print("width=", width)
   if (width<4):
      #print("Broken fence / Not enough pixel fitted for width problem")   
      return []     
    #not cosidered -1 and +2    
   Render_Fence = fence[int(Y_Peak-(width/2)-1):int(Y_Peak+(width/2)+2),int(X_Peak-(width/2)-1):int(X_Peak+(width/2)+2)]
   #print("Render_fence is=", Render_Fence)
   if (len(Render_Fence)==0):
       #print("Broken fence / render fence is zero size") 
       return []
   if (len(Render_Fence[0])==0):
       #print("Broken fence / render fence is zero shape") 
       return []
   Rad_Render_Fence = (len(fence)-width)/2
   if (X_Peak + Rad_Render_Fence) > len(fence):
       #print("Broken fence / Over radius")
       return []
   if (Y_Peak + Rad_Render_Fence)  > len(fence):
       #print("Broken fence / Over radius")
       return []
   for x in range(len(Render_Fence[0])):
        for y in range(len(Render_Fence)):
          if (Render_Fence[y][x] > Border):
            if (x==0) or (x==len(Render_Fence)):
                #print("Broken fence / Over radius")
                return []
            if (y==0) or (y==len(Render_Fence)):
                #print("Broken fence / Over radius")
                return [] 
   return Render_Fence        
            


def finalcenteroid():
    data, Mean_value, sigma = stdev()
    fenceborder= 40
    compset= int(fenceborder/3)+1
    Border = Mean_value + 3*sigma
    #print ( "Border [ADU]:" , Border)            
    BX = ([])
    BY = ([])
    finalcenterox= ([])
    finalcenteroy= ([])
    c= 0
    for x in range (24 , 875): #(24,2024)
        for y in range(24 , 2024): #(24,2024)
            if (data[y][x]>Border):
                if (data[y-1][x]<Border):
                    if (data[y][x-1]<Border):
                        #print("satisfactional result as follow: x , y =", x , y )
                        fence= data[y-compset:y-compset+40,x-compset:x-compset+40]
                        #print(fence)
                        TMX, TMY = Peak_Fence(fence)
                        #print(TMX, TMY)
                        newbornfence= Fence(fence, Mean_value, sigma)
                        #print("newborn fence=", newbornfence)
                        if (len(newbornfence)>0):
                            if len(newbornfence[0]>0):
                                X_Peak , Y_Peak = Peak_Fence(newbornfence)
                                x_centro , y_centro = place_centero(newbornfence, Border)
                                isit = "true"
                                if ((x-compset+TMX) in BX) and ((y-compset+TMX) in BY):
                                    c= c+1
                                    #print("gave up this star")
                                    isit= "false"
                                if(x_centro % 2== 0.0):
                                    if((y_centro % 2 ==0.0)):
                                        #print("gave up this star")
                                        isit= "false"
                                if ((x-compset+TMX-X_Peak+x_centro)==x-compset+TMX):
                                    if ((y-compset+TMY-Y_Peak+y_centro)==y-compset+TMY):
                                        #print("gave up this star")
                                        isit= "false"
                                if(isit=="true"):
                                    BX.append(x-compset+TMX)
                                    BY.append(y-compset+TMY)
                                    finalcenterox.append(x-compset+TMX-X_Peak+x_centro)
                                    finalcenterox= list(dict.fromkeys(finalcenterox))
                                    finalcenteroy.append(y-compset+TMY-Y_Peak+y_centro)
                                    finalcenteroy= list(dict.fromkeys(finalcenteroy))
                                    c=c+1
                        else:
                            c= c+0
    for x in range (900 , 2024): #(24,2024)
        for y in range(24 , 2024): #(24,2024)
            if (data[y][x]>Border):
                if (data[y-1][x]<Border):
                    if (data[y][x-1]<Border):
                        #print("satisfactional result as follow: x , y =", x , y )
                        fence= data[y-compset:y-compset+40,x-compset:x-compset+40]
                        #print(fence)
                        TMX, TMY = Peak_Fence(fence)
                        #print(TMX, TMY)
                        newbornfence= Fence(fence, Mean_value, sigma)
                        #print("newborn fence=", newbornfence)
                        if (len(newbornfence)>0):
                            if len(newbornfence[0]>0):
                                X_Peak , Y_Peak = Peak_Fence(newbornfence)
                                x_centro , y_centro = place_centero(newbornfence, Border)
                                isit = "true"
                                if ((x-compset+TMX) in BX) and ((y-compset+TMX) in BY):
                                    c= c+1
                                    #print("gave up this star")
                                    isit= "false"
                                if(x_centro % 2== 0.0):
                                    if((y_centro % 2 ==0.0)):
                                        #print("gave up this star")
                                        isit= "false"
                                if ((x-compset+TMX-X_Peak+x_centro)==x-compset+TMX):
                                    if ((y-compset+TMY-Y_Peak+y_centro)==y-compset+TMY):
                                        #print("gave up this star")
                                        isit= "false"
                                if(isit=="true"):
                                    BX.append(x-compset+TMX)
                                    BY.append(y-compset+TMY)
                                    finalcenterox.append(x-compset+TMX-X_Peak+x_centro)
                                    finalcenterox= list(dict.fromkeys(finalcenterox))
                                    finalcenteroy.append(y-compset+TMY-Y_Peak+y_centro)
                                    finalcenteroy= list(dict.fromkeys(finalcenteroy))
                                    c=c+1
                        else:
                            c= c+0
    #return finalcenterox, finalcenteroy           
    #print("x centeroid",finalcenterox)
    #print("y centeroid",finalcenteroy)  

    plt.figure(15,figsize=(15,15))
    plt.title('Centroid', fontsize=20)
    plt.scatter(finalcenterox,finalcenteroy, s= 50, facecolors = 'none', edgecolors = 'red')
    #plt.scatter(rad, ded, marker='+')
    plt.xlabel('pixel', fontsize=16)
    plt.ylabel('pixel', fontsize=16)
    #plt.show()
      
    plt.figure(15,figsize=(15,15))
    plt.title('Reverse Gray Calibrated Intensity Plot', fontsize=20)
    plt.imshow(data, cmap='gray_r', norm=LogNorm(),)
    plt.xlabel('pixel', fontsize=16)
    plt.ylabel('pixel', fontsize=16)
    plt.show()

    plt.figure(15,figsize=(15,15))
    plt.title('Reverse Gray Calibrated Intensity Plot', fontsize=20)
    plt.imshow(data, cmap='gray_r', norm=LogNorm(), )
    plt.xlabel('pixel', fontsize=16)
    plt.ylabel('pixel', fontsize=16)
    plt.show()

     

       
    return finalcenterox, finalcenteroy


    
#print("final centroid location x, y =",finalcenteroid())    
#print(finalcenteroid())

centx, centy = finalcenteroid()



#print("final center for x=", centx)
#print("final center for y=", centy)

#print("data", data)


#######################################
#usnbo part ---------------------------------------
#function: query Vizier USNO-B1 catalog

#Data files
fits1='NGC7331-S001-R001-C001-r.fts' 
x = fits.getdata(fits1)
hdr = fits.getheader(fits1)
file = fits.open('NGC7331-S001-R001-C001-r.fts' )
#print(file[0].header)

def usno(radeg,decdeg,fovam): # RA/Dec in decimal degrees/J2000.0 FOV in arc min. 

    import numpy as np
    import urllib.request as url
    
    #url format for USNO
    str1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=USNO-B1' 
    str2 = '&-c.ra={:4.6f}&-c.dec={:4.6f}&-c.bm={:4.7f}/{:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)
    #final URL: make sure it does not have any spaces or carriage returns/line feeds when copy-pasting
    sr = str1+str2 

    # Read from the webpage, storing the page's contents in 's'. 
    f = url.urlopen(sr)
    s = f.read() 
    f.close()

    #column interpretation compass
    namecol, RAcol, DECcol, rband = 0, 1, 2, 12
    null1, null2 = '     ', ''
    
    #split webpage string into lines
    sl = s.splitlines() 
    sl = sl[45:-1] # get rid of header 
    name = np.array([]) 
    rad = np.array([]) # RA in degrees 
    ded = np.array([]) # DEC in degrees
    rmag = np.array([]) # rmage 
    #get data from each line
    for k in sl:
        kw = k.decode().split('\t')
        if kw[0] != '':  
            name = np.append(name,kw[namecol])
            rad = np.append(rad,float(kw[RAcol])) 
            ded = np.append(ded,float(kw[DECcol]))
            # deal with case where no mag is reported
            if (kw[rband] != null1) and (kw[rband] != null2):
                rmag = np.append(rmag,float(kw[rband]))
            else:
                rmag = np.append(rmag,np.nan)
                
    #return data
    return name,rad,ded,rmag

#main function
if __name__ == '__main__':

    ras= file[0].header['ra']
    des = file[0].header['dec']
    radeg = 15*( float ( ras [ 0 : 2 ] ) + (float ( ras [3:5])/60.) +( float ( ras [6:])/3600.))
    dsgn = np.sign ( float ( des [ 0 : 3 ] ) )
    dedeg = float ( des [ 0 : 3 ] ) + dsgn*float ( des [4:6])/60. + dsgn*float ( des [7:])/3600.
            #print ’RA, DEC [ deg ] =’, radeg , dedeg
    #get stars in a 5'x5' square around (RA, DEC = 30, 30)
    name,rad,ded,rmag = usno(radeg,dedeg,34.55)
    #print(rmag) #there's missing mags, but you don't need mags anyways.
    #print("len rmag=",len(rmag))
    #plot positions
    #plt.figure(15,figsize=(12,12))
    mask=np.where(rmag<15)[0]
    rad=rad[mask]
    ded=ded[mask]
    plt.title("USNO-B1 Stars")
    plt.scatter(rad*100, ded*100, marker='+')
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    plt.tight_layout()
    plt.show()


'''
#convert RA & DEC coordinates into [ deg ]
def radec_deg( ras , des ):
            radeg = 15*( float ( ras [ 0 : 2 ] ) + float ( ras [3:5])/60. + float ( ras [6:])/3600.)
            dsgn = np . sign ( float ( des [ 0 : 3 ] ) )
            dedeg = float ( des [ 0 : 3 ] ) + dsgn*float ( des [4:6])/60. + dsgn*float ( des [7:])/3600.
            #print ’RA, DEC [ deg ] =’, radeg , dedeg
            return radeg , dedeg
        
'''        
        

#Converting coordinates
def coordinatecon(RA, DEC):
    X=0
    Y=0
    RA=np.radians(RA)
    DEC=np.radians(DEC)
    a0 = np.radians(radeg)
    d0 = np.radians(dedeg)
    X=-((np.cos(DEC))*(np.sin(RA-a0)))/(np.cos(d0)*np.cos(DEC)*np.cos(RA-a0)+np.sin(DEC)*np.sin(d0))
    Y=-((np.sin(d0)*np.cos(DEC)*np.cos(RA-a0)-(np.cos(d0)*np.sin(DEC))))/(np.cos(d0)*np.cos(DEC)*np.cos(RA-a0)+np.sin(DEC)*np.sin(d0))
    P=0.018
    XO=1024
    YO=1024
    F=3454
    X=F*(X/P)+XO
    Y=F*(Y/P)+YO
    
    
    return X,Y



X,Y = coordinatecon(rad,ded)
#print("rad=",rad)
#print("ded=",ded)
#print("Converted rad and ded",coordinatecon(rad,ded))

plt.figure(15,figsize=(12,12))
plt.title("USNO-B1 Stars")
plt.scatter(X,Y, marker='+')
plt.xlabel("X")
plt.ylabel("Y")
plt.tight_layout()
plt.show()
 #comtime
 
 
 
plt.figure(15,figsize=(12,12))
plt.title("SHOW TIME")
plt.scatter(X,Y, marker='+')
plt.scatter(centx,centy, s= 50, facecolors = 'none', edgecolors = 'red')
plt.xlabel("X")
plt.ylabel("Y")
plt.tight_layout()
plt.show()


#print(len(centx))
#print(len(centy))
#print(len(X))
#print(len(Y))
#print(len(rad))


xstar=[]
ystar=[]
xlist=[]
ylist=[] 
xoffset=[]
yoffset=[]
#offsetx= (centx-X)
#offsety= (centy-Y)

for i in range(len(centx)):
    centym=centy[i]
    centxm=centx[i]
    #centym=centy[i]
    for j in range(len(rad)):
        Ym=Y[j]
        Xm=X[j]
        #Ym=Y[j]
        sepration= np.sqrt(((centxm-Xm)**2)+((centym-Ym)**2))
        if sepration < 18:
            xstar.append(centxm)
            ystar.append(centym)
            xlist.append(Xm)
            ylist.append(Ym)
            xoffset.append(Xm-centxm)
            yoffset.append(Ym-centym)




plt.figure(15,figsize=(12,12))
plt.title("Overlaping the plots")
plt.xlabel("x pixel")
plt.ylabel("y pixel")  
plt.ylim(0, 2000)     
plt.xlim(0, 2000)     
plt.scatter(xlist, ylist,c='red', marker='^', label="centroid positions")
plt.scatter(xstar, ystar,c='blue',marker='v', label="star positions")
plt.legend(loc='best')
plt.show()
            


plt.figure(15,figsize=(12,12))
plt.title("yoffset the plots")
plt.xlabel("x pixel")
plt.ylabel("y pixel")        
plt.scatter(xlist, xoffset,c='red', marker='v', label="x")
plt.scatter(ylist, yoffset,c='blue',marker='^', label="y")
plt.legend(loc='best')
plt.show()
              

          
            
            


#Plate solution . Returns Plate Constants .
def platesolve ( xpix , ypix , xcent , ycent ):
    xpixels = xpix [ : ] # make independant copy of xpixels & ypixels
    ypixels = ypix [ : ] # ( pixel coordinates from USNO)
    ppr = 3454.00/0.018 #nominal value of pixels per radian , f /p
    for j in range ( len ( xpix )): #change from x & y [ pixel coord ] to X & Y [ standard coord ] 
        xpix [ j ] = ( xpix [ j ]-1024.0)/ppr
        ypix [ j ] = ( ypix [ j ]-1024.0)/ppr
    a , d = np.zeros(len(xpix)) , np.zeros(len(ypix))
    B, D = np.zeros((len(xpix),3)) , np.zeros((len(ypix),3))
    for k in range (len(xcent)):
        a[k] = xcent[k]
        for l in range(len(ycent)):
            d[l] = ycent[l]
            for m in range (len(xpix)):
                B[m][0]=ppr*xpix[m]
                B[m][1]=ppr*ypix[m]
                B[m][2]=1
            for n in range(len(ypix)):
                D[n][0]=ppr*xpix[n]
                D[n][1]=ppr*ypix[n]
                D[n][2]=1
            BtBinv = np.linalg.inv(np.dot(np.transpose (B) ,B))#for x pixels
            Bta = np.dot(np.transpose (B) ,a)
            cx = np.dot(BtBinv , Bta)
            DtDinv = np.linalg.inv(np.dot(np.transpose (D) ,D))#for y pixels
            Dtd = np.dot (np.transpose (D) ,d)
            cy = np.dot (DtDinv ,Dtd)
            #print("Plate Constants:")
            #print("a[1][1] , a[1][2] , x_0 =", cx)
            #print("a[2][1] , a[2][2] , y_0 =", cy)
            T = [[ppr*cx[0] ,ppr*cx[1] , cx[2]] , [ppr*cy[0] , ppr*cy[1] , cy[2]], [ 0 , 0 , 1 ]]
            #print ("\nT Matrix :\n " , T)
            ETS = np.sqrt (np.abs(np.linalg.det(T)))
            #print ("\nEffective Telescope Scale : " , ETS, "[ pix ]/[ rad ] vs " , ppr)
            chix2 = np.dot((np.transpose (a-(np.dot (B, cx )))) ,( a-(np.dot (B, cx ))))
            chiy2 = np.dot((np.transpose (d-(np.dot (D, cy )))) ,(d-(np.dot (D, cy ))))
            #print ("\nChiˆ2 x =", chix2 /( len ( xpix)-3)) # /dof for reduced chi−squared values ?
            #print ("\nChiˆ2 y =", chiy2 /( len ( ypix)-3))
            return T,ETS, xpixels , ypixels
        

        
        
    
T=platesolve(xlist,ylist,xstar,ystar)[0]
x=np.array([xstar[:],ystar[:],np.ones(len(ystar))[:]])
xlist= np.array([xlist])
#print(np.linalg.det(T))  
Tinverse=np.linalg.inv(T)
Tinverse=np.array(Tinverse)
#print("x",x)
#print("Tinverse",Tinverse)
#new=Tinverse.dot(x)
new=np.dot(Tinverse,x)
#print("new=",new)


