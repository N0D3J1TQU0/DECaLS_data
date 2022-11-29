import astropy.stats as st
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, join
import astropy.units as u
import numpy as np
from scipy.stats import chisquare
from astropy.coordinates import SkyCoord
from lmfit import Model
import warnings
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=68.3,Om0=0.299) #Bocquet et al. 2015 cosmology

def linear(x, A, B):
    return (A*x) + B
lmodel = Model(linear)

def model_dist(m,b,xlst,ylst):
    lst = []
    for i in range(len(xlst)):
        x0 = xlst[i]
        y0 = ylst[i]
        Mneg = -1/m
        Bneg = y0 - (Mneg*x0)
        x1 = (b-Bneg)/(Mneg-m)
        y1 = (m*x1)+b
        if y0>y1:
            lst += [np.sqrt(((x0-x1)**2) + ((y0-y1)**2))]
        else:
            lst += [(-1)*np.sqrt(((x0-x1)**2) + ((y0-y1)**2))]
    return lst

def col_evo_weights(z):
    #[gr,ri,iz,rz]
    if z<0.3:
        w = [1,0.33,0,0.33]
    if z>=0.3 and z<0.5:
        w = [0,1,0.33,0.66]
    if z>=0.5 and z<0.65:
        w = [0,0.66,0.33,1]
    if z>=0.65 and z<0.8:
        w = [0,0,0.66,1]
    if z>=0.8:
        w = [0,0,1,0.66]
    return w


####################input#########
RA= 65.7490   #17.8446#80.5159
DEC= -46.1436 #-55.3138#-50.4394
brick_lst = ["0657m462","0657m460"]#"0660m462","0654m460"]#["0176m552"]#["0804m505"]#["0229m135"]#["0036m305"]#,"2101p027","2103p027","2102p030"]#["1978m012"]#["3423m445","3413m445","3409m445"]
cluster = "SPT-CLJ0422-4608"#"SPT-CLJ0111-5518"#"SPT-CLJ0522-5026"
r200 = 2.79*u.arcmin#3.23*u.arcmin#3.51*u.arcmin
cols = ["release","brickid","brickname","objid","ra","dec","type","flux_g","flux_r","flux_i","flux_z","flux_ivar_g","flux_ivar_r","flux_ivar_i","flux_ivar_z","mw_transmission_g","mw_transmission_r","mw_transmission_i","mw_transmission_z"]
s = 2   #sigmas for sigmaclipping
r2cut = 0.5   #how many r200
##################################
cols = ['release','brickid','brickname','objid','brick_primary','ra','dec','ra_ivar','dec_ivar','type','flux_g','flux_r','flux_i','flux_z','flux_ivar_g','flux_ivar_r','flux_ivar_i','flux_ivar_z','mw_transmission_g','mw_transmission_r','mw_transmission_i','mw_transmission_z']


brick_lst =  ["0657m462","0657m460","0654m460"]
for i in range(len(brick_lst)):   #load catalogs
    brick = brick_lst[i]
    if i==0:
        crt = Table.read(cluster+"/"+brick+"/tractor/tractor-"+brick+".fits",format="fits")
        crt = crt[crt["brick_primary"]==False] 
        crt = crt[cols]
    else:
        icrt = Table.read(cluster+"/"+brick+"/tractor/tractor-"+brick+".fits",format="fits")
        icrt = icrt[icrt["brick_primary"]==False]
        icrt = icrt[cols]
        crt = vstack([crt,icrt])

with warnings.catch_warnings():  # Ignore warnings
    warnings.simplefilter('ignore')
    crt["m_g"] = 22.5-2.5*np.log10(crt["flux_g"]/crt["mw_transmission_g"])  #estimate magnitudes
    crt["m_r"] = 22.5-2.5*np.log10(crt["flux_r"]/crt["mw_transmission_r"])
    crt["m_i"] = 22.5-2.5*np.log10(crt["flux_i"]/crt["mw_transmission_i"])
    crt["m_z"] = 22.5-2.5*np.log10(crt["flux_z"]/crt["mw_transmission_z"])
    #crt["sig_mag_err"] = -2.5*np.log10((crt["flux_ivar_i"]**(-(1/2)))/crt["mw_transmission_i"])

crt = crt[np.isnan(crt["m_g"])==False]    #remove nan values
crt = crt[np.isnan(crt["m_r"])==False]
crt = crt[np.isnan(crt["m_i"])==False]
crt = crt[np.isnan(crt["m_z"])==False]

x_err = (crt["flux_ivar_g"]**(-(1/2)))/crt["mw_transmission_g"]   #estimate mag errors
x = crt["flux_g"]/crt["mw_transmission_g"]
crt["mag_g_err"] = (x_err/(x*np.log(10)))*2.5
x_err = (crt["flux_ivar_r"]**(-(1/2)))/crt["mw_transmission_r"]  
x = crt["flux_r"]/crt["mw_transmission_r"]
crt["mag_r_err"] = (x_err/(x*np.log(10)))*2.5
x_err = (crt["flux_ivar_i"]**(-(1/2)))/crt["mw_transmission_i"]
x = crt["flux_i"]/crt["mw_transmission_i"]
crt["mag_i_err"] = (x_err/(x*np.log(10)))*2.5 
x_err = (crt["flux_ivar_z"]**(-(1/2)))/crt["mw_transmission_z"]
x = crt["flux_z"]/crt["mw_transmission_z"]
crt["mag_z_err"] = (x_err/(x*np.log(10)))*2.5 

###################################################################

ra, dec = 65.6147,-46.12468
d = SkyCoord(ra*u.degree,dec*u.degree)     #find Mpc distance and pos angle of each galaxy from the cluster center
catalog = SkyCoord(crt["ra"],crt["dec"])
lst = d.separation(catalog).to(u.arcmin)
crt["PRJ_SEP"] = lst
#crt = crt[crt["PRJ_SEP"]<=r200]
point = crt[crt["PRJ_SEP"]<=1]

point = point[point["brickname"]!="0654m460"]

nextlst = []
for i in range(len(point)):
    ira = point["ra"][i]
    idec = point["dec"][i]
    d = SkyCoord(ira*u.degree,idec*u.degree)     #find Mpc distance and pos angle of each galaxy from the cluster center
    catalog = SkyCoord(point["ra"],point["dec"])
    lst = d.separation(catalog).to(u.arcsec)
    lst.sort()
    nextlst += [lst[1].value]    #take nearest neighbor
point["nearest"] = nextlst
point.sort("nearest")

plt.clf()   #see spot taken
plt.scatter(crt["ra"],crt["dec"])
plt.scatter(point["ra"],point["dec"],color="red")
plt.show()

lek = []         #separate by brickname
for i in range(len(brick_lst)):
    lek += [point[point["brickname"]==brick_lst[i]]]

plt.clf()    #see data overlap
#plt.scatter(crt["ra"],crt["dec"])
plt.scatter(lek[0]["ra"],lek[0]["dec"],edgecolor="black",marker="v",alpha=0.9,label=brick_lst[0])
plt.scatter(lek[1]["ra"],lek[1]["dec"],edgecolor="black",marker="s",alpha=0.6,label=brick_lst[1])
plt.scatter(lek[2]["ra"],lek[2]["dec"],edgecolor="black",marker="^",alpha=0.9,label=brick_lst[2])
plt.legend()
plt.show()

plt.clf()
plt.hist(point["nearest"],range=[0,10],bins=50)
plt.xlabel("nearest separation [arcsec]")
plt.show()

repoint = point[point["nearest"]<0.2]
repoint.sort("nearest")
repoint[:6].sort("ra")

first = repoint[::2]
second = repoint[1::2]
band = "r"
first["magoffset"] = first["m_"+band]-second["m_"+band]
second["magoffset"] = first["m_"+band]-second["m_"+band]


################################################
band = "r"
band2 = "i"
plt.clf()
plt.scatter(first["nearest"],first["m_"+band],marker="v",label=first["brickname"][0])
plt.scatter(second["nearest"],second["m_"+band],marker="^",label=second["brickname"][0])
plt.errorbar(first["nearest"],first["m_"+band],yerr=first["mag_"+band+"_err"],fmt="none",color="black",barsavobe=True,capsize=5,alpha=0.5)
plt.errorbar(second["nearest"],second["m_"+band],yerr=second["mag_"+band+"_err"],fmt="none",color="black",barsavobe=True,capsize=5,alpha=0.5)
plt.ylabel("mag "+band+"-band")
plt.xlabel("Nearest neighbor sep [arcsec]")
plt.gca().invert_yaxis()
plt.legend(loc="upper right",frameon=False)
plt.show()

plt.clf()
#plt.scatter(first["m_"+band],first["m_"+band]-second["m_"+band],marker="v",label="first - second")
plt.errorbar(first["m_"+band],first["m_"+band]-second["m_"+band],yerr=first["mag_"+band+"_err"],fmt="ro",ecolor="black",capsize=5)
plt.ylabel("m$_{"+band+",1}$ - m$_{"+band+",2}$")
plt.xlabel("m$_{"+band+"}$")
plt.text(18,1,s="rms = "+str((np.sum(x*x for x in first["m_"+band]-second["m_"+band])/len(first["m_"+band]-second["m_"+band]))**(1/2)))
#plt.gca().invert_yaxis()
#plt.legend(loc="upper left",frameon=False)
plt.show()


plt.clf()
plt.scatter(np.linspace(1,len(first)+1,len(first)),first["m_r"],marker="v",label="first")
plt.scatter(np.linspace(1,len(second)+1,len(second)),second["m_r"],marker="^",label="second")
plt.errorbar(np.linspace(1,len(first)+1,len(first)),first["m_r"],yerr=first["mag_r_err"],fmt="none",color="black",barsavobe=True,capsize=5,alpha=0.5)
plt.errorbar(np.linspace(1,len(second)+1,len(second)),second["m_r"],yerr=second["mag_r_err"],fmt="none",color="black",barsavobe=True,capsize=5,alpha=0.5)
plt.ylabel("mag r-band")
plt.gca().invert_yaxis()
plt.legend(loc="lower left")
plt.show()


plt.clf()    #see data overlap
#plt.scatter(crt["ra"],crt["dec"])
plt.scatter(lek[0]["ra"],lek[0]["dec"],edgecolor="black",facecolor="none",marker="s",label=brick_lst[0])
plt.scatter(lek[1]["ra"],lek[1]["dec"],color="black",marker="+",label=brick_lst[1])
plt.scatter(lek[2]["ra"],lek[2]["dec"],color="black",marker="x",label=brick_lst[2])
plt.errorbar(lek[0]["ra"],lek[0]["dec"],xerr=lek[0]["ra_ivar"]**(-2),yerr=lek[0]["dec_ivar"]**(-2),fmt="none",color="red",barsavobe=True,capsize=0.0,alpha=0.5)
plt.errorbar(lek[1]["ra"],lek[1]["dec"],xerr=lek[1]["ra_ivar"]**(-2),yerr=lek[1]["dec_ivar"]**(-2),fmt="none",color="red",barsavobe=True,capsize=0.0,alpha=0.5)
plt.errorbar(lek[2]["ra"],lek[2]["dec"],xerr=lek[2]["ra_ivar"]**(-2),yerr=lek[2]["dec_ivar"]**(-2),fmt="none",color="red",barsavobe=True,capsize=0.0,alpha=0.5)
plt.legend()
plt.show()





















