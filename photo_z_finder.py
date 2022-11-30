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
    #w = [0.0,0.0,1.0,0.0]

    #if z<0.35:
    #    w = [1.33,0.33,0,0.33]
    #if z>=0.35 and z<0.75:
    #    w = [0,0.66,0.25,0.875]
    #if z>=0.75:
    #    w = [0,0.33,1,0.66]

    if z<0.36:
        w = [1.33,0.33,0,0.33]
    if z>=0.36 and z<0.6:
        w = [0,1,0.25,1]
    if z>=0.6 and z<0.76:
        w = [0,0.5,0.25,1.25]
    if z>=0.76:
        w = [0,0,1.25,0.33]
    return w

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

bcgz20 =  Table.read("bayliss_bcg_radec.latex",format="latex")
####################input#########
init_RA= [65.7490,356.1481, 27.7898, 26.1795, 90.0614, 359.7075, 72.9661, 58.5612, 69.929, 54.4573, 17.8446, 23.9753, 80.5159, 315.0969, 93.0249, 87.5504]
init_DEC= [-46.1436,-42.41, -56.911, -48.1281, -43.8879, -61.4862, -49.8796, -59.0733, -53.5038, -49.4738, -55.3138, -59.0814, -50.4394, -57.1643, -43.2992, -50.3236]
init_brick_lst = [["0657m462","0657m460"],["3560m427","3561m422","3561m425","3564m422","3564m425","3564m427"],["0276m570","0278m567","0278m572","0280m570"],["0259m480","0260m482","0262m480","0263m482"],["0898m437","0898m440","0900m442","0901m437","0901m440","0903m442"],["3597m612","3597m615","3597m617"],["0729m500","0730m497","0730m502","0733m500"],["0586m590","0586m592"],["0696m532","0696m535","0696m537","0700m532","0700m537","0701m535"],["0543m492","0543m495","0545m497","0546m495","0547m492"],["0175m550","0176m552","0177m555","0180m550","0181m552"],["0238m590","0240m592","0243m590"],["0804m502","0804m505","0805m507","0808m502","0808m505"],["3148m575","3150m570","3151m572","3153m575","3154m570"],["0928m430","0929m432","0929m435","0932m430","0932m432","0932m435"],["0874m502","0874m505","0876m500","0878m502","0878m505"]]
init_cluster = ['SPT-CLJ0422-4608','SPT-CLJ2344-4224', 'SPT-CLJ0151-5654', 'SPT-CLJ0144-4807', 'SPT-CLJ0600-4353', 'SPT-CLJ2358-6129', 'SPT-CLJ0451-4952', 'SPT-CLJ0354-5904', 'SPT-CLJ0439-5330', 'SPT-CLJ0337-4928', 'SPT-CLJ0111-5518', 'SPT-CLJ0135-5904', 'SPT-CLJ0522-5026', 'SPT-CLJ2100-5708', 'SPT-CLJ0612-4317', 'SPT-CLJ0550-5019']
init_r200 = [2.79,5.44,5.54,5.27, 5.4,4.92,4.34,4.61,4.23,3.59,3.23,3.57,3.51, 3.6,3.73, 2.9]*u.arcmin
cols = ["release","brickid","brickname","objid","ra","dec","type","flux_g","flux_r","flux_i","flux_z","flux_ivar_g","flux_ivar_r","flux_ivar_i","flux_ivar_z","mw_transmission_g","mw_transmission_r","mw_transmission_i","mw_transmission_z"]
s = 2.0   #sigmas for sigmaclipping
r2cut = 0.5   #how many r200
##################################
pz_lst = []
for RA, DEC, brick_lst, cluster, r200 in zip(init_RA,init_DEC,init_brick_lst,init_cluster,init_r200):
    print("============"+cluster+"=============")
    for i in range(len(brick_lst)):   #load catalogs
        brick = brick_lst[i]
        if i==0:
            crt = Table.read(cluster+"/"+brick+"/tractor/tractor-"+brick+".fits",format="fits")
            crt = crt[crt["brick_primary"]==True]
            crt = crt[cols]
        else:
            icrt = Table.read(cluster+"/"+brick+"/tractor/tractor-"+brick+".fits",format="fits")
            icrt = icrt[icrt["brick_primary"]==True]
            icrt = icrt[cols]
            crt = vstack([crt,icrt])

    ########################estimate_magnitudes########################
    with warnings.catch_warnings():  # Ignore warnings
        warnings.simplefilter('ignore')
        crt["m_g"] = 22.5-2.5*np.log10(crt["flux_g"]/crt["mw_transmission_g"])  #estimate magnitudes
        crt["m_r"] = 22.5-2.5*np.log10(crt["flux_r"]/crt["mw_transmission_r"])
        crt["m_i"] = 22.5-2.5*np.log10(crt["flux_i"]/crt["mw_transmission_i"])
        crt["m_z"] = 22.5-2.5*np.log10(crt["flux_z"]/crt["mw_transmission_z"])
        #crt["sig_mag_err"] = -2.5*np.log10((crt["flux_ivar_i"]**(-(1/2)))/crt["mw_transmission_i"])

        #crt["m_g"] = (22.5-2.5*np.log10(crt["flux_g"]))#*crt["mw_transmission_r"]  #estimate magnitudes
        #crt["m_r"] = (22.5-2.5*np.log10(crt["flux_r"]))#*crt["mw_transmission_r"]
        #crt["m_i"] = 22.5-2.5*np.log10(crt["flux_i"])
        #crt["m_z"] = 22.5-2.5*np.log10(crt["flux_z"])
        
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
    ################cut#####################
    ra, dec = RA, DEC
    d = SkyCoord(ra*u.degree,dec*u.degree)     #find Mpc distance and pos angle of each galaxy from the cluster center
    catalog = SkyCoord(crt["ra"],crt["dec"])
    lst = d.separation(catalog).to(u.arcmin)
    crt["PRJ_SEP"] = lst
    crt = crt[crt["type"]!="PSF"]      #remove stars
    crt = crt[crt["type"]!="DUP"]      #remove gaia sources
    grcrt = crt.copy()
    gcrt = crt[crt["PRJ_SEP"]<=r200] 
    crt = crt[crt["PRJ_SEP"]<=r2cut*r200]    #0.5r200 cut --> hennig17
    #crt = crt[crt["flux_ivar_i"]**(-(1/2))<10**(0.1/2.5)]     #sig_mag < 0.1 cut
    crt["N_ID"] = np.linspace(0,len(crt)-1,len(crt),dtype="int")   #give identifier label
    #########ds9reg#########
    t = open(cluster+"/"+cluster+".reg","w")
    t.write("# Region file format: DS9 version 4.1\n")
    t.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
    t.write("fk5\n")
    for i in range(len(crt)):
        t.write("circle("+str(crt["ra"][i])+","+str(crt["dec"][i])+",7.920\")\n")
    t.close()
    pth = cluster+"/"+brick_lst[0]+"/coadd/"+brick_lst[0]+"/"
    bnum = 0
    print("#----------"+str(r2cut)+"R200 data / ds9 region")
    print("ds9 -rgb -red "+pth+"legacysurvey-"+brick_lst[bnum]+"-image-i.fits.fz -linear -scale 99.5 -blue "+pth+"legacysurvey-"+brick_lst[bnum]+"-image-g.fits.fz -linear -scale 99.5 -green "+pth+"legacysurvey-"+brick_lst[bnum]+"-image-r.fits.fz -linear -scale 99.5 -region "+cluster+"/"+cluster+".reg &")
    ########################
    cmr = Table.read("/home/n0d3j1tqu0/Programs/easyGalaxy/cmr_z_DES_griz_zf3_0p4decay_chab_new",format="ascii")
    cmr = cmr[1:]   #remove some nan values in the model
    #cmr = cmr[:290]
    cmr = cmr[:100]
    zcl_lst = []
    for i in range(len(cmr)):
        gcol = list(cmr[cmr.colnames[2::6]][i])    #extract g, r, i, z magnitudes at 6 different levels 3, 2, 1, 0.5, 0.4, 0.3
        rcol = list(cmr[cmr.colnames[3::6]][i])
        icol = list(cmr[cmr.colnames[4::6]][i])
        zcol = list(cmr[cmr.colnames[5::6]][i])
     
        mr = [rcol,icol,zcol,zcol,zcol]
        mb = [gcol,rcol,icol,rcol,gcol]
        zlcut = [cmr["col16"][i],cmr["col17"][i],cmr["col18"][i],cmr["col18"][i]]#,cmr["col18"][i]]
        zrlab = ["r","i","z","z"]#,"z"]
        zblab = ["g","r","i","r"]#,"g"]
        res_lst = []
        lab_lst = []
        sig_lst = []
        for mredder, mbluer, lcut, rlab, blab in zip(mr,mb,zlcut,zrlab,zblab):
            lcrt = crt[crt["m_"+rlab]<=lcut+2]
            lcrt = lcrt[lcrt["m_"+rlab]>=lcut-4]
            lcrt = lcrt[lcrt["mag_"+rlab+"_err"]<0.1]    #sig_mag < 0.1 hennig17
            redder = np.array(lcrt["m_"+rlab])
            bluer = np.array(lcrt["m_"+blab])
            rederr = np.array(lcrt["mag_"+rlab+"_err"])
            bluerr = np.array(lcrt["mag_"+blab+"_err"])
        
            if len(lcrt)<10:
                if len(lcrt)!=0:
                    lcrt = Table(lcrt[0])
                    lcrt.remove_row(0)
                sig_lst += [0.0]
                res_lst += [lcrt]
                lab_lst += [blab+rlab]
                continue

            #make rs model line by fitting to the six points
            params = lmodel.make_params(A=((mbluer[2]-mredder[2])-(mbluer[1]-mredder[1]))/(mredder[2]-mredder[1]),B=0) #initial params
            result = lmodel.fit(np.array(mbluer)-np.array(mredder), params, x=mredder)
            Ares = np.array(result.params)[0]
            Bres = np.array(result.params)[1]
            real_yrcs = bluer-redder
            model_yrcs = linear(redder,Ares,Bres)
            md_list = model_dist(Ares,Bres,redder,real_yrcs)    #take ortogonal distance to model
            devs = np.array(md_list)
            theta = np.arctan(Ares)
            sig_proj2 = (np.sqrt(bluerr**2 + rederr**2)*np.cos(theta))**2 + (rederr*np.sin(theta))**2
            sig_int2 = 0.03**2            #hennig17
            sig_col2 = sig_proj2 + sig_int2   #estimate ortogonal sigma2
            clip = Table([lcrt["N_ID"],devs,np.sqrt(sig_proj2),real_yrcs,model_yrcs],names=["N_ID","devs","devs_err","real_yrcs","model_yrcs"])
            clip = clip[abs(clip["devs"])<=0.22]      #work only with data within +-0.22 from model  (lopez-cruz##)
            if len(clip)<10:
                if len(clip)!=0:
                    clip = Table(clip[0])
                    clip.remove_row(0)
                sig_lst += [0.0]
                res_lst += [clip]
                lab_lst += [blab+rlab]
                continue
            #######################testing_model###########################
            #---initial_biweight_iteration---#
            rcs = 0    #initial values (these values are not too important)
            sig = np.mean(np.sqrt(sig_col2))  
            for l in range(3):
                rcs = st.biweight_location(devs,M=np.array([rcs]))    #improve initial guess before 3sigma clipping (beers et al 1990)
            cont = 0
            #---3sigmaclipping
            while True:
                try:
                    rcs, sig = weighted_avg_and_std(clip["devs"],weights=clip["devs_err"]**(-1)) 
                except:
                    pass
                #rcs = st.biweight_location(clip["devs"],M=np.array([rcs]))  
                #sig = st.biweight_scale(clip["devs"],M=np.array([sig]))   #---Sigma Biweight
                rej = clip[abs(clip["devs"]-rcs) >= (s*sig)]   #----rejected_galaxies
                if len(rej) == 0:      #---when there are no rejected galaxies end the iteration and calculate sigma error
                    #print("NUMBER OF ITERATIONS: "+str(cont))
                    break
                clip = clip[abs(clip["devs"]-rcs) < (s*sig)] #---if len(rej)!=0 then cut in s*sigma and start over again
                cont += 1 
            if len(clip)<10:
                if len(clip)!=0: 
                    clip = Table(clip[0])
                    clip.remove_row(0)
                sig_lst += [0.0]
                res_lst += [clip]
                lab_lst += [blab+rlab]
                continue
            sig_lst += [sig]
            res_lst += [clip]
            lab_lst += [blab+rlab]
            ######################################################
        #---chi2
        chi_lst = []
        for j in range(len(res_lst)):
            clip = res_lst[j]
            if len(clip)<10:
                chi_lst += [[0.0,lab_lst[j]]]    #add chi2 0.0 if len(clip)<10
                continue
            c2 = np.sum(((clip["real_yrcs"]-clip["model_yrcs"])**2)/clip["model_yrcs"])/(len(clip)-1-2)
            chi_lst += [[c2,lab_lst[j]]]
        chi_lst = Table(np.array(chi_lst),names=["chi2","band"])
        chi_lst["chi2"] = Table.Column(chi_lst["chi2"],dtype="float")
        cweight = col_evo_weights(cmr["col1"][i])
        while len(cweight)<len(chi_lst):    #add weight 0 for other filters ()
            cweight += [0]
        tablon = Table([cweight,list(chi_lst["band"])],names=["cweight","band"])
        chitab = join(chi_lst,tablon,keys="band") 
        chitab = chitab[chitab["chi2"]!=0.0]     #chi2=0.0 are len(clip)<10, not interested
        if len(chitab)!=0:    #for bands with estimated chi2, get weighted average
            if np.sum(chitab["cweight"])==0.0:
                chi2 = 0.0
            else:
                chi2 = np.average(chitab["chi2"],weights=chitab["cweight"])    
        else:    #if no bands with chi2, then make chi2=0.0
            chi2 = 0.0
        
        #---Guille-thing
        #york = ["gr","ri","iz","rz","gz"]
        #grand_lst = []
        #for j in range(len(crt)):
        #    idi = crt["N_ID"][j]
        #    lst = []
        #    plox = ''
        #    for r in range(len(res_lst)):
        #        resrow = res_lst[r][res_lst[r]["N_ID"]==idi]
        #        if len(resrow)==0:
        #            lst += [99]
        #            continue
        #        plox += york[r]
        #        lst += [(resrow["model_yrcs"][0]-resrow["real_yrcs"][0])**2]
        #    grand_lst += [[idi,np.sqrt(np.sum(lst)),plox]]

        #grand_tab = Table(np.array(grand_lst),names=["N_ID","res_sum","bands"])
        #grand_tab["N_ID"] = Table.Column(grand_tab["N_ID"],dtype="int")
        #grand_tab["res_sum"] = Table.Column(grand_tab["res_sum"],dtype="float")
        #glakoon = join(crt,grand_tab,"N_ID")
        #thing = np.sum(glakoon["res_sum"]**2)/len(glakoon)-6
        ########################################################
        zcl_lst += [[cmr["col1"][i],chi2,]]
    ztab = Table(np.reshape(zcl_lst,[len(zcl_lst),2]))

    plt.clf()
    iztab = ztab[ztab["col1"]!=99.0]
    plt.scatter(iztab["col0"],iztab["col1"],s=5,label="combined chi2")
    #iztab = ztab[ztab["col2"]!=99.0]
    #plt.scatter(iztab["col0"],iztab["col2"],s=5,label="r-i")
    #iztab = ztab[ztab["col3"]!=99.0]
    #plt.scatter(iztab["col0"],iztab["col3"],s=5,label="i-z")
    #iztab = ztab[ztab["col4"]!=99.0]
    #plt.scatter(iztab["col0"],iztab["col4"],s=5,label="r-z")
    #iztab = ztab[ztab["col5"]!=99.0]
    #plt.scatter(iztab["col0"],iztab["col5"],s=5,label="g-z")
    plt.xlabel("z")
    plt.ylabel(r"$\chi^{2}$/DOF")
    plt.legend()
    plt.savefig(cluster+"/"+brick_lst[0]+"_chi2.jpg")
    plt.close()
    #plt.show()
    
    zt = ztab[ztab["col1"]!=0.0]
    pz = zt[zt["col1"]==np.min(zt["col1"])]["col0"][0]
    print("photo-z: "+str(pz))

    #######################SELECT-RS#######################
    wid = col_evo_weights(pz)    #select proper band for this redshift
    band = ["gr","ri","iz","rz","gz"][wid.index(np.max(wid))]
    rlab = band[1]
    blab = band[0]

    i = list(cmr["col1"]).index(pz)
    gcol = list(cmr[cmr.colnames[2::6]][i])    #extract g, r, i, z magnitudes at 6 different levels 3, 2, 1, 0.5, 0.4, 0.3
    rcol = list(cmr[cmr.colnames[3::6]][i])
    icol = list(cmr[cmr.colnames[4::6]][i])
    zcol = list(cmr[cmr.colnames[5::6]][i])

    mr = [rcol,icol,zcol,zcol,zcol]
    mb = [gcol,rcol,icol,rcol,gcol]
    zlcut = [cmr["col16"][i],cmr["col17"][i],cmr["col18"][i],cmr["col18"][i],cmr["col18"][i]]
    zrlab = ["r","i","z","z","z"]
    zblab = ["g","r","i","r","g"]
    res_lst = []
    lab_lst = []
    for mredder, mbluer, lcut, rlab, blab in zip(mr,mb,zlcut,zrlab,zblab):
        lcrt = crt[crt["m_"+rlab]<=lcut+2]
        lcrt = lcrt[lcrt["m_"+rlab]>=lcut-4]
        lcrt = lcrt[lcrt["mag_"+rlab+"_err"]<0.1]    #sig_mag < 0.1 hennig17
        redder = np.array(lcrt["m_"+rlab])
        bluer = np.array(lcrt["m_"+blab])
        rederr = np.array(lcrt["flux_ivar_"+rlab]**(-(1/2)))
        bluerr = np.array(lcrt["flux_ivar_"+blab]**(-(1/2)))
        if len(lcrt)<10:
            if len(lcrt)!=0:
                lcrt = Table(lcrt[0])
                lcrt.remove_row(0)
            res_lst += [lcrt]
            lab_lst += [blab+rlab]
            continue
        #make rs model line by fitting to the six points
        params = lmodel.make_params(A=((mbluer[2]-mredder[2])-(mbluer[1]-mredder[1]))/(mredder[2]-mredder[1]),B=0) #initial params
        result = lmodel.fit(np.array(mbluer)-np.array(mredder), params, x=mredder)
        Ares = np.array(result.params)[0]
        Bres = np.array(result.params)[1]
        real_yrcs = bluer-redder
        model_yrcs = linear(redder,Ares,Bres)
        md_list = model_dist(Ares,Bres,redder,real_yrcs)    #take ortogonal distance to model
        devs = np.array(md_list)
        theta = np.arctan(Ares)
        sig_proj2 = (np.sqrt(bluerr**2 + rederr**2)*np.cos(theta))**2 + (rederr*np.sin(theta))**2
        sig_int2 = 0.03**2            #hennig17
        sig_col2 = sig_proj2 + sig_int2   #estimate ortogonal sigma2
        clip = Table([lcrt["N_ID"],devs,np.sqrt(sig_proj2),redder,real_yrcs,model_yrcs],names=["N_ID","devs","devs_err","redder","real_yrcs","model_yrcs"])
        clip = clip[abs(clip["devs"])<=0.22]      #work only with data within +-0.22 from model  (lopez-cruz##)
        if len(clip)<10:
            if len(clip)!=0:
                clip = Table(clip[0])
                clip.remove_row(0)
            res_lst += [clip]
            lab_lst += [blab+rlab]
            continue
        #######################testing_model###########################
        #---initial_biweight_iteration---#
        rcs = 0    #initial values (these values are not too important)
        sig = np.mean(np.sqrt(sig_col2))
        for l in range(3):
            rcs = st.biweight_location(devs,M=np.array([rcs]))    #improve initial guess before 3sigma clipping (beers et al 1990)
        cont = 0
        #---3sigmaclipping
        while True:
            try:
                rcs, sig = weighted_avg_and_std(clip["devs"],weights=clip["devs_err"]**(-1))
            except:
                pass
            #rcs = st.biweight_location(clip["devs"],M=np.array([rcs]))  
            #sig = st.biweight_scale(clip["devs"],M=np.array([sig]))   #---Sigma Biweight
            rej = clip[abs(clip["devs"]-rcs) >= (s*sig)]   #----rejected_galaxies
            if len(rej) == 0:      #---when there are no rejected galaxies end the iteration and calculate sigma error
                #print("NUMBER OF ITERATIONS: "+str(cont))
                break
            clip = clip[abs(clip["devs"]-rcs) < (s*sig)] #---if len(rej)!=0 then cut in s*sigma and start over again
            cont += 1
        if len(clip)<10:
            if len(clip)!=0:
                clip = Table(clip[0])
                clip.remove_row(0)
            res_lst += [clip]
            lab_lst += [blab+rlab]
            continue
        res_lst += [clip]
        lab_lst += [blab+rlab]

    #------------------------------------get clipped data used to get the z-model
    wid = col_evo_weights(pz)
    clip = res_lst[wid.index(np.max(wid))]
    band = ["gr","ri","iz","rz","gz"][wid.index(np.max(wid))]
    rlab = band[1]
    blab = band[0]
    lcut = [cmr["col16"][i],cmr["col17"][i],cmr["col18"][i],cmr["col18"][i],cmr["col18"][i]][wid.index(np.max(wid))]
    lcrt = crt[crt["m_"+rlab]<=lcut+2]
    lcrt = lcrt[lcrt["m_"+rlab]>=lcut-4]
    lcrt = lcrt[lcrt["mag_"+rlab+"_err"]<0.1]    #sig_mag < 0.1 hennig17
    redder = np.array(lcrt["m_"+rlab])
    bluer = np.array(lcrt["m_"+blab])
    mredder = [rcol,icol,zcol,zcol,zcol][wid.index(np.max(wid))]
    mbluer = [gcol,rcol,icol,rcol,gcol][wid.index(np.max(wid))]
    params = lmodel.make_params(A=((mbluer[2]-mredder[2])-(mbluer[1]-mredder[1]))/(mredder[2]-mredder[1]),B=0) #initial params
    result = lmodel.fit(np.array(mbluer)-np.array(mredder), params, x=mredder)
    Ares = np.array(result.params)[0]
    Bres = np.array(result.params)[1]
    #-------------------------------------

    #########################PLOTTING###########################
    #plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['xtick.major.size'] = 7
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['ytick.major.size'] = 7
    plt.rcParams['ytick.minor.size'] = 3

    plt.clf()     #plotting _zmodel.jpg
    fig, ax = plt.subplots()#subplot_kw={'aspect': 'equal'})
    ax.scatter(crt["m_"+rlab],crt["m_"+blab]-crt["m_"+rlab],s=5,c="black",label=str(r2cut)+"R200 galaxies")
    ax.scatter(clip["redder"],clip["real_yrcs"],s=5,c="red",label="clipped RS")
    #ax.scatter(bcg_row["MOF_BDF_MAG_I_CORRECTED"],bcg_row["MOF_BDF_MAG_R_CORRECTED"]-bcg_row["MOF_BDF_MAG_I_CORRECTED"],s=70,marker="^",color="none",edgecolor="black",label="BCG",zorder=5)
    x = np.linspace(min(redder),max(redder),len(redder))
    ax.plot(x,linear(x,Ares,Bres),color="red",label="RED SEQUENCE MODEL, z="+str(pz))
    ax.plot(x,linear(x,Ares,Bres)+0.22,ls="--",color="red")
    ax.plot(x,linear(x,Ares,Bres)-0.22,ls="--",color="red")
    #red_gal = Table(crt[0])
    #red_gal.remove_row(0)
    #blue_gal = Table(crt[0])
    #blue_gal.remove_row(0)
    #for i in range(len(crt)):
    #    if crt["m_r"][i]-crt["m_i"][i] <= linear(crt["m_i"][i],Ares,Bres)+0.22 and crt["m_r"][i]-crt["m_i"][i] >= linear(crt["m_i"][i],Ares,Bres)-0.22:
    #        red_gal.add_row(crt[i])
    #    if crt["m_r"][i]-crt["m_i"][i] < linear(crt["m_i"][i],Ares,Bres)-0.22:
    #        blue_gal.add_row(crt[i])
    #ax.scatter(red_gal["m_i"],red_gal["m_r"]-red_gal["m_i"],s=14,marker="s",color="red",label="RED")
    #ax.scatter(blue_gal["m_i"],blue_gal["m_r"]-blue_gal["m_i"],s=12,marker="s",color="blue",label="BLUE")
    ax.set_ylim([-0.5,2])
    ax.set_xlim([15,27])
    ax.tick_params(top=True,right=True)
    ax.tick_params(which="minor",top=True,right=True)
    plt.xlabel(rlab)
    plt.ylabel(blab+"-"+rlab)
    plt.title(cluster)
    plt.minorticks_on()
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(cluster+"/"+brick_lst[0]+"_zmodel.jpg")
    plt.close()
    #plt.show()
    #############################################################
    #-------------------------------------get R200 data within +-0.22 mags from z-model
    gredder = np.array(gcrt["m_"+rlab])   
    gbluer = np.array(gcrt["m_"+blab])
    gmredder = [rcol,icol,zcol,zcol,zcol][wid.index(np.max(wid))]
    gmbluer = [gcol,rcol,icol,rcol,gcol][wid.index(np.max(wid))]
    params = lmodel.make_params(A=((gmbluer[2]-gmredder[2])-(gmbluer[1]-gmredder[1]))/(gmredder[2]-gmredder[1]),B=0) #initial params
    result = lmodel.fit(np.array(gmbluer)-np.array(gmredder), params, x=gmredder)
    gAres = np.array(result.params)[0]
    gBres = np.array(result.params)[1]
    greal_yrcs = gbluer-gredder
    gmodel_yrcs = linear(gredder,gAres,gBres)
    md_list = model_dist(gAres,gBres,gredder,greal_yrcs)    #take ortogonal distance to model
    devs = np.array(md_list)
    gcrt["devs"] = devs
    redgal = gcrt[abs(gcrt["devs"])<=0.22]
    redgal.sort("m_"+rlab)                  #select mark BCGs as true
    truecat = np.zeros(len(redgal))
    truecat[:2] = True
    truecat[2:] = False
    redgal["is_bcg"] = Table.Column(truecat,dtype="bool")    #create column denoting BCG
    redgal.write(cluster+"/"+cluster+"_redsequence.cat",format="ascii",overwrite=True)   #save table
    #########ds9reg#########
    t = open(cluster+"/"+cluster+"_bcg.reg","w")
    t.write("# Region file format: DS9 version 4.1\n")
    t.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
    t.write("fk5\n")
    t.write("circle("+str(RA)+","+str(dec)+","+str(r200.to(u.arcsec).value)+"\") # color=red text={R200}\n")
    for r in range(len(redgal)):
        if redgal["is_bcg"][r]==True:
            t.write("circle("+str(redgal["ra"][r])+","+str(redgal["dec"][r])+",5.000\") # color=red text={m_"+rlab+" = "+str(np.around(redgal["m_"+rlab][r],2))+"}\n")
        else: 
            t.write("circle("+str(redgal["ra"][r])+","+str(redgal["dec"][r])+",5.000\")\n")
    t.close()
    pth = cluster+"/"+brick_lst[0]+"/coadd/"+brick_lst[0]+"/"
    bnum = 0
    print("#----------R200 data with BCGs / ds9 region")
    print("ds9 -rgb -red "+pth+"legacysurvey-"+brick_lst[bnum]+"-image-i.fits.fz -linear -scale 99.5 -blue "+pth+"legacysurvey-"+brick_lst[bnum]+"-image-g.fits.fz -linear -scale 99.5 -green "+pth+"legacysurvey-"+brick_lst[bnum]+"-image-r.fits.fz -linear -scale 99.5 -region "+cluster+"/"+cluster+"_bcg.reg &")
    ########################
    ##########################PLOTTING############################
    #---------------------------------------
    plt.clf()     #plotting _cmd.jpg
    fig, ax = plt.subplots()#subplot_kw={'aspect': 'equal'})
    ax.scatter(gcrt["m_"+rlab],gcrt["m_"+blab]-gcrt["m_"+rlab],s=5,c="black",label="R200 galaxies")
    ax.scatter(redgal["m_"+rlab],redgal["m_"+blab]-redgal["m_"+rlab],s=5,c="red",label="RED GALAXIES")
    bcgtab = redgal[redgal["is_bcg"]==True]
    bcgtab.sort("m_"+rlab)
    for r in range(len(bcgtab)):
        ax.scatter(bcgtab["m_"+rlab][r],bcgtab["m_"+blab][r]-bcgtab["m_"+rlab][r],s=20,c="none",edgecolor="black",marker=["^","v"][r],label="BCG "+str(r+1))
    #-------15SPT
    if cluster!="SPT-CLJ0422-4608":
        bz20 = bcgz20[bcgz20["SPT-CL"]==cluster.split("L")[1]] 
        ra, dec = bz20["RA"][0], bz20["DEC"][0]
        d = SkyCoord(ra*u.degree,dec*u.degree)     #find Mpc distance and pos angle of each galaxy from the cluster center
        catalog = SkyCoord(grcrt["ra"],grcrt["dec"])
        lst = d.separation(catalog).to(u.arcmin)
        grcrt["PRJ_SEP"] = lst
        fbz20 = grcrt[grcrt["PRJ_SEP"]==np.min(grcrt["PRJ_SEP"])]
        ax.scatter(fbz20["m_"+rlab][0],fbz20["m_"+blab][0]-fbz20["m_"+rlab][0],s=20,c="none",edgecolor="yellow",marker="s",label="BCG Z20")
    #------------
    #ax.scatter(bcg_row["MOF_BDF_MAG_I_CORRECTED"],bcg_row["MOF_BDF_MAG_R_CORRECTED"]-bcg_row["MOF_BDF_MAG_I_CORRECTED"],s=70,marker="^",color="none",edgecolor="black",label="BCG",zorder=5)
    x = np.linspace(min(redgal["m_"+rlab]),max(redgal["m_"+rlab]),len(redgal))
    ax.plot(x,linear(x,gAres,gBres),color="red",label="RED SEQUENCE MODEL, z="+str(pz))
    ax.plot(x,linear(x,gAres,gBres)+0.22,ls="--",color="red")
    ax.plot(x,linear(x,gAres,gBres)-0.22,ls="--",color="red")
    ax.set_ylim([-0.5,2])
    ax.set_xlim([15,27])
    ax.tick_params(top=True,right=True)
    ax.tick_params(which="minor",top=True,right=True)
    plt.xlabel(rlab)
    plt.ylabel(blab+"-"+rlab)
    plt.title(cluster)
    plt.minorticks_on()
    plt.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(cluster+"/"+brick_lst[0]+"_cmd.jpg")
    plt.close() 
    #plt.show()
    ##############################################################
    pz_lst += [pz]

pztab = Table([init_cluster,pz_lst])
pztab.write("photo-z",format="ascii")

#crt = Table.read("Datatable.latex",format="latex")
#crt["col0"] = ("SPT-CL"+" SPT-CL".join(crt["SPT-CL"])).split(" ")
#crt = join(pztab,crt)
#crt.sort("z")
#crt.rename_column("col1","photo-z")
#plt.clf()
#ax_scatter = plt.axes([0.1,0.3,0.82,0.56])
#ax_scatter.plot(np.linspace(0,1),np.linspace(0,1),color="black")
#ax_scatter.scatter(crt["z"],crt["photo-z"])
##plt.errorbar(crt["z"],crt["photo-z"],xerr=crt["z_err"],fmt="none",color="black",capsize=1)
#ax_scatter.set_xlim([0.2,0.8])
#ax_scatter.set_ylim([0.2,0.8])
#ax_scatter.set_xlabel("spec-z")
#ax_scatter.set_ylabel("photo-z")
#ax_scatter.get_xaxis().set_visible(False)
#ax_res = plt.axes([0.1,0.1,0.82,0.2])
#ax_res.scatter(crt["z"],(crt["photo-z"]-crt["z"])/(1+crt["z"]))
#ax_res.plot(np.linspace(0,1,10),np.zeros(10),color="black")
#ax_res.set_xlim([0.2,0.8])
#ax_res.set_xlabel("spec-z")
#ax_res.set_ylabel(r"z$_{p} - $z$_{s}$/(1+z$_{s}$)",fontsize=8)
#ax_res.text(0.6,0.1,s=r"$\Delta$z = "+str(np.around(np.median((crt["photo-z"]-crt["z"])/(1+crt["z"])),6)))
#ax_res.tick_params(right=True,labelright=True,labelleft=False,left=False)
#plt.show()

for i in range(len(init_cluster)):
     cluster = init_cluster[i]
     if i == 0:
         bcgtab = Table.read(cluster+"/"+cluster+"_redsequence.cat",format="ascii")
         bcgtab = Table(bcgtab[bcgtab["is_bcg"]=="True"])
     else:
         ibcgtab = Table.read(cluster+"/"+cluster+"_redsequence.cat",format="ascii")
         ibcgtab = Table(ibcgtab[ibcgtab["is_bcg"]=="True"])
         for j in range(len(ibcgtab)):
             bcgtab.add_row(ibcgtab[j])

bcg1 = []
bcg2 = []
for i in range(len(pztab)):
    ibcg = bcgtab[2*i:(2*i+2)]
    band = ["gr","ri","iz","rz"][col_evo_weights(pztab["col1"][i]).index(np.max(col_evo_weights(pztab["col1"][i])))]
    rlab = band[1]
    ibcg.sort("m_"+rlab)
    bcg1 += [str(ibcg["ra"][0])+","+str(ibcg["dec"][0])]
    bcg2 += [str(ibcg["ra"][1])+","+str(ibcg["dec"][1])]
pztab["bcg1"] = bcg1
pztab["bcg2"] = bcg2

#-----15SPT
bcg = Table.read("bayliss_bcg_radec.latex",format="latex")
bcg["col0"] = ('SPT-CL'+' SPT-CL'.join(bcg["SPT-CL"])).split(" ")
kal = join(pztab,bcg)
kal.sort("col1")
kal.remove_column("SPT-CL")
lst = []
for i in range(len(kal)):
    catalog = SkyCoord(np.array(','.join(list(kal[i][kal.colnames[2:4]])).split(","),dtype="float")[::2]*u.deg, np.array(','.join(list(kal[i][kal.colnames[2:4]])).split(","),dtype="float")[1::2]*u.deg)
    d = SkyCoord(kal["RA"][i]*u.deg,kal["DEC"][i]*u.deg)
    d.separation(catalog)
    lst += [np.min(d.separation(catalog).to(u.arcsec)).value]
kal["offset"] = lst
#---------
len(kal[kal["offset"]<1])/len(kal)






