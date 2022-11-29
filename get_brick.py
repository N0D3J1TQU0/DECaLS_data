import os
import sys
import numpy as np
from paramiko import SSHClient  #pip install paramiko
from scp import SCPClient   #pip install scp

CNAME = ['SPT-CLJ2344-4224', 'SPT-CLJ0151-5654', 'SPT-CLJ0144-4807', 'SPT-CLJ0600-4353', 'SPT-CLJ2358-6129', 'SPT-CLJ0451-4952', 'SPT-CLJ0354-5904', 'SPT-CLJ0439-5330', 'SPT-CLJ0337-4928', 'SPT-CLJ0111-5518', 'SPT-CLJ0135-5904', 'SPT-CLJ0522-5026', 'SPT-CLJ2100-5708', 'SPT-CLJ0612-4317', 'SPT-CLJ0550-5019']
SPTBRICKS = [["3560m427","3561m422","3561m425","3564m422","3564m425","3564m427",],["0276m570","0278m567","0278m572","0280m570"],["0259m480","0260m482","0262m480","0263m482"]]

CRA = np.array([356.1481, 27.7898, 26.1795, 90.0614, 359.7075, 72.9661, 58.5612, 69.929, 54.4573, 17.8446, 23.9753, 80.5159, 315.0969, 93.0249, 87.5504])

CDEC = np.array([-42.41, -56.911, -48.1281, -43.8879, -61.4862, -49.8796, -59.0733, -53.5038, -49.4738, -55.3138, -59.0814, -50.4394, -57.1643, -43.2992, -50.3236])

#ispt = 2
######################################
USER="irisuser"
HOST="cori.nersc.gov"
PASS="irispassword"
size1=0.3
outdir_lst = CNAME
init_RA=np.array([CRA-size1,CRA+size1]).T#[150,210]   #at least one decimal, can be float or list
init_DEC=np.array([CDEC-size1,CDEC+size1]).T#[-30,-25]  #at least one decimal, can be float or list
PATH_lst=["/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/","/global/cscratch1/sd/comparat/dr10/","/global/cscratch1/sd/mxhf/dr10c-test/","/global/cscratch1/sd/jsnigula/dr10c/","/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/"]
#coadds_lst=["blobmodel","ccds","chi2","depth","galdepth","image","invvar","maskbits","model","nexp","psfsize","resid","wise","wisemodel","wiseresid"]
#brick_lst = ["0657m460","0658m462","0657m463"]
#brick_lst = ["0262m530","0266m530","0263m532","0261m535","0260m527","0265m527"]
#brick_lst = ["0660m462"]

image=False
######################################
ssh = SSHClient()   #connect to DB
ssh.load_system_host_keys()
MFA=str(input("MFA key: "))
ssh.connect(hostname=HOST, username=USER, password=PASS+MFA)

for outdir, RA, DEC in zip(outdir_lst,init_RA,init_DEC):
    print(outdir+" - RA"+str(RA)+" DEC"+str(DEC))
    RA = list(RA)
    DEC = list(DEC)
    #convert to brick string
    try:
        len(brick_lst)
        area=False
        if type(brick_lst)==str:
            gonakillmyself[0]
    except:
        brick_lst = "0"
        if np.shape((RA,DEC))!=(2,):
            area=True
        else:
            area=False

    if str(brick_lst)=="0" and np.shape((RA,DEC))==(2,):
        RAtxt = str(int(np.around(RA*10,0)))
        lst=''
        for i in range(4 - len(RAtxt)):
            lst+="0"
        RAtxt=lst+RAtxt
        DECtxt=str(int(np.around(DEC*10,0)))
        sign=DECtxt.split("-")
        if len(sign)==2:
            sign="m"
            sgn="-"
        else:
            sign="p"
            sgn="+"
        lst=''
        for i in range(3 - len(DECtxt.split(sgn)[-1])):
            lst+="0"
        DECtxt=lst+DECtxt.split(sgn)[-1]
        brick_lst = [RAtxt+sign+DECtxt]
    if str(brick_lst)=="0" and np.shape((RA,DEC))!=(2,):
        RAtxt_inf = str(int(np.around(RA[0]*10,0))) 
        RAtxt_sup = str(int(np.around(RA[1]*10,0)))
        RAtxt_lst = np.linspace(int(RAtxt_inf),int(RAtxt_sup),int(RAtxt_sup)-int(RAtxt_inf)+1)
        RAtxt_temp = []
        for r in range(len(RAtxt_lst)):
            RAtxt = str(int(RAtxt_lst[r]))
            lst=''
            for i in range(4 - len(RAtxt)):
                lst+="0"
            RAtxt=lst+RAtxt
            RAtxt_temp += [RAtxt]
        RAtxt_lst = RAtxt_temp
        DECtxt_inf = str(int(np.around(DEC[0]*10,0)))
        DECtxt_sup = str(int(np.around(DEC[1]*10,0)))
        DECtxt_lst = np.linspace(int(DECtxt_inf),int(DECtxt_sup),int(DECtxt_sup)-int(DECtxt_inf)+1)
        DECtxt_temp = []
        for r in range(len(DECtxt_lst)):
            DECtxt = str(int(DECtxt_lst[r]))
            sign=DECtxt.split("-")
            if len(sign)==2:
                sign="m"
                sgn="-"
            else:
                sign="p"
                sgn="+"
            lst=''
            for i in range(3 - len(DECtxt.split(sgn)[-1])):
                lst+="0"
            DECtxt=sign+lst+DECtxt.split(sgn)[-1]
            DECtxt_temp += [DECtxt]
        DECtxt_lst = DECtxt_temp    
        brick_lst = []
        for i in range(len(RAtxt_lst)):
            for j in range(len(DECtxt_lst)):
                brick_lst += [RAtxt_lst[i]+DECtxt_lst[j]]

    #ssh = SSHClient()   #connect to DB
    #ssh.load_system_host_keys()
    #MFA=str(input("MFA key: "))
    #ssh.connect(hostname=HOST, username=USER, password=PASS+MFA)

    if area:     #locate available bricks in range
        ratoxt = np.linspace(int(RAtxt_inf),int(RAtxt_sup),int(RAtxt_sup)-int(RAtxt_inf)+1)
        #np.linspace(int(RAtxt_inf),int(RAtxt_sup),int((int(RAtxt_sup)-int(RAtxt_inf))/10)+1)/10
        diril = []
        print("Locating bricks...")
        for i in range(len(ratoxt)):
            ratxt = str(int(ratoxt[i]/10))
            while len(ratxt)<3:
                ratxt = "0"+ratxt
            for PATH in PATH_lst:
                DECtxt_lst.sort()
                stdin, stdout, stderr = ssh.exec_command("ls "+PATH+"tractor/"+ratxt+"/tractor*"+sign+"["+DECtxt_lst[0][1]+"-"+DECtxt_lst[-1][1]+"]*")   #make ls in the RA directory
                diri = str(stdout.read())[2:-1]   #remove b' at start and ' at end of string
                diri = diri.split("\\n")
                diri = diri[:-1]    #empty element at the end caused by split
                if len(diri)==0 or diri[0]=='':
                    continue
                for d in diri:
                    diril += ["tractor-"+d.split("/tractor-")[-1]]
        print("Crossmatching...")
        brck_tmp = []
        for brck in brick_lst:
            for inbrck in diril:
                if ".fits" in inbrck and brck in inbrck:
                    brck_tmp += [brck]
        brick_lst = list(np.unique(brck_tmp))

    try:                  #create directory to save data
        os.makedirs(outdir)
    except:
        # directory already exists
        pass

    t = open(outdir+"/bricks_info.txt",'w')
    t.write("LIST OF BIRCKS:\n")
    for ibrick in brick_lst:
        print("Input brick set to: "+ibrick)
        RAtxt=ibrick[:4]
        DECtxt=ibrick[5:]
        sign=ibrick[4]
        abort=0
        for r in range(len(PATH_lst)):
            PATH = PATH_lst[r]
            #search nearest brick
            stdin, stdout, stderr = ssh.exec_command("ls "+PATH+"tractor/"+RAtxt[:-1])   #make ls in the RA directory
            diri = str(stdout.read()).split("\\n")
            diri[0] = diri[0][2:]    #first element has a "b'", remove it 
            if len(brick_lst)>1:
                cont = 0
                for l in range(len(diri)):
                    if ibrick in diri[l]:
                        print("Brick found in "+PATH)
                        cont = 1
                        break
                if cont == 1:
                    nbrick = ibrick[:]
                    break
                print("No brick found in "+PATH)
                
            else:
                cont = 0
                if len(str(int(RAtxt)-10))<4:    #correct zeros before number
                    lst=''
                    for i in range(4 - len(str(int(RAtxt)-10))):
                        lst+="0"
                    RAtxt1=lst+str(int(RAtxt)-10)
                else:
                    RAtxt1 = str(int(RAtxt)-10)
                if RAtxt[:-1]=="000":  #exception in case we are in inferior limit. Just repeat RAtxt
                    RAtxt1 = RAtxt
                stdin, stdout, stderr = ssh.exec_command("ls "+PATH+"tractor/"+RAtxt1[:-1])   #make ls in the RA-1 directory
                diri1 = str(stdout.read()).split("\\n")
                diri1[0] = diri1[0][2:]    #first element has a "b'", remove it 

                if len(str(int(RAtxt)+10))<4:   #correct zeros before number
                    lst=''
                    for i in range(4 - len(str(int(RAtxt)+10))):
                        lst+="0"
                    RAtxt2=lst+str(int(RAtxt)+10)
                else:
                    RAtxt2 = str(int(RAtxt)+10)
                if RAtxt[:-1]=="359":  #exception in case we are in superior limit. Just repeat RAtxt
                    RAtxt2 = RAtxt
                stdin, stdout, stderr = ssh.exec_command("ls "+PATH+"tractor/"+RAtxt2[:-1])   #make ls in the RA-1 directory
                diri2 = str(stdout.read()).split("\\n")
                diri2[0] = diri2[0][2:]    #first element has a "b'", remove it
                
                diritot = diri+diri1+diri2   #im going to search also in near RA directories
                diril = []    #look tractor files
                for i in range(len(diritot)):
                    if sign in diritot[i] and "tractor" in diritot[i]:
                        diril += [diritot[i]]
                if len(diril)==0:
                    print("No brick found in "+PATH)
                    cont = 0
                    continue
                lst = []      #take ra and dec from bricks in diril and take the nearest to RAtxt,DECtxt
                for i in range(len(diril)):
                    ratxt = diril[i].split("-")[1].split(sign)[0]
                    dectxt = diril[i].split(sign)[1].split(".")[0]
                    radif = int(ratxt) - int(RAtxt)
                    decdif = int(dectxt) - int(DECtxt)
                    lst += [np.sqrt(radif**2 + decdif**2)]
                nbrick = diril[np.array(lst).argmin()].split("-")[1].split(".")[0]     #true nearest brick in DB
                signt = 1
                if nbrick[4]=="m":
                    signt = -1
                signo = 1
                if sign=="m":
                    signo = -1
                if abs((int(RAtxt)/10) - (int(nbrick[:4])/10))<0.25 and abs((signo*(int(DECtxt)/10)) - (signt*(int(nbrick[5:])/10)))<0.25:
                    print("Brick found in "+PATH)
                    cont = 1
                    break
                print("No brick found in "+PATH)
        if cont == 0:
            print("No brick found. Aborting...")
            t.write(ibrick+" - missing\n")
            if ibrick==brick_lst[len(brick_lst)-1]:
                break 
            continue
            #print(PATH_lst[r+90])
        else:
            print("Nearest brick found: "+nbrick)   #so, sometimes the nearest brick could be in other RA directory (the ### ones). This should get the true nearest brick
            RAtxt = nbrick[:4]   #correct RAtxt in case changed the RA directory
       
        try:    #create directories to download data
            os.makedirs(outdir+"/"+nbrick)
            os.makedirs(outdir+"/"+nbrick+"/tractor")
            if image:
                os.makedirs(outdir+"/"+nbrick+"/coadd")
        except:
            pass

        # Define progress callback that prints the current percentage completed for the file
        #def progress(filename, size, sent):
        #    sys.stdout.write("%s's progress: %.2f%%   \r" % (filename, float(sent)/float(size)*100) )
        # SCPCLient takes a paramiko transport and progress callback as its arguments.
        scp = SCPClient(ssh.get_transport())
        #download data
        for i in ["tractor","coadd"]:
            try: 
                if i=="tractor":
                    for j in ["tractor.fits","brick.sha256sum"]:
                        print("Downloading "+j+"...")    
                        if "dr9" in PATH and j=="brick.sha256sum":
                            print("no "+str(j))
                            continue
                        scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+j.split(".")[0]+"-"+nbrick+"."+j.split(".")[1], outdir+"/"+nbrick+"/"+i+"/", recursive=True)
                        if j=="tractor.fits":
                            t.write(nbrick+" - "+PATH+"\n")
                if image and i=="coadd":
                    print("Downloading coadd...")
                    scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+nbrick, outdir+"/"+nbrick+"/"+i+"/", recursive=True)
            except:
                ssh = SSHClient()   #connect to DB
                ssh.load_system_host_keys()
                MFA=str(input("MFA key: "))
                ssh.connect(hostname=HOST, username=USER, password=PASS+MFA)
                scp = SCPClient(ssh.get_transport())
                if i=="tractor":
                    for j in ["tractor.fits","brick.sha256sum"]:
                        print("Downloading "+j+"...")
                        try:
                            scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+j.split(".")[0]+"-"+nbrick+"."+j.split(".")[1], outdir+"/"+nbrick+"/"+i+"/", recursive=True)
                            if j=="tractor.fits":
                                t.write(nbrick+" - "+PATH+"\n")
                        except:
                            print("no "+j)
                            continue
                if image and i=="coadd":
                    print("Downloading coadd...")
                    try:
                        scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+nbrick, outdir+"/"+nbrick+"/"+i+"/", recursive=True)
                    except:
                        print("no coadd")
                        continue
    #download individual coadds stated in coadd_lst
    #        for j in coadds_lst:
    #            for k in ["","-g","-i","-r","-z","-W1","-W2","-W3","-W4"]:
    #                for r in ["fits.fz","fits","jpg"]:
    #                    try:
    #                        scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+brick+"/legacysurvey-"+brick+"-"+j+k+"."+r, brick+"/"+i+"/", recursive=True)         
    #                        print("Downloading "+j+k+"."+r)
    #                    except:
    #                        print(j+k+"."+r+": No such file or directory")
    #                        continue
    t.close()
    del brick_lst
