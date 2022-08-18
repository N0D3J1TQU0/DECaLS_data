import os
import sys
import numpy as np
from paramiko import SSHClient  #pip install paramiko
from scp import SCPClient   #pip install scp

######################################
USER="user"
HOST="cori.nersc.gov"
PASS="irispassword"
RA=[150,210]   #can be float or list with min and max values
DEC=[-30,-20]  #can be float or list with min and max values
PATH_lst=["/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/","/global/cscratch1/sd/comparat/dr10/","/global/cscratch1/sd/mxhf/dr10c-test/","/global/cscratch1/sd/jsnigula/dr10c/","/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/"]
#use brick_lst to search and download specific bricks, brick_lst has priority over RA,DEC (i.e. RA,DEC ignored)
#brick_lst = ['1686p007', '1688p007', '1691p007', '1693p007', '1683p010', '1686p010', '1688p010', '1691p010', '1693p010', '1696p010', '1681p012', '1683p012', '1686p012', '1688p012', '1691p012', '1693p012', '1696p012', '1681p015', '1683p015', '1686p015', '1688p015', '1691p015', '1693p015', '1696p015', '1681p017', '1683p017', '1686p017', '1688p017', '1691p017', '1693p017', '1696p017', '1683p020', '1686p020', '1688p020', '1691p020', '1693p020', '1696p020', '1686p022', '1688p022', '1691p022', '1693p022']
brick_lst = ['1973m040', '1976m040', '1978m040', '1981m040', '1971m037', '1973m037', '1976m037', '1978m037', '1981m037', '1984m037', '1971m035', '1973m035', '1976m035', '1978m035', '1981m035', '1984m035', '1968m032', '1971m032', '1973m032', '1976m032', '1978m032', '1981m032', '1984m032', '1971m030', '1973m030', '1976m030', '1978m030', '1981m030', '1983m030', '1986m030', '1971m027', '1973m027', '1976m027', '1978m027', '1981m027', '1983m027', '1971m025', '1973m025', '1976m025', '1978m025', '1981m025', '1983m025']
#brick_lst = ['2064m125', '2067m125', '2069m125', '2072m125', '2062m122', '2064m122', '2067m122', '2069m122', '2072m122', '2074m122', '2061m120', '2064m120', '2066m120', '2069m120', '2071m120', '2074m120', '2077m120', '2061m117', '2063m117', '2066m117', '2068m117', '2071m117', '2074m117', '2076m117', '2061m115', '2063m115', '2066m115', '2068m115', '2071m115', '2074m115', '2076m115', '2060m112', '2063m112', '2066m112', '2068m112', '2071m112', '2073m112', '2076m112', '2065m110', '2068m110', '2070m110', '2073m110']

image=False    #download coadds?
######################################

#convert to brick string
try:
    len(brick_lst)
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

ssh = SSHClient()   #connect to DB
ssh.load_system_host_keys()
MFA=str(input("MFA key: "))
ssh.connect(hostname=HOST, username=USER, password=PASS+MFA)

if area:
    ratoxt = np.linspace(int(RAtxt_inf),int(RAtxt_sup),int((int(RAtxt_sup)-int(RAtxt_inf))/10)+1)/10
    diril = []
    print("Locating bricks...")
    for i in range(len(ratoxt)):
        ratxt = str(int(ratoxt[i]))
        for PATH in PATH_lst:
            DECtxt_lst.sort()
            stdin, stdout, stderr = ssh.exec_command("ls "+PATH+"tractor/"+ratxt+"/tractor*"+sign+"["+DECtxt_lst[0][1]+"-"+DECtxt_lst[-1][1]+"]*")   #make ls in the RA directory
            diri = str(stdout.read()).split("\\n")
            if len(diri)==0:
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

t = open("bricks_info.txt",'w')
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
        os.makedirs(nbrick)
        os.makedirs(nbrick+"/tractor")
        if image:
            os.makedirs(nbrick+"/coadd")
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
                    scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+j.split(".")[0]+"-"+nbrick+"."+j.split(".")[1], nbrick+"/"+i+"/", recursive=True)
                    if j=="tractor.fits":
                        t.write(nbrick+" - "+PATH+"\n")
            if image and i=="coadd":
                print("Downloading coadd...")
                scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+nbrick, nbrick+"/"+i+"/", recursive=True)
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
                        scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+j.split(".")[0]+"-"+nbrick+"."+j.split(".")[1], nbrick+"/"+i+"/", recursive=True)
                        if j=="tractor.fits":
                            t.write(nbrick+" - "+PATH+"\n")
                    except:
                        print("no "+j)
                        continue
            if image and i=="coadd":
                print("Downloading coadd...")
                try:
                    scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+nbrick, nbrick+"/"+i+"/", recursive=True)
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
