import os
import sys
import numpy as np
from paramiko import SSHClient  #pip install paramiko
from scp import SCPClient   #pip install scp

######################################
USER="user"
HOST="cori.nersc.gov"
MFA=str(input("MFA key: "))
PASS="irispassword"+MFA
RA=172.977   #at least one decimal    #172.977,-19.928 for Abell1300 :D
DEC=-19.928  #at least one decimal 
PATH="/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/"
#coadds_lst=["blobmodel","ccds","chi2","depth","galdepth","image","invvar","maskbits","model","nexp","psfsize","resid","wise","wisemodel","wiseresid"]
######################################

#convert to brick string
RAtxt = str(int(np.around(RA*10,0)))
if len(RAtxt)<=3:
    lst=''
    for i in range(4 - len(RAtxt)):
        lst+="0"
    RAtxt=lst+RAtxt
DECtxt=str(int(np.around(DEC*10,0)))
sign=DECtxt.split("-")
if len(sign)==2:
    sign="m"
else:
    sign="p"
if sign=="m" and len(DECtxt)<=3:
    lst=''
    for i in range(3 - len(DECtxt.split("-")[1])):
        lst+="0"
    DECtxt=lst+DECtxt.split("-")[1]
elif sign=="m" and len(DECtxt)>3:
    DECtxt=DECtxt.split("-")[1]
if sign=="p" and len(DECtxt)<=2:
    lst=''
    for i in range(3 - len(DECtxt)):
        lst+="0"
    DECtxt=lst+DECtxt
print("Input brick set to: "+RAtxt+sign+DECtxt)

ssh = SSHClient()   #connect to DB
ssh.load_system_host_keys()
ssh.connect(hostname=HOST, username=USER, password=PASS)

#search nearest brick
stdin, stdout, stderr = ssh.exec_command("ls "+PATH+"tractor/"+RAtxt[:-1])   #make ls in the RA directory
diri = str(stdout.read()).split("\\n")
diri[0] = diri[0][2:]    #first element has a "b'", remove it 

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
lst = []      #take ra and dec from bricks in diril and take the nearest to RAtxt,DECtxt
for i in range(len(diril)):
    ratxt = diril[i].split("-")[1].split(sign)[0]
    dectxt = diril[i].split(sign)[1].split(".")[0]
    radif = int(ratxt) - int(RAtxt)
    decdif = int(dectxt) - int(DECtxt)
    lst += [np.sqrt(radif**2 + decdif**2)]
brick = diril[np.array(lst).argmin()].split("-")[1].split(".")[0]     #true nearest brick in DB
print("Nearest brick found: "+brick)   #so, sometimes the nearest brick could be in other RA directory (the ### ones). This should get the true nearest brick
RAtxt = brick[:4]   #correct RAtxt in case changed the RA directory

try:    #create directories to download data
    os.makedirs(brick)
    os.makedirs(brick+"/tractor")
    os.makedirs(brick+"/coadd")
except:
    pass

# Define progress callback that prints the current percentage completed for the file
#def progress(filename, size, sent):
#    sys.stdout.write("%s's progress: %.2f%%   \r" % (filename, float(sent)/float(size)*100) )
# SCPCLient takes a paramiko transport and progress callback as its arguments.
scp = SCPClient(ssh.get_transport())
#download data
for i in ["tractor","coadd"]:
    if i=="tractor":
        for j in ["tractor.fits","brick.sha256sum"]:
            print("Downloading "+j+"...")    
            scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+j.split(".")[0]+"-"+brick+"."+j.split(".")[1], brick+"/"+i+"/", recursive=True)
    if i=="coadd":
        print("Downloading coadd...")
        scp.get(PATH+i+"/"+RAtxt[:-1]+"/"+brick, brick+"/"+i+"/", recursive=True)
        
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

