import numpy as np
from mpl_toolkits import mplot3d
from scipy import *
from pylab import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy.matlib as matlib
from shutil import copy
from os import mkdir
import shutil
import os
import logging
from replace_data import replace_data

import warnings
warnings.filterwarnings("ignore")

logging.basicConfig(level=logging.DEBUG)

_log = logging.getLogger(__name__)



def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  dx = dx0+np.arange(1.,n+1.,1.)*a
  return dx

H = 200.
h0 = 140.
om = 2.*np.pi/12.42/3600.
#N0=5.0e-3
u0 = 0.08
f0 = 1.0e-4
Ah = 2.0e-2
Kh = 2.0e-2
Cd = 1.0e-3
N0 = 1.2e-2
is2D=0
if is2D==0:
    b = 500
    a2 = 40**2*200000/b # a2=40**2*500/b whenever b changes a remains. # another set changing a, a2=40**2*1000/b 
    runname = 'TideU0%02dsinN0%03dH%03dho%03dGauRidNewParHeada%06db%03dAh%04dKh%04d_freeslipBotSide_Cdqdt%03d/' % (100*u0,1000*N0,H,h0,a2,b,10000*Ah,10000*Kh,1000*Cd)
else:
    runname = 'TideU0%02dsinN0%03dH%03dho%03dGauRid2DAh%04dKh%04d_freeslipBotSide_Cdqdt%03d/' % (100*u0,1000*N0,H,h0,10000*Ah,10000*Kh,1000*Cd)


comments = ''
outdir0= '../results/' + runname   ########### change everytime compile !!!!!!!!!!!!

indir =outdir0+'/indata/'


# reset f0 in data
shutil.copy('data', 'dataF')
#replace_data('dataF', 'f0', '%1.3e'%f0)


#### Set up the output directory
backupmodel=1
if backupmodel:
    try:
        mkdir(outdir0)
    except:
        import datetime
        import time
        ts = time.time()
        st=datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d%H%M%S')
        shutil.move(outdir0[:-1],outdir0[:-1]+'.bak'+st)
        mkdir(outdir0)
        
        _log.info(outdir0+' Exists')
    outdir=outdir0
    try:
        mkdir(outdir)
    except:
        _log.info(outdir+' Exists')
    outdir=outdir+'input/'
    try:
        mkdir(outdir)
    except:
        _log.info(outdir+' Exists')
    try:
        mkdir(outdir+'/figs/')
    except:
        pass

    copy('./gendataU0_Gaustopo_constN0.py',outdir)
else:
      outdir=outdir+'input/'

outdir_p = '/home/jxchang/projects/def-jklymak/jxchang/HighRes2/input/'
print('outdir_p:'+str(outdir_p))
copy('./gendataU0_Gaustopo_constN0.py',outdir_p)

## Copy some other files
_log.info( "Copying files")

try:
    shutil.rmtree(outdir+'/../code/')
except:
    _log.info("code is not there anyhow")
shutil.copytree('../code', outdir+'/../code/')
shutil.copytree('../python', outdir+'/../python/')

try:
    shutil.rmtree(outdir+'/../build/')
except:
    _log.info("build is not there anyhow")
_log.info(outdir+'/../build/')
mkdir(outdir+'/../build/')

# copy any data that is in the local indata
shutil.copytree('../indata/', outdir+'/../indata/')

shutil.copy('../build/mitgcmuv', outdir+'/../build/mitgcmuv')
#shutil.copy('../build/mitgcmuvU%02d'%u0, outdir+'/../build/mitgcmuv%02d'%u0)
shutil.copy('../build/Makefile', outdir+'/../build/Makefile')
shutil.copy('dataF', outdir+'/data')
shutil.copy('dataF', outdir_p+'/data')
shutil.copy('eedata', outdir)
shutil.copy('eedata', outdir_p)
shutil.copy('data.kl10', outdir)
shutil.copy('data.kl10', outdir_p)
#shutil.copy('data.btforcing', outdir)
try:
    shutil.copy('data.kpp', outdir)
    shutil.copy('data.kpp', outdir_p)
except:
    pass
#shutil.copy('data.rbcs', outdir)
try:
    shutil.copy('data.obcs', outdir)
    shutil.copy('data.obcs', outdir_p)
except:
    pass
try:
    shutil.copy('data.diagnostics', outdir)
    shutil.copy('data.diagnostics', outdir_p)
except:
    pass
try:
    shutil.copy('data.pkg', outdir+'/data.pkg')
    shutil.copy('data.pkg', outdir_p+'/data.pkg')
except:
    pass
try:
    shutil.copy('data.rbcs', outdir+'/data.rbcs')
    shutil.copy('data.rbcs', outdir_p+'/data.rbcs')
except:
    pass

_log.info("Done copying files")



# These must match ../code/SIZE.h
ny = 2*60
nx = 24*62
nz = 40

_log.info('nx %d ny %d', nx, ny)


############## Make the grids #############

# y direction:
yt = 3e3
# x direction
xt = 120e3

nmid = 1200
dx0 = 25.
nleft = int((nx-nmid)/2)
print(nleft)
nright = int((nx-nmid)/2)
print(nright)
xleft = xt/2-dx0*nmid/2
xright = xt/2-dx0*nmid/2
print(xleft)
dx = np.zeros(nx)
ai =  np.zeros(nleft)
for i in range(1,nleft+1):
        ai[i-1]=dx0*math.pow(1.0273943318816,i)
print(ai)
print(sum(ai))
dxleft=np.flipud(ai)
dxright=ai
#dxleft = np.flipud(lininc(nleft,xt/2,dx0))
#dxright = lininc(nright,xt/2,dx0)

dx[0:nleft]=dxleft
dx[(nleft):(nleft+nmid)]=dx0
dx[(nleft+nmid):]=dxright
xg=np.cumsum(np.insert(dx,0,0))
xg = xg-xg[int(np.floor(nx/2))]

xc=0.5*(xg[:-1]+xg[1:])
_log.info('XCoffset=%1.4f'%xc[0])
print(sum(dx[0:nleft]))
print(sum(dx[(nleft):(nleft+nmid)]))
print(sum(dx[(nleft+nmid):]))
print(dx[-10:])

# save dx 
with open(indir+"/delXvar.bin", "wb") as f:
	dx.tofile(f)
f.close()

dy = zeros(ny)+yt/ny
yg = np.cumsum(np.insert(dy,0,0))
yg = yg-yg[int(np.floor(ny/2))]
yc = 0.5*(yg[:-1]+yg[1:])
_log.info('dx %f dy %f', dx[0], dy)

# save dy
with open(indir+"/delYvar.bin", "wb") as f:
        dy.tofile(f)
f.close()

xx, yy = np.meshgrid(xc, yc)

# plot
if 1:
    plot(xg[:-1]/1000.,dx)
    xlim([-40,40])
    ylim([0,1000])
    dxmin=min(dx)
    title('dx_{min}=%2.3fm' %dxmin)
    savefig(outdir+'/figs/dx.png')
if 1:
    clf()
    plot(yg[:-1],dy)
    title('dy')
    savefig(outdir+'figs/dy.png')

# Ocean domain extends from (xo,yo) to (xeast,ynorth)
# (i.e. the ocean spans nx-2, ny-2 grid cells)
# out-of-box-config: xo=yo=0, dx=dy=20 km, ocean extent (0,0)-(1200,1200) km
# model domain includes a land cell surrounding the ocean domain
# The full model domain cell centers are located at:
#    XC[0,:] = -10, +10, ..., +1210 (km)
#    YC[:,0] = -10, +10, ..., +1210 (km)
# and full model domain cell corners are located at:
#    XG[0,:] = -20, 0, ..., 1200 [, 1220] (km)
#    YG[:,0] = -20, 0, ..., 1200 [, 1220] (km)
# where the last value in brackets is not included in the MITgcm grid variable
# and reflects the eastern and northern edge of the model domain respectively.
# See section 2.11.4 of the MITgcm users manual.

# topo
#    BASE   1500   2000   3500    4500   5000   7000   9000
#a^2=7500, 28000, 50000, 150000, 250000,400000,600000,1000000
toposill = h0*np.exp(-(xx/700)**2)

H = 200.
if is2D == 0:
    toposill=-H+toposill
    toposill[toposill<-H]=-H
    print(np.shape(toposill))

    topowallN=10*(b-(xx**2/a2)+(yy-1487))
    topowallN[topowallN<0]=-H
    topowallN[topowallN>0]=0
    print(np.shape(topowallN))

    topowallS=10*(b-(xx**2/a2)-(yy+1487))
    topowallS[topowallS<0]=-H
    topowallS[topowallS>0]=0
    print(np.shape(topowallS))

    topo=np.zeros(np.shape(xx))
    for i in range(len(xc)):
        for j in range(len(yc)):
            topo[j,i]=max(toposill[j,i],topowallN[j,i],topowallS[j,i])
    
    topo[0,:]=0.
    topo[-1,:]=0.

else:
    toposill=-H+toposill
    toposill[toposill<-H]=-H
    topo=np.zeros(np.shape(xx))
    topo=toposill
    topo[0,:]=0.
    topo[-1,:]=0.


print('size of rep topo' + repr(shape(topo)))

# Walls (surrounding domain); generate bathymetry file
#h[:, [0,-1]] = 0   # set ocean depth to zero at east and west walls
#h[[0,-1], :] = 0   # set ocean depth to zero at south and north walls

print('check wall if 0')
print(topo[0,500:800])
print(topo[-1,500:800])

# plot
if 1:
    clf()
    subplot(2,1,1)
    plot(xx/1.e3,topo)
    # xlim([-20.,20.])
    subplot(2,1,2)
    plot(yy/1.e3,topo)
    savefig(outdir+'/figs/topo.png')

if 1:
    clf()
    ax=axes(projection='3d')
    ax.plot_surface(xx, yy, topo,cmap='viridis', edgecolor='none')
    ax.set_title('Surface topo')
    savefig(outdir+'/figs/topo3d.png')


with open(indir+"/topo.bin", "wb") as f:
        topo.tofile(f)
f.close()
# dz:
# dz is from the surface down (right?).  Its saved as positive.

dz=zeros(nz)+H/nz

with open(indir+"/delZvar.bin", "wb") as f:
        dz.tofile(f)
f.close()
print(dz)
print(shape(dz))
z=cumsum(dz)

# temperature profile...
from scipy import integrate
g=9.8
alpha= 2e-4
rhoNil=999.8

T0=35+cumsum(N0**2/g/alpha*(-dz))

with open(indir+"/TRef.bin", "wb") as f:
        T0.tofile(f)
f.close()
print(T0)
print('shpe of T0' + str(np.shape(T0)) )

# save T0 over whole domain
TT0 = np.tile(T0,[nx,ny,1]).T
with open(indir+"/T0.bin", "wb") as f:
        TT0.tofile(f)
print(TT0)
print('shpe of TT0' + str(np.shape(TT0)) )

# plot:
if 1:
    clf()
    plot(T0,z,'-.')
    savefig(outdir+'/figs/TO.png')

# Forcing for boundaries
dt=3720.
time = arange(0,12.*3720.,dt)  # JODY: 12* 3720s (dt) = 12.4hr
print(time/3600./12.4)
om = 2*pi/12.40/3600;
uw = u0*np.sin(om*time)
ue = u0*np.sin(om*time)
# plot:
if 1:
    clf()
    plot(time/3600./12.4,ue,label='Ue')
    plot(time/3600/12.4,uw,label='Uw')
    legend()
    xlabel('t/T')
    ylabel('Vel')
    title('%d' % time[-1])
    savefig(outdir+'/figs/Vels.png')

# try time,nz,ny...

uen=zeros((shape(time)[0],nz,ny))
for j in range(0,ny):
  for i in range(0,nz):
    uen[:,i,j]=ue
#print(uen)

uwn=zeros((shape(time)[0],nz,ny))
print(shape(uwn))
for j in range(0,ny):
  for i in range(0,nz):
    uwn[:,i,j]=uw
#print(uwn)

with open(indir+"/Ue.bin","wb") as f:
  uen.tofile(f)

with open(indir+"/Uw.bin", "wb") as f:
  uwn.tofile(f)

t=zeros((shape(time)[0],nz,ny))
for j in range(0,ny):
        for i in range(0,nz):
                for k in range(0,shape(time)[0]):
                        t[k,i,j]=T0[i]
print(shape(t))
with open(indir+"/Te.bin", "wb") as f:
        t.tofile(f)
f.close()
with open(indir+"/Tw.bin", "wb") as f:
        t.tofile(f)
f.close()


_log.info('Writing info to README')

## Copy some other files

############ Save to README
with open('README','r') as f:
    data=f.read()
with open('README','w') as f:
    import datetime
    import time
    ts = time.time()
    st=datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    f.write( st+'\n')
    f.write( outdir+'\n')
    f.write(comments+'\n\n')
    f.write(data)


_log.info('All Done!')

_log.info('Archiving to home directory')

try:
    shutil.rmtree('../archive/'+runname)
except:
    pass

shutil.copytree(outdir0+'/input/', '../archive/'+runname+'/input')
shutil.copytree(outdir0+'/python/', '../archive/'+runname+'/python')
shutil.copytree(outdir0+'/code', '../archive/'+runname+'/code')

try:
    shutil.copytree(outdir0+'/input/', outdir_p+'../archive/'+runname+'/input')
except:
    pass

try:
    shutil.copytree(outdir0+'/python/', outdir_p+'../archive/'+runname+'/python')
except:
    pass

try:
    shutil.copytree(outdir0+'/code', outdir_p+'../archive/'+runname+'/code')
except:
    pass

exit()



