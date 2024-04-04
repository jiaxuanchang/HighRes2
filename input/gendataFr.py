import numpy as np
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


logging.basicConfig(level=logging.DEBUG)

_log = logging.getLogger(__name__)



def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  dx = dx0+np.arange(1.,n+1.,1.)*a
  return dx

Fr=0.05
H =200.
h0=140.
om = 2.*np.pi/12.42/3600.
N0=5.0e-3
u0=Fr*N0*h0
f0 = 1.0e-4

runname = 'TideFr%04dN0linH%03dho%03d/' % (10000*Fr,H,h0)
comments = ''
outdir0= '../results/' + runname   ########### change everytime compile !!!!!!!!!!!!

indir =outdir0+'/indata/'


# reset f0 in data
shutil.copy('data', 'dataF')
replace_data('dataF', 'f0', '%1.3e'%f0)


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

    copy('./gendataFr.py',outdir)
else:
      outdir=outdir+'input/'

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
shutil.copy('eedata', outdir)
shutil.copy('data.kl10', outdir)
#shutil.copy('data.btforcing', outdir)
try:
    shutil.copy('data.kpp', outdir)
except:
    pass
#shutil.copy('data.rbcs', outdir)
try:
    shutil.copy('data.obcs', outdir)
except:
    pass
try:
    shutil.copy('data.diagnostics', outdir)
except:
    pass
try:
    shutil.copy('data.pkg', outdir+'/data.pkg')
except:
    pass
try:
    shutil.copy('data.rbcs', outdir+'/data.rbcs')
except:
    pass

_log.info("Done copying files")



# These must match ../code/SIZE.h
ny = 1
nx = 4*100
nz = 40

_log.info('nx %d ny %d', nx, ny)


############## Make the grids #############

# y direction:
dy = 1000.0
# x direction
xt = 80e3

nmid = 200
dx0 = 20.
nleft = int((nx-nmid)/2)
print(nleft)
nright = int((nx-nmid)/2)
dx = np.zeros(nx)
dxleft = np.flipud(lininc(nleft,xt/2,dx0))
dxright = lininc(nright,xt/2,dx0)
dx[0:nleft]=dxleft
dx[(nleft):(nleft+nmid)]=dx0
dx[(nleft+nmid):]=dxright
x=np.cumsum(dx)
x = x-x[int(np.floor(nx/2))]
_log.info('XCoffset=%1.4f'%x[0])

_log.info('dx %f dy %f', dx[0], dy)

# save dx 
with open(indir+"/delXvar.bin", "wb") as f:
	dx.tofile(f)
f.close()

# plot
if 1:
    plot(x/1000.,dx)
    xlim([-40,40])
    ylim([0,1000])
    dxmin=min(dx)
    title('dx_{min}=%2.3fm' %dxmin)
    savefig(outdir+'/figs/dx.png')

# topo
sigma = 2000. # m

topo = 600*exp(-x*x/(sigma**2))-600+h0
#topo = h0*exp(-x*x/(3000**2))
print(shape(topo))
print(topo)
topo[topo<0.]=0.
topo=-H+topo
topo[topo<-H]=-H
print(topo)

# plot
if 1:
    clf()
    plot(x/1.e3,topo)
    # xlim([-20.,20.])
    savefig(outdir+'/figs/topo.png')


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
import gsw
from scipy.io import loadmat
from scipy.interpolate import interp1d

alpha = 2e-4*1000
x = loadmat('density_bk.mat')
rho = x['Sigr']
p = x['grid_p']*100
TT= 35-(rho-1022)/alpha
p=p.flatten()
TT=TT.flatten()
id=np.argwhere(~np.isnan(TT))
fT = interp1d(p[id].flatten(),TT[id].flatten(),fill_value='extrapolate')
T0 = fT(z)

g = 9.8
#T0 = 35+cumsum(N0**2/g/alpha*(-dz))

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
    plot(T0,z)
    savefig(outdir+'/figs/TO.png')

# Forcing for boundaries
dt=3720.
time = arange(0,12.*3720.,dt)  # JODY: 12* 3720s (dt) = 12.4hr
print(time/3600./12.4)
om = 2*pi/12.40/3600;
uw = u0*np.cos(om*time)
ue = u0*np.cos(om*time)
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

exit()



