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


logging.basicConfig(level=logging.DEBUG)

_log = logging.getLogger(__name__)



om = 2.*np.pi/12.42/3600.
#N0=5.0e-3
u0 = 0.24
f0 = 1.0e-4
Ah = 2.0e-2
Kh = 2.0e-2
N0 = 1.2e-2
Cd = 1.0e-3

runname = 'TideU0%02dsinN0%03dAh%04dKh%04d_freeslipBotSide_Cdqdt%03d/' % (100*u0,1000*N0,10000*Ah,10000*Kh,1000*Cd)
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

    copy('./gendataU0constN0.py',outdir)
else:
      outdir=outdir+'input/'

outdir_p = '/home/jxchang/projects/def-jklymak/jxchang/HighRes2/input/'
print('outdir_p:'+str(outdir_p))
copy('./gendataU0constN0.py',outdir_p)

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
shutil.copy('topo.bin', outdir_p+'../indata')
shutil.copy('topo.bin', outdir+'../indata')

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
ny = 2*48
nx = 12*50
nz = 40

_log.info('nx %d ny %d', nx, ny)


############## Make the grids #############

dx=50
dy=50

x=np.arange(-15000,15050,dx)
y=np.arange(-1600,3250,dy)

_log.info('dx %f dy %f', dx, dy)

xx, yy = np.meshgrid(x, y)

dx=np.diff(x)
dy=np.diff(y)
print(dx)

# save dx
with open(indir+"/delXvar.bin", "wb") as f:
        dx.tofile(f)
f.close()

# save dy
with open(indir+"/delYvar.bin", "wb") as f:
        dy.tofile(f)
f.close()

# plot
if 1:
    plot(x[:-1]/1000.,dx)
    #xlim([-20,20])
    #ylim([0,80])
    title('dx')
    savefig(outdir+'/figs/dx.png')
if 1:
    clf()
    plot(y[:-1],dy)
    title('dy')
    savefig(outdir+'figs/dy.png')


# dz:
# dz is from the surface down (right?).  Its saved as positive.

H = 452
dz=zeros(nz)+H/nz

with open(indir+"/delZvar.bin", "wb") as f:
	dz.tofile(f)
f.close()
print(dz)
print(shape(dz))
z=cumsum(dz)

# temperature profile...

alpha = 2e-4
g = 9.8
T0 = 50+cumsum(N0**2/g/alpha*(-dz))
print('T0 profile:' + str(T0))

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
STOP

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



