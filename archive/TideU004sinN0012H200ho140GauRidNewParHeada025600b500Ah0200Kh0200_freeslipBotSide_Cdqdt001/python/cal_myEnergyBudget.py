from xmitgcm import open_mdsdataset
import xgcm
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import pandas as pd
import numpy as np
import math
from scipy import integrate
import xarray as xr
import string

currentDirectory = os.getcwd()
data_dir = currentDirectory[:-7] + '/input/'
print(data_dir)

ds1 = open_mdsdataset(data_dir, geometry='cartesian', endian='<',prefix=['energyvars','statevars','statevars2d'])
ds2 = open_mdsdataset(data_dir, geometry='cartesian', endian='<',prefix=['energymvars'])
t = 0
grid = xgcm.Grid(ds1, periodic=False)
print(grid)

f0 = 1.e-4
g = 9.8
rhoNil=999.8
rhoConst=rhoNil

om=2*np.pi/12.4/3600
alpha = 2e-4
beta = 0e-4
nz = 200
#dz=H/nz 

tR_fname="../indata/TRef.bin"
tRef = np.fromfile(tR_fname)
refSalt=35.
refTemp=tRef[0]
print('refTemp='+ str(refTemp))
print('tRef='+ str(tRef))

rho2=rhoNil*(1-(alpha*(tRef-refTemp)))
rhoS=np.roll(rho2,1)
N2=g/rhoNil*(rho2-rhoS)/ds1['drF'].values
#print('N2='+str(N2))

N2[0]=g/rhoNil*(rhoS[0]-rho2[0])/ds1['drF'][0]
#print(N2)

xmin = 34000
xmax = 50000
numcolt=21
numcolv=21
                                        
ttlen=len(ds1.time)
print('the length of time:' + str(ttlen) )

time1=ds1.coords['time']
time2=ds2.coords['time']
xc=ds1.coords['XC']
xg=ds1.coords['XG']
yc=ds1.coords['YC']
yg=ds1.coords['YG']
z=ds1.coords['Z']

PS = 0.5*(ds1['PHIHYD'].roll(YC=1).values+ds1['PHIHYD'].values) 
PW = 0.5*(ds1['PHIHYD'].roll(XC=1).values+ds1['PHIHYD'].values)

#print(np.shape(PS))
#print('PHIHYD=' + str(ds1['PHIHYD'].isel(YC=0,time=100,Z=2).values))
#print('PS='+ str(PS[100,2,0,:]))
print('PW='+ str(PW[100,:,0,1]))
#print(ds1['PHIHYD'].roll(XC=1)+ds1['PHIHYD'])
#print(ds1['PHIHYD'].roll(XC=1).values[100,2,0,:])
#print(ds1['PHIHYD'].values[100,2,0,:])

#when use interp-> can't with .values???  can't work
#US=xr.DataArray(ds1['UVEL'].interp(XG=xc,YC=yg,kwargs={'fill_value':'extrapolate'}).data,coords=[time1,z,yg,xc],dims=['time','Z','YG','XC'])
#print(US.values)
#VW=(ds1['VVEL'].interp(XC=xg,YG=yc, kwargs={'fill_value': 'None'})).values()

# this US.VW.UC.VC only correct in 2d cases
US=0.5*(ds1['UVEL'].roll(XG=-1).values+ds1['UVEL'].values)
VW=0.5*(ds1['VVEL'].roll(XC=1).values+ds1['VVEL'].values) #VW(when X=0, incorrec
UC=US
VC=ds1['VVEL'].values


#print('U:' + str(ds1['UVEL'].isel(YC=0,time=100,Z=2).values))
#print('US:' +str(US[100,2,0,:]))
#print('V:' +str(ds1['VVEL'].isel(YG=0,time=100,Z=2).values))
#print('VW:' + str(VW[100,2,0,:]))

#depth mean pressure
P0S=np.sum(PS*(ds1['drF']*ds1['hFacS']).values,axis=1)
P0W=np.sum(PW*(ds1['drF']*ds1['hFacW']).values,axis=1)
#print('P0W=' +str(P0W[100,0,:]))
#print(np.shape(P0W))

#depth mean velocity
U0W=((ds1['UVEL']*ds1['drF']*ds1['hFacW']).sum('Z')).values
U0S= np.sum(US*(ds1['drF']*ds1['hFacS']).values,axis=1)
U0C= np.sum(UC*(ds1['drF']*ds1['hFacS']).values,axis=1)
V0W= np.sum(VW*(ds1['drF']*ds1['hFacW']).values,axis=1)
V0S=((ds1['VVEL']*ds1['drF']*ds1['hFacS']).sum('Z')).values
V0C= np.sum(UC*(ds1['drF']*ds1['hFacS']).values,axis=1)
#print(np.shape(U0W))
#print(np.shape(U0S))

ZW=((ds1['drF']*ds1['hFacW']).sum('Z')).values
ZS=((ds1['drF']*ds1['hFacS']).sum('Z')).values
ZC=((ds1['drF']*ds1['hFacC']).sum('Z')).values
#print(np.shape(ZW))
#print('ZW=' +str(ZW[0,:]))
#print(ds1)

for i in range(400):
    if ZW[0,i]!=0:
        P0W[:,0,i]=P0W[:,0,i]/ZW[0,i]
        U0W[:,0,i]=U0W[:,0,i]/ZW[0,i]
        V0W[:,0,i]=V0W[:,0,i]/ZW[0,i]
    else:
        P0W[:,0,i]=0
        U0W[:,0,i]=0
        V0W[:,0,i]=0

#print(P0W[100,0,:])

for i in range(400):
    if ZS[0,i]!=0:
        P0S[:,0,i]=P0S[:,0,i]/ZS[0,i]
        U0S[:,0,i]=U0S[:,0,i]/ZS[0,i]
        V0S[:,0,i]=V0S[:,0,i]/ZS[0,i]
    else:
        P0S[:,0,i]=0
        U0S[:,0,i]=0
        V0S[:,0,i]=0
#print(P0S[100,0,:])
#print(U0W[100,0,:])
for i in range(400):
    if ZC[0,i]!=0:
        U0C[:,0,i]=U0C[:,0,i]/ZC[0,i]
        V0C[:,0,i]=V0C[:,0,i]/ZC[0,i]
    else:
        U0C[:,0,i]=0
        V0C[:,0,i]=0
#SIZE OF DEPTH MEAN VARIABLE:(TIME,Y,X)

ds1['PS']=xr.DataArray(PS,coords=[time1,z,yg,xc],dims=['time','Z','YG','XC'])
ds1['PW']=xr.DataArray(PW,coords=[time1,z,yc,xg],dims=['time','Z','YC','XG'])
ds1['P0S']=xr.DataArray(P0S,coords=[time1,yg,xc],dims=['time','YG','XC'])
ds1['P0W']=xr.DataArray(P0W,coords=[time1,yc,xg],dims=['time','YC','XG'])
ds1['US']=xr.DataArray(US,coords=[time1,z,yg,xc],dims=['time','Z','YG','XC'])
ds1['UC']=xr.DataArray(UC,coords=[time1,z,yc,xc],dims=['time','Z','YC','XC'])
ds1['U0S']=xr.DataArray(U0S,coords=[time1,yg,xc],dims=['time','YG','XC'])
ds1['U0W']=xr.DataArray(U0W,coords=[time1,yc,xg],dims=['time','YC','XG'])
ds1['U0C']=xr.DataArray(U0C,coords=[time1,yc,xc],dims=['time','YC','XC'])
ds1['VW']=xr.DataArray(VW,coords=[time1,z,yc,xg],dims=['time','Z','YC','XG'])
ds1['VC']=xr.DataArray(VC,coords=[time1,z,yc,xc],dims=['time','Z','YC','XC'])
ds1['V0S']=xr.DataArray(V0S,coords=[time1,yg,xc],dims=['time','YG','XC'])
ds1['V0W']=xr.DataArray(V0W,coords=[time1,yc,xg],dims=['time','YC','XG'])
ds1['V0C']=xr.DataArray(V0C,coords=[time1,yc,xc],dims=['time','YC','XC'])
ds1['N2']=xr.DataArray(N2,coords=[z],dims=['Z'])
ds1['tRef']=xr.DataArray(tRef,coords=[z],dims=['Z'])

ds1['upW']=(ds1['UVEL']-ds1['U0W'])*ds1['maskW']
ds1['upS']=(ds1['US']-ds1['U0S'])*ds1['maskS']
ds1['upC']=(ds1['UC']-ds1['U0C'])*ds1['maskC']
ds1['vpW']=(ds1['VW']-ds1['V0W'])*ds1['maskW']
ds1['vpS']=(ds1['VVEL']-ds1['V0S'])*ds1['maskS']
ds1['vpC']=(ds1['VC']-ds1['V0C'])*ds1['maskC']

t=100
if 0:
    plt.clf()
    f, ax = plt.subplots(1, 3, figsize=(11,9) , sharey=True)
    ds1['upW'].isel(time=t,YC=0).plot(ax=ax[0],cmap='RdBu_r',vmax=1,vmin=-1,cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('upW')
    ds1['upS'].isel(time=t,YG=0).plot(ax=ax[1],cmap='RdBu_r',vmax=1,vmin=-1,cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('upS')
    ds1['upC'].isel(time=t,YC=0).plot(ax=ax[2],cmap='RdBu_r',vmax=1,vmin=-1,cbar_kwargs={"label": "", "aspect": 40})
    ax[2].set_title('upC')
    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
            size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myup.png')

if 0:
    plt.clf()
    f, ax = plt.subplots(1, 3, figsize=(11,9) , sharey=True)
    ds1['vpW'].isel(time=t,YC=0).plot(ax=ax[0],cmap='RdBu_r',vmax=1,vmin=-1,cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('vpW')
    ds1['vpS'].isel(time=t,YG=0).plot(ax=ax[1],cmap='RdBu_r',vmax=1,vmin=-1,cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('vpS')
    ds1['vpC'].isel(time=t,YC=0).plot(ax=ax[2],cmap='RdBu_r',vmax=1,vmin=-1,cbar_kwargs={"label": "", "aspect": 40})
    ax[2].set_title('vpC')
    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
            size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myvp.png')


#follow Jody's Diagnostic, let's pretend wVEL wC=wW=wS
ds1['wW']=xr.DataArray(0.5*(ds1['WVEL'].roll(Zl=1).values+ds1['WVEL'].values),coords=[time1,z,yc,xg],dims=['time','Z','YC','XG'])
ds1['wS']=xr.DataArray(0.5*(ds1['WVEL'].roll(Zl=1).values+ds1['WVEL'].values),coords=[time1,z,yg,xc],dims=['time','Z','YG','XC'])
ds1['wC']=xr.DataArray(0.5*(ds1['WVEL'].roll(Zl=1).values+ds1['WVEL'].values),coords=[time1,z,yc,xc],dims=['time','Z','YC','XC'])
print(ds1['wC'])

#presure work
ds1['uPbc']=(ds1['PW']-ds1['P0W'])*ds1['upW']*ds1['drF']*ds1['hFacW']*ds1['maskW']
ds1['vPbc']=(ds1['PS']-ds1['P0S'])*ds1['vpS']*ds1['drF']*ds1['hFacS']*ds1['maskS']
#print(ds1['uPbc'])
#print(ds1['vPbc'])

#kinetic energy
ds1['kEpW'] = 0.5*((ds1['upW']*ds1['upW']+ds1['vpW']*ds1['vpW'] + ds1['wW']*ds1['wW'])*ds1['maskW']) #(6.9)
ds1['kEpS'] = 0.5*((ds1['upS']*ds1['upS']+ds1['vpS']*ds1['vpS'] + ds1['wS']*ds1['wS'])*ds1['maskS'])
ds1['kEpC'] = 0.5*((ds1['upC']*ds1['upC']+ds1['vpC']*ds1['vpC'] + ds1['wC']*ds1['wC'])*ds1['maskC'])

#potential energy
dRho=rhoNil-rhoConst
ds1['rhoC']=(rhoNil*(-alpha*(ds1['THETA']-ds1['tRef'])+beta*(ds1['SALT']-refSalt)))*ds1['maskC']+dRho
#print(ds1['THETA'][100,:,0,198].values)
#print(ds1['rhoC'][100,:,0,198].values)
ds1['Ep']=g*g*ds1['rhoC']*ds1['rhoC']/rhoNil/rhoNil/ds1['N2']
myEbc=ds1['kEpC']+ds1['Ep']
#print(kEpC[100,:,0,198].values)
#print(Ep[100,:,0,198].values)
#print(myEbc[100,:,0,198].values)

if 0:
    plt.clf()
    f, ax = plt.subplots(1, 3, figsize=(11,9) , sharey=True)
    #plt.figure(figsize=(15,6))
    
    ds1['kEpW'].sum('Z').plot(ax=ax[0],cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('kEpW')
    ds1['kEpC'].sum('Z').plot(ax=ax[1],cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('kEpC')
    ds1['kEpS'].sum('Z').plot(ax=ax[2],cbar_kwargs={"label": "", "aspect": 40})
    ax[2].set_title('kEpS')

    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
                size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myKE.png')
    #plt.show()

if 0:
    plt.clf()
    f, ax = plt.subplots(1, 3, figsize=(11,9) , sharey=True)
    
    ds1['kEpC'].sum('Z').plot(ax=ax[0],cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('kEpC')
    ds1['Ep'].sum('Z').plot(ax=ax[1],cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('Ep')
    myEbc.sum('Z').plot(ax=ax[2],cbar_kwargs={"label": "", "aspect": 40})
    ax[2].set_title('Ebc')

    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
                size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myKE&PE&Ebc.png')
    #plt.show()


if 0:
    plt.clf()
    f, ax = plt.subplots(1, 2, figsize=(7.3,9) , sharey=True)
    
    ds1['uPbc'].sum('Z').plot(ax=ax[0],cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('uPbc')
    ds1['vPbc'].sum('Z').plot(ax=ax[1],cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('vPbc')

    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
                size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myuPbc&vPbc.png')
    #plt.show()


ds1['EpW']=xr.DataArray(ds1['Ep'].values,coords=[time1,z,yc,xg],dims=['time','Z','YC','XG'])
ds1['EpS']=xr.DataArray(ds1['Ep'].values,coords=[time1,z,yg,xc],dims=['time','Z','YG','XC'])

ds1['hkEpW']=(ds1['U0W']*ds1['upW']+ds1['V0W']*ds1['vpW'])*ds1['maskW']
ds1['hkEpS']=(ds1['U0S']*ds1['upS']+ds1['V0S']*ds1['vpS'])*ds1['maskS']

ds1['uEbc']=ds1['UVEL']*(ds1['kEpW']+ds1['hkEpW']+ds1['EpW'])*ds1['drF']*ds1['hFacW']*ds1['maskW']
ds1['vEbc']=ds1['VVEL']*(ds1['kEpS']+ds1['hkEpS']+ds1['EpS'])*ds1['drF']*ds1['hFacS']*ds1['maskS']


# calculate Conversion term no use differentiate -> use finite difference make coords fromXG to XC!!!!! W1.W2.ah0
W1=grid.diff(ZW*ds1['U0W'],'X', boundary='extrapolate')#+(ZS*ds1['U0S']).differentiate('YG')
W2=grid.diff(ds1['U0W'],'X', boundary='extrapolate')#+ds1['V0S'].differentiate('YG')
ds1['W']=-W1-W2*ZC
print(ds1['W'].dims)

upupW=ds1['upW']*ds1['upW']*ds1['drF']*ds1['hFacW']*ds1['maskW']
upvpS=ds1['upS']*ds1['vpS']*ds1['drF']*ds1['hFacS']*ds1['maskS']
upvpW=ds1['upW']*ds1['vpW']*ds1['drF']*ds1['hFacW']*ds1['maskW']
vpvpS=ds1['vpS']*ds1['vpS']*ds1['drF']*ds1['hFacS']*ds1['maskS']

conv1=ds1['rhoC']*g*ds1['W']*ds1['drF']*ds1['hFacC']/rhoNil
ah0=ds1['U0C']*grid.diff(upupW,'X', boundary='extrapolate')+ds1['V0C']*grid.diff(upvpW,'X', boundary='extrapolate')#+U0C*(upvS.differentiate('YG'))+V0C*(upvpS.differentiate('YG'))
print(conv1.dims)
print(ah0.dims)

ds1['Conv']=conv1+ah0
print(ds1['Conv'].dims)

if 0:
    plt.clf()
    f, ax = plt.subplots(1, 2, figsize=(7.3,9) , sharey=True)
    
    ds1['hkEpW'].sum('Z').plot(ax=ax[0],cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('hkEpW')
    ds1['hkEpS'].sum('Z').plot(ax=ax[1],cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('hkEpS')

    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
                size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myhkEp.png')


if 0:
    plt.clf()
    f, ax = plt.subplots(1, 2, figsize=(7.3,9) , sharey=True)
    
    ds1['uEbc'].sum('Z').plot(ax=ax[0],cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('uEbc')
    ds1['vEbc'].sum('Z').plot(ax=ax[1],cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('vEbc')

    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
                size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myuEbc&vEbc.png')


if 1:
    plt.clf()
    f, ax = plt.subplots(1, 2, figsize=(7.3,9) , sharey=True)
    
    ds2['SDIAG10'].isel(YG=0).plot(ax=ax[0],cbar_kwargs={"label": "", "aspect": 40})
    ax[0].set_title('His Conv')
    ds1['Conv'].sum('Z').plot(ax=ax[1],cbar_kwargs={"label": "", "aspect": 40})
    ax[1].set_title('My Conv')

    for n, axs in enumerate(ax):
        axs.text(-0.1, 1, string.ascii_lowercase[n], transform=axs.transAxes,
                size=20, weight='bold')
    plt.tight_layout()
    plt.savefig('./figs/myConv.png')
