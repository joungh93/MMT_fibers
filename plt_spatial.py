#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 20:15:40 2019

@author: jlee
"""


import time
start_time = time.time()

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pandas as pd
from PIL import Image
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord


# ----- Reading known GC data ----- #

# NH10 GCs
di_NH10 = "/data/jlee/HSCv6/M81/Bell/diagram/NH10/"
df_NH10 = pd.read_pickle(di_NH10+"df_NH10.pkl")

# NW18 GCs
di_NW18 = "/data/jlee/HSCv6/M81/Bell/diagram/LeeNW+18/"
df_NW18 = pd.read_pickle(di_NW18+"df_NW18.pkl")

NWGC = (df_NW18['NH_Obj'] == "none")


# ----- Reading new GC data ----- #

## Bell fields
dir_bel = "/data/jlee/HSCv6/M81/Bell/Phot2/diagram/"
dir_Cat = "/data/jlee/HSCv6/M81/Bell/Phot2/HSC_sephot/Catalogs/"

# # Before visual inspection
# df_bel = pd.read_pickle(dir_Cat+"df_sep_G.pkl")
# for i in np.arange(3):
# 	exec(f"gc{i+1:d}_bel = pd.read_pickle('"+dir_bel+f"src_pr{i+1:d}GC.pkl')[0]")
# 	exec(f"n_gc{i+1:d}_bel = np.sum(gc{i+1:d}_bel)")

# After visual inspection
df_bel = pd.read_pickle(dir_Cat+"df_sep_G.pkl")
for i in np.arange(4):
	dir_vis = dir_bel+f"vis{i+1:d}/"
	dv = pd.read_csv(dir_vis+"vis_results.csv") 
	if (i < 3):
		exec(f"gc{i+1:d}_bel = dv['ID'][dv['vis_flag'] == 1]")
		exec(f"n_gc{i+1:d}_bel = len(gc{i+1:d}_bel)")
	else:
		ysc_bel = dv['ID'][dv['vis_flag'] == 1]
		n_ysc_bel = len(ysc_bel)

## Okamoto fields
dir_oka = "/data2/jlee/HSCv6/M81/Okamoto/F1234/Phot/diagram/"
dir_Cat = "/data2/jlee/HSCv6/M81/Okamoto/F1234/Phot/HSC_sephot/Catalogs/"

# # Before visual inspection
# df_oka = pd.read_pickle(dir_Cat+"df_sep_G.pkl")
# for i in np.arange(3):
# 	exec(f"gc{i+1:d}_oka = pd.read_pickle('"+dir_oka+f"src_pr{i+1:d}GC.pkl')[0]")
# 	exec(f"n_gc{i+1:d}_oka = np.sum(gc{i+1:d}_oka)")

# After visual inspection
df_oka = pd.read_pickle(dir_Cat+"df_sep_G.pkl")
for i in np.arange(4):
	dir_vis = dir_oka+f"vis{i+1:d}/"
	dv = pd.read_csv(dir_vis+"vis_results.csv") 
	if (i < 3):
		exec(f"gc{i+1:d}_oka = dv['ID'][dv['vis_flag'] == 1]")
		exec(f"n_gc{i+1:d}_oka = len(gc{i+1:d}_oka)")
	else:
		ysc_oka = dv['ID'][dv['vis_flag'] == 1]
		n_ysc_oka = len(ysc_oka)


# ----- Figure 1 : Spatial distribution ----- #
bkg_img = "/data/jlee/HSCv6/M81/MMT_2020A/okamoto_1.JPG"
bkg = Image.open(bkg_img)

fig1 = plt.figure(1, figsize=(13,10))
ax1 = fig1.add_subplot(1,1,1)
ax1.set_position([0.13,0.10,0.63,0.82])
ax1.set_xticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
ax1.set_xticklabels([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5], fontsize=22.5)
ax1.set_yticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5])
ax1.set_yticklabels([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5], fontsize=22.5)
ax1.set_xlabel(r'$\Delta$R.A. [deg]', fontsize=22.5)
ax1.set_ylabel(r'$\Delta$Decl. [deg]', fontsize=22.5)
ax1.set_xlim([1.3,-0.8])
ax1.set_ylim([-0.9,1.2])
ax1.tick_params(width=2.0, length=12.0)
plt.minorticks_on()
ax1.tick_params(width=2.0, length=8.0, which='minor')
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(2.0)
# -------------------- #

# # ----- Background image (red RGB map) ----- #
# ax1.imshow(bkg, extent=[1.2,-0.3,-0.9,1.0])

# # ----- Field-of-view ----- #
# x0_fov, y0_fov = [0.2679, 0.6935, 0.0], [0.6166, -0.3314, 0.0]

# for i in np.arange(len(x0_fov)):
# 	p1 = plt.Circle((x0_fov[i], y0_fov[i]), radius=0.85, color='gray',
# 		            linewidth=2.75, linestyle='--', alpha=0.85, fill=False)
# 	ax1.add_artist(p1)

# ----- Centers of M81, M82, and NGC3077 ----- #
ra0_m81, dec0_m81 = 148.88822, 69.06529
c0 = SkyCoord(ra=ra0_m81, dec=dec0_m81, unit=u.deg)

x0_m82, y0_m82 = 0.0291, 0.6140
x0_n3077, y0_n3077 = 0.69358, -0.33137
x0_h9, y0_h9 = 0.17706, -0.01946
x0_bk3n, y0_bk3n = -0.15646, -0.09645

ax1.plot(0.0, 0.0, '*', ms=25.0, color='orange', mec='k', mew=2.0, zorder=10)
ax1.plot(x0_m82, y0_m82, '*', ms=25.0, color='orange', mec='k', mew=2.0, zorder=10)
ax1.plot(x0_n3077, y0_n3077, '*', ms=25.0, color='orange', mec='k', mew=2.0, zorder=10)
ax1.plot(x0_h9, y0_h9, '*', ms=17.5, color='orange', mec='k', mew=1.5, zorder=10)
ax1.plot(x0_bk3n, y0_bk3n, '*', ms=17.5, color='orange', mec='k', mew=1.5, zorder=10)


# ----- Configuration 1 ----- #
cfg1 = pd.read_csv('./match_cfg1.csv')
ra1, dec1 = np.zeros(len(cfg1)), np.zeros(len(cfg1))
for i in np.arange(len(cfg1)):
	rah = float(cfg1['ra_y'].values[i].split(':')[0])
	ram = float(cfg1['ra_y'].values[i].split(':')[1])
	ras = float(cfg1['ra_y'].values[i].split(':')[2])
	ra1[i] = 15.0*(rah + ram/60.0 + ras/3600.0)

	ded = float(cfg1['dec_y'].values[i].split(':')[0])
	dem = float(cfg1['dec_y'].values[i].split(':')[1])
	des = float(cfg1['dec_y'].values[i].split(':')[2])
	dec1[i] = ded + dem/60.0 + des/3600.0
dra1, ddec1 = (ra1-ra0_m81)*np.cos(dec1*np.pi/180.0), (dec1-dec0_m81)

gcc1  = cfg1.loc[cfg1['object'].str.startswith('GC1')]
gc21  = cfg1.loc[cfg1['object'].str.startswith('GC2')]
gc31  = cfg1.loc[cfg1['object'].str.startswith('GC3')]
sc1   = cfg1.loc[cfg1['object'].str.startswith('SC4')]
ysc1  = cfg1.loc[cfg1['object'].str.startswith('ySC')]
kgc1  = cfg1.loc[cfg1['object'].str.startswith('kGC')]
snr1  = cfg1.loc[cfg1['object'].str.startswith('SNR')]
pnh1  = cfg1.loc[cfg1['object'].str.startswith('PNH')]
x1    = cfg1.loc[cfg1['object'].str.startswith('Xra')]
star1 = cfg1.loc[cfg1['object'].str.startswith('Sta')]
gal1  = cfg1.loc[cfg1['object'].str.startswith('Gal')]


# ----- Configuration 2 ----- #
cfg2 = pd.read_csv('./match_cfg2.csv')
ra2, dec2 = np.zeros(len(cfg2)), np.zeros(len(cfg2))
for i in np.arange(len(cfg2)):
	rah = float(cfg2['ra_y'].values[i].split(':')[0])
	ram = float(cfg2['ra_y'].values[i].split(':')[1])
	ras = float(cfg2['ra_y'].values[i].split(':')[2])
	ra2[i] = 15.0*(rah + ram/60.0 + ras/3600.0)

	ded = float(cfg2['dec_y'].values[i].split(':')[0])
	dem = float(cfg2['dec_y'].values[i].split(':')[1])
	des = float(cfg2['dec_y'].values[i].split(':')[2])
	dec2[i] = ded + dem/60.0 + des/3600.0
dra2, ddec2 = (ra2-ra0_m81)*np.cos(dec2*np.pi/180.0), (dec2-dec0_m81)

gcc2  = cfg2.loc[cfg2['object'].str.startswith('GC1')]
gc22  = cfg2.loc[cfg2['object'].str.startswith('GC2')]
gc32  = cfg2.loc[cfg2['object'].str.startswith('GC3')]
sc2   = cfg2.loc[cfg2['object'].str.startswith('SC4')]
ysc2  = cfg2.loc[cfg2['object'].str.startswith('ySC')]
kgc2  = cfg2.loc[cfg2['object'].str.startswith('kGC')]
snr2  = cfg2.loc[cfg2['object'].str.startswith('SNR')]
pnh2  = cfg2.loc[cfg2['object'].str.startswith('PNH')]
x2    = cfg2.loc[cfg2['object'].str.startswith('Xra')]
star2 = cfg2.loc[cfg2['object'].str.startswith('Sta')]
gal2  = cfg2.loc[cfg2['object'].str.startswith('Gal')]


# ----- Configuration 3 ----- #
cfg3 = pd.read_csv('./match_cfg3.csv')
ra3, dec3 = np.zeros(len(cfg3)), np.zeros(len(cfg3))
for i in np.arange(len(cfg3)):
	rah = float(cfg3['ra_y'].values[i].split(':')[0])
	ram = float(cfg3['ra_y'].values[i].split(':')[1])
	ras = float(cfg3['ra_y'].values[i].split(':')[2])
	ra3[i] = 15.0*(rah + ram/60.0 + ras/3600.0)

	ded = float(cfg3['dec_y'].values[i].split(':')[0])
	dem = float(cfg3['dec_y'].values[i].split(':')[1])
	des = float(cfg3['dec_y'].values[i].split(':')[2])
	dec3[i] = ded + dem/60.0 + des/3600.0
dra3, ddec3 = (ra3-ra0_m81)*np.cos(dec3*np.pi/180.0), (dec3-dec0_m81)

gcc3  = cfg3.loc[cfg3['object'].str.startswith('GC1')]
gc23  = cfg3.loc[cfg3['object'].str.startswith('GC2')]
gc33  = cfg3.loc[cfg3['object'].str.startswith('GC3')]
sc3   = cfg3.loc[cfg3['object'].str.startswith('SC4')]
ysc3  = cfg3.loc[cfg3['object'].str.startswith('ySC')]
kgc3  = cfg3.loc[cfg3['object'].str.startswith('kGC')]
snr3  = cfg3.loc[cfg3['object'].str.startswith('SNR')]
pnh3  = cfg3.loc[cfg3['object'].str.startswith('PNH')]
x3    = cfg3.loc[cfg3['object'].str.startswith('Xra')]
star3 = cfg3.loc[cfg3['object'].str.startswith('Sta')]
gal3  = cfg3.loc[cfg3['object'].str.startswith('Gal')]


# ----- MMT FOV ----- #
df_fld = pd.read_csv('target1.fld', header=None, skiprows=2, sep='\t',
	                 usecols=(0,1), names=('ra', 'dec'))
mmt_fld = []
for i in np.arange(len(df_fld)):
	ra1 = Angle(df_fld['ra'].values[i], unit=u.hour).deg
	dec1 = Angle(df_fld['dec'].values[i], unit=u.deg).deg
	c1R = SkyCoord(ra=ra1, dec=dec0_m81, unit=u.deg)
	c1D = SkyCoord(ra=ra0_m81, dec=dec1, unit=u.deg)
	dra = c1R.separation(c0).deg
	ddec = c1D.separation(c0).deg
	if (ra1 < ra0_m81):
		dra *= -1
	if (dec1 < dec0_m81):
		ddec *= -1
	mmt_fld.append([dra, ddec])

# mmt_fld = [[-0.15,0.05], [0.3,0.6], [0.69358, -0.33137]]
# mmt_fld = [[-0.1400,0.0533], [0.2850,0.5983], [0.6800, -0.3247]]
# mmt_fld = [[-0.20,0.0], [0.25,0.60], [0.6800, -0.3247]]

for i in np.arange(np.shape(mmt_fld)[0]):
	p1 = plt.Circle((mmt_fld[i][0], mmt_fld[i][1]), radius=0.51, color='k',
		            linewidth=3.25, linestyle='-', alpha=0.5, fill=False)
	ax1.add_artist(p1)

# ax1.text(-0.60, 0.25, 'F1', fontsize=27.5, fontweight='bold', color='k')
ra_f1, dec_f1 = ra0_m81 + mmt_fld[0][0]/np.cos((dec0_m81+mmt_fld[0][1])*np.pi/180.0), dec0_m81+mmt_fld[0][1]
print("F1 center : (%02d:%02d:%05.2f, +%02d:%02d:%05.2f)" \
	  %(Angle(ra_f1, u.deg).hms[0], Angle(ra_f1, u.deg).hms[1], Angle(ra_f1, u.deg).hms[2],
	  	Angle(dec_f1, u.deg).dms[0], Angle(dec_f1, u.deg).dms[1], Angle(dec_f1, u.deg).dms[2]))

# ax1.text(0.85, 1.0, 'F2', fontsize=27.5, fontweight='bold', color='k')
ra_f2, dec_f2 = ra0_m81 + mmt_fld[1][0]/np.cos((dec0_m81+mmt_fld[1][1])*np.pi/180.0), dec0_m81+mmt_fld[1][1]
print("F2 center : (%02d:%02d:%05.2f, +%02d:%02d:%05.2f)" \
	  %(Angle(ra_f2, u.deg).hms[0], Angle(ra_f2, u.deg).hms[1], Angle(ra_f2, u.deg).hms[2],
	  	Angle(dec_f2, u.deg).dms[0], Angle(dec_f2, u.deg).dms[1], Angle(dec_f2, u.deg).dms[2]))

# ax1.text(1.15, -0.9, 'F3', fontsize=27.5, fontweight='bold', color='k')
ra_f3, dec_f3 = ra0_m81 + mmt_fld[2][0]/np.cos((dec0_m81+mmt_fld[2][1])*np.pi/180.0), dec0_m81+mmt_fld[2][1]
print("F3 center : (%02d:%02d:%05.2f, +%02d:%02d:%05.2f)" \
	  %(Angle(ra_f3, u.deg).hms[0], Angle(ra_f3, u.deg).hms[1], Angle(ra_f3, u.deg).hms[2],
	  	Angle(dec_f3, u.deg).dms[0], Angle(dec_f3, u.deg).dms[1], Angle(dec_f3, u.deg).dms[2]))


# ----- Stars ----- #
ax1.plot(dra1[star1.index.values], ddec1[star1.index.values], 'o', ms=3.0, color='gray', alpha=0.5)
ax1.plot(dra2[star2.index.values], ddec2[star2.index.values], 'o', ms=3.0, color='gray', alpha=0.5)
ax1.plot(dra3[star3.index.values], ddec3[star3.index.values], 'o', ms=3.0, color='gray', alpha=0.5)
sym_star, = ax1.plot(-100.0,-100.0, 'o', ms=3.0, color='gray', alpha=0.5, label='Field stars')

# ----- Galaxies ----- #
ax1.plot(dra1[gal1.index.values], ddec1[gal1.index.values], 'o', ms=3.0, color='dodgerblue', alpha=0.7)
ax1.plot(dra2[gal2.index.values], ddec2[gal2.index.values], 'o', ms=3.0, color='dodgerblue', alpha=0.7)
ax1.plot(dra3[gal3.index.values], ddec3[gal3.index.values], 'o', ms=3.0, color='dodgerblue', alpha=0.7)
sym_gal, = ax1.plot(-100.0,-100.0, 'o', ms=3.0, color='dodgerblue', alpha=0.5, label='Background galaxies')

# ----- PNe + HII regions ----- #
ax1.plot(dra1[pnh1.index.values], ddec1[pnh1.index.values], 'v', ms=7.5, color='green', alpha=0.8)
ax1.plot(dra2[pnh2.index.values], ddec2[pnh2.index.values], 'v', ms=7.5, color='green', alpha=0.8)
ax1.plot(dra3[pnh3.index.values], ddec3[pnh3.index.values], 'v', ms=7.5, color='green', alpha=0.8)
sym_pnh, = ax1.plot(-100.0,-100.0, 'v', ms=7.5, color='green', alpha=0.8, label='PNe/HII regions')

# ----- SNRs ----- #
ax1.plot(dra1[snr1.index.values], ddec1[snr1.index.values], '^', ms=7.5, color='brown', alpha=0.8)
ax1.plot(dra2[snr2.index.values], ddec2[snr2.index.values], '^', ms=7.5, color='brown', alpha=0.8)
ax1.plot(dra3[snr3.index.values], ddec3[snr3.index.values], '^', ms=7.5, color='brown', alpha=0.8)
sym_snr, = ax1.plot(-100.0,-100.0, '^', ms=7.5, color='brown', alpha=0.8, label='SNRs')

# ----- X-ray sources ----- #
ax1.plot(dra1[x1.index.values], ddec1[x1.index.values], 'x', ms=9.0, mec='blueviolet', mew=2.5, alpha=0.9)
ax1.plot(dra1[x2.index.values], ddec1[x2.index.values], 'x', ms=9.0, mec='blueviolet', mew=2.5, alpha=0.9)
ax1.plot(dra1[x3.index.values], ddec1[x3.index.values], 'x', ms=9.0, mec='blueviolet', mew=2.5, alpha=0.9)
sym_x, = ax1.plot(-100.0,-100.0, 'x', ms=9.0, mec='blueviolet', mew=2.5, alpha=0.9, label='X-ray sources')


# ----- Young star clusters ----- #
ax1.plot(dra2[ysc1.index.values], ddec2[ysc1.index.values], 'o', ms=8.0, color='cyan', mec='k', mew=1.0, alpha=0.9)
ax1.plot(dra2[ysc2.index.values], ddec2[ysc2.index.values], 'o', ms=8.0, color='cyan', mec='k', mew=1.0, alpha=0.9)
ax1.plot(dra2[ysc3.index.values], ddec2[ysc3.index.values], 'o', ms=8.0, color='cyan', mec='k', mew=1.0, alpha=0.9)
sym_ysc, = ax1.plot(-100.0,-100.0, 'o', ms=8.0, color='cyan', mec='k', mew=1.0, alpha=0.9, label='YSCs (Lim+13)')


# ---- Known GCs ----- #
ax1.plot(dra1[kgc1.index.values], ddec1[kgc1.index.values], 's', ms=6.0, color='blue', alpha=0.8)
ax1.plot(dra2[kgc2.index.values], ddec2[kgc2.index.values], 's', ms=6.0, color='blue', alpha=0.8)
ax1.plot(dra3[kgc3.index.values], ddec3[kgc3.index.values], 's', ms=6.0, color='blue', alpha=0.8)
sym_kgc, = ax1.plot(-100.0,-100.0, 's', ms=6.0, color='blue', alpha=0.8, label='Known GCs')


# ----- New GCs ----- #
ax1.plot(dra1[gcc1.index.values], ddec1[gcc1.index.values], 'o', ms=8.0, fillstyle='none', mec='red', mew=1.5)
ax1.plot(dra2[gcc2.index.values], ddec2[gcc2.index.values], 'o', ms=8.0, fillstyle='none', mec='red', mew=1.5)
ax1.plot(dra3[gcc3.index.values], ddec3[gcc3.index.values], 'o', ms=8.0, fillstyle='none', mec='red', mew=1.5)
sym_gcc, = ax1.plot(-100.0,-100.0, 'o', ms=8.0, fillstyle='none', mec='red', mew=1.5, label='Rank 1 GCs')

ax1.plot(dra1[gc21.index.values], ddec1[gc21.index.values], 'o', ms=8.0, fillstyle='none', mec='darkorange', mew=1.5)
ax1.plot(dra2[gc22.index.values], ddec2[gc22.index.values], 'o', ms=8.0, fillstyle='none', mec='darkorange', mew=1.5)
ax1.plot(dra3[gc23.index.values], ddec3[gc23.index.values], 'o', ms=8.0, fillstyle='none', mec='darkorange', mew=1.5)
sym_gc2, = ax1.plot(-100.0,-100.0, 'o', ms=8.0, fillstyle='none', mec='darkorange', mew=1.5, label='Rank 2 GCs')

ax1.plot(dra1[gc31.index.values], ddec1[gc31.index.values], 'o', ms=8.0, fillstyle='none', mec='violet', mew=1.5)
ax1.plot(dra2[gc32.index.values], ddec2[gc32.index.values], 'o', ms=8.0, fillstyle='none', mec='violet', mew=1.5)
ax1.plot(dra3[gc33.index.values], ddec3[gc33.index.values], 'o', ms=8.0, fillstyle='none', mec='violet', mew=1.5)
sym_gc3, = ax1.plot(-100.0,-100.0, 'o', ms=8.0, fillstyle='none', mec='violet', mew=1.5, label='Rank 3 GCs')

ax1.plot(dra1[sc1.index.values], ddec1[sc1.index.values], 'D', ms=8.0, fillstyle='none', mec='blueviolet', mew=1.5)
ax1.plot(dra2[sc2.index.values], ddec2[sc2.index.values], 'D', ms=8.0, fillstyle='none', mec='blueviolet', mew=1.5)
ax1.plot(dra3[sc3.index.values], ddec3[sc3.index.values], 'D', ms=8.0, fillstyle='none', mec='blueviolet', mew=1.5)
sym_sc, = ax1.plot(-100.0,-100.0, 'D', ms=8.0, fillstyle='none', mec='blueviolet', mew=1.5, label='New YSCs')


# # ----- Isolated GCs ----- #
# ax1.plot(-0.1851, 0.4562, 'D', ms=8.0, color='magenta', mew=1.0)
# # ax1.plot(-0.1927, 0.5893, 'D', ms=8.0, color='magenta', mew=1.0)
# sym_igc, = ax1.plot(-100.0,-100.0, 'D', ms=8.0, color='magenta', mew=1.0, label='IGCs (Jang+12)')


# Counting
print('M81 Group galaxies : M81, M82, NGC 3077, HolmIX, BK3N')
print('Rank 1 GCs : {0:d}'.format(len(gcc1.index.values)+len(gcc2.index.values)+len(gcc3.index.values)))
print('Rank 2 GCs : {0:d}'.format(len(gc21.index.values)+len(gc22.index.values)+len(gc23.index.values)))
print('Rank 3 GCs : {0:d}'.format(len(gc31.index.values)+len(gc32.index.values)+len(gc33.index.values)))
print('Young star clusters : {0:d}'.format(len(sc1.index.values)+len(sc2.index.values)+len(sc3.index.values)))
print('Known GCs : {0:d}'.format(len(kgc1.index.values)+len(kgc2.index.values)+len(kgc3.index.values)))
print('Young star clusters (Lim+13) : {0:d}'.format(len(ysc1.index.values)+len(ysc2.index.values)+len(ysc3.index.values)))
print('SNRs : {0:d}'.format(len(snr1.index.values)+len(snr2.index.values)+len(snr3.index.values)))
print('PNe/HII regions : {0:d}'.format(len(pnh1.index.values)+len(pnh2.index.values)+len(pnh3.index.values)))
print('X-ray sources : {0:d}'.format(len(x1.index.values)+len(x2.index.values)+len(x3.index.values)))
print('Background galaxies : {0:d}'.format(len(gal1.index.values)+len(gal2.index.values)+len(gal3.index.values)))
print('Field stars : {0:d}'.format(len(star1.index.values)+len(star2.index.values)+len(star3.index.values)))


# ----- Legends ----- #
ax1.legend(handles=[sym_gcc, sym_gc2, sym_gc3, sym_sc, sym_kgc, sym_ysc, sym_snr, sym_pnh, sym_x, sym_gal, sym_star],
	       fontsize=15.0, loc=(1.025, 0.55),
	       handlelength=0, numpoints=1, frameon=True, borderpad=0.7, framealpha=0.8, edgecolor='gray')

# ----- Scale bar ----- #
ax1.arrow(-0.15, -0.8, -0.46477, 0.0, width=0.02, head_width=0., head_length=0.,
          fc='k', ec='k', alpha=0.8)
ax1.text(-0.25, -0.75, '30 kpc', fontsize=18.0, fontweight='bold', color='k')

# # ----- Galaxy name ----- #
# ax1.text(-0.06, -0.12, 'M81', fontsize=18.0, fontweight='bold', color='darkgreen')
# ax1.text(0.325, 0.6, 'M82', fontsize=18.0, fontweight='bold', color='darkgreen')
# ax1.text(1.0, -0.475, 'NGC 3077', fontsize=18.0, fontweight='bold', color='darkgreen')

# ----- Title ----- #
plt.suptitle("MMT configurations for the M81 Group", fontsize=27.5, fontweight='bold')


plt.savefig("fig_spatial.png", dpi=300)
plt.savefig("fig_spatial.pdf", dpi=300)
plt.close()


# ---- Printing the running time ----- #
print("--- %s seconds ---" % (time.time() - start_time))