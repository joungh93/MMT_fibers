#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 3 14:46:34 2020

@author: jlee
"""


import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord


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


# ----- Writing catalog ----- #
f = open('target1.cat','w')
f.write('ra\tdec\tobject\trank\ttype'+'\n')
f.write('--\t---\t------\t----\t----'+'\n')


# Galaxy centers & diffuse light region (if needed)
def write_gal(ra, dec, galname, priority, add_region=True):
	ra0, dec0 = ra, dec
	sra0 = Angle(ra0, u.deg).to_string(unit=u.hour, sep=':')
	if (len(sra0.split(':')[0]) == 1):
		prefix = '0'
	elif (len(sra0.split(':')[0]) > 1):
		prefix = ''
	sra = [prefix+sra0.split(':')[0], sra0.split(':')[1], '%.3f' %(float(sra0.split(':')[2]))]
	sra0 = ':'.join(sra)

	sdec0 = Angle(dec0, u.deg).to_string(unit=u.deg, sep=':')
	sdec = [sdec0.split(':')[0], sdec0.split(':')[1], '%.2f' %(float(sdec0.split(':')[2]))]
	sdec0 = ':'.join(sdec)

	f.write(sra0+'\t'+sdec0+'\t'+galname+'-C'+'\t'+str(priority)+'\t'+'TARGET'+'\n')

	if add_region:
		c0 = SkyCoord(ra=ra0, dec=dec0, unit=u.deg)
		direction = ['N', 'E', 'S', 'W']
		separation = 30./3600 * u.deg
		for i in np.arange(len(direction)):
			position_angle = 90.*i*u.deg
			c1 = c0.directional_offset_by(position_angle, separation)
			ra1, dec1 = c1.ra.deg, c1.dec.deg

			sra1 = Angle(ra1, u.deg).to_string(unit=u.hour, sep=':')
			if (len(sra1.split(':')[0]) == 1):
				prefix = '0'
			elif (len(sra1.split(':')[0]) > 1):
				prefix = ''
			sra = [prefix+sra1.split(':')[0], sra1.split(':')[1], '%.3f' %(float(sra1.split(':')[2]))]
			sra1 = ':'.join(sra)

			sdec1 = Angle(dec1, u.deg).to_string(unit=u.deg, sep=':')
			sdec = [sdec1.split(':')[0], sdec1.split(':')[1], '%.2f' %(float(sdec1.split(':')[2]))]
			sdec1 = ':'.join(sdec)

			f.write(sra1+'\t'+sdec1+'\t'+galname+'-'+direction[i]+'\t'+str(priority)+'\t'+'TARGET'+'\n')
				
# write_gal(ra, dec, galname, priority, add_region=True)
ra0_m81, dec0_m81 = 148.88822, 69.06529
ra0_m82, dec0_m82 = 148.96969, 69.67938
ra0_n3077, dec0_n3077 = 150.82946, 68.73392
ra0_h9, dec0_h9 = 149.383333, 69.045833
ra0_bk3n, dec0_bk3n = 148.452244, 68.968841

write_gal(ra0_m81, dec0_m81, 'M81', 1, add_region=True)
write_gal(ra0_m82, dec0_m82, 'M82', 1, add_region=True)
write_gal(ra0_n3077, dec0_n3077, 'NGC3077', 1, add_region=True)
write_gal(ra0_h9, dec0_h9, 'HolmIX', 1, add_region=False)
write_gal(ra0_bk3n, dec0_bk3n, 'BK3N', 1, add_region=False)


# GC candidates
def write_gcc(priority, gc_array_bell, gc_array_okamoto):
	ra_gc = np.concatenate([df_bel.loc[gc_array_bell]['ALPHA_J2000'].values,
		                    df_oka.loc[gc_array_okamoto]['ALPHA_J2000'].values], axis=0)
	dec_gc = np.concatenate([df_bel.loc[gc_array_bell]['DELTA_J2000'].values,
		                     df_oka.loc[gc_array_okamoto]['DELTA_J2000'].values], axis=0)
	sra_gc = Angle(ra_gc, u.deg).to_string(unit = u.hour, sep=':')
	sdec_gc = Angle(dec_gc, u.deg).to_string(unit = u.degree, sep=':')

	x_gc = np.concatenate([df_bel.loc[gc_array_bell]['X_IMAGE'].values,
		                   df_oka.loc[gc_array_okamoto]['X_IMAGE'].values], axis=0)
	y_gc = np.concatenate([df_bel.loc[gc_array_bell]['Y_IMAGE'].values,
		                   df_oka.loc[gc_array_okamoto]['Y_IMAGE'].values], axis=0)
	dist = 0.168*np.sqrt((x_gc-17142)**2 + (y_gc-17142)**2)/3600.

	for i in np.arange(len(ra_gc)):
		if (float(sra_gc[i].split(':')[0]) < 10.0):
			sra = ['0'+sra_gc[i].split(':')[0], sra_gc[i].split(':')[1], '%06.3f' %(float(sra_gc[i].split(':')[2]))]
		else:
			sra = [sra_gc[i].split(':')[0], sra_gc[i].split(':')[1], '%06.3f' %(float(sra_gc[i].split(':')[2]))]
		sra = ':'.join(sra)	

		sdec = [sdec_gc[i].split(':')[0], sdec_gc[i].split(':')[1], '%05.2f' %(float(sdec_gc[i].split(':')[2]))]
		sdec = ':'.join(sdec)

		if (priority <= 3):
			obj = "GC"
			if ((priority == 2) & (dist[i] > 40./60)):
				rank = priority + 50
			else:
				rank = priority
		if (priority == 4):
			obj = "SC"
			if (dist[i] > 40./60):
				rank = priority + 50
			else:
				rank = priority

		f.write(sra+'\t'+sdec+'\t'+obj+f"{priority:d}-{i+1:04d}"+'\t'+f'{rank:d}'+'\t'+'TARGET'+'\n')


write_gcc(1, gc1_bel, gc1_oka)
write_gcc(2, gc2_bel, gc2_oka)
write_gcc(3, gc3_bel, gc3_oka)
write_gcc(4, ysc_bel, ysc_oka)


# Known GCs (216)
df_kGC = pd.read_pickle(dir_bel+"./kGC.pkl")
# Mi_kGC = df_kGC['imag'] - 27.80
# sel_kGC = (Mi_kGC < -8.5)

sra_kGC = Angle(df_kGC['RA'], u.deg).to_string(unit = u.hour, sep=':')
sdec_kGC = Angle(df_kGC['Decl'], u.deg).to_string(unit = u.degree, sep=':')

for i in np.arange(len(df_kGC)):
	if (float(sra_kGC[i].split(':')[0]) < 10.0):
		sra = ['0'+sra_kGC[i].split(':')[0], sra_kGC[i].split(':')[1], '%06.3f' %(float(sra_kGC[i].split(':')[2]))]
	else:
		sra = [sra_kGC[i].split(':')[0], sra_kGC[i].split(':')[1], '%06.3f' %(float(sra_kGC[i].split(':')[2]))]
	sra = ':'.join(sra)	

	sdec = [sdec_kGC[i].split(':')[0], sdec_kGC[i].split(':')[1], '%05.2f' %(float(sdec_kGC[i].split(':')[2]))]
	sdec = ':'.join(sdec)
	
	f.write(sra+'\t'+sdec+'\t'+'kGC-{0:04d}'.format(i+1)+'\t'+'10'+'\t'+'TARGET'+'\n')


# NGC 3077 SNRs (Leonidaki+13)
df_L13 = pd.read_pickle('./df_L13.pkl')
for i in np.arange(len(df_L13)):
	f.write(df_L13['RA'].values[i]+'\t'+df_L13['Decl'].values[i]+'\t'+'SNR-L13-{0:02d}'.format(i+1)+'\t'+'20'+'\t'+'TARGET'+'\n')


# Young star clusters (Lim et al. (2013))
df_Lim13_cut = pd.read_pickle('./df_Lim13_cut.pkl')
for i in np.arange(len(df_Lim13_cut)):
	sra = '%02d' %(df_Lim13_cut['RAh'].values[i])+':'+'%02d' %(df_Lim13_cut['RAm'].values[i])+':'+'%06.3f' %(df_Lim13_cut['RAs'].values[i])
	sdec = '%02d' %(df_Lim13_cut['DEd'].values[i])+':'+'%02d' %(df_Lim13_cut['DEm'].values[i])+':'+'%05.2f' %(df_Lim13_cut['DEs'].values[i])
	f.write(sra+'\t'+sdec+'\t'+'ySC-Lim13-{0:03d}'.format(i+1)+'\t'+'25'+'\t'+'TARGET'+'\n')


# M81/M82 SNRs + M81 HII regions (Lee et al. (2015))
df_L15 = pd.read_pickle('./df_L15.pkl')
for i in np.arange(len(df_L15)):
	f.write(df_L15['RA'].values[i]+'\t'+df_L15['Decl'].values[i]+'\t'+'SNR-M81-{0:02d}'.format(i+1)+'\t'+'30'+'\t'+'TARGET'+'\n')

sra_HII, sdec_HII = '09:56:14.29', '68:50:20.7'
f.write(sra_HII+'\t'+sdec_HII+'\t'+'PNH-L15-{0:02d}'.format(i+1)+'\t'+'30'+'\t'+'TARGET'+'\n')


# M81 HII regions (Patterson et al. (2012))
df_P12 = pd.read_pickle('./df_P12.pkl')
for i in np.arange(len(df_P12)):
	f.write(df_P12['RA'].values[i]+'\t'+df_P12['Decl'].values[i]+'\t'+'PNH-P12-{0:02d}'.format(i+1)+'\t'+'35'+'\t'+'TARGET'+'\n')


# M82 PNe (Johnson+09)
df_J09 = pd.read_pickle('./df_J09.pkl')
for i in np.arange(len(df_J09)):
	sra = '%02d' %(df_J09['RAh'].values[i])+':'+'%02d' %(df_J09['RAm'].values[i])+':'+'%06.3f' %(df_J09['RAs'].values[i])
	sdec = '%02d' %(df_J09['DEd'].values[i])+':'+'%02d' %(df_J09['DEm'].values[i])+':'+'%05.2f' %(df_J09['DEs'].values[i])
	f.write(sra+'\t'+sdec+'\t'+'PNH-J09-{0:02d}'.format(i+1)+'\t'+'35'+'\t'+'TARGET'+'\n')


# M81 PNe + HII regions (Stanghellini+10)
df_S10 = pd.read_pickle('./df_S10.pkl')
for i in np.arange(len(df_S10)):
	f.write(df_S10['RA'].values[i]+'\t'+df_S10['Decl'].values[i]+'\t'+'PNH-S10-{0:02d}'.format(i+1)+'\t'+'40'+'\t'+'TARGET'+'\n')


# X-ray sources (Liu & Bregman (2005))
df_LB05 = pd.read_pickle('./df_LB05.pkl')
for i in np.arange(len(df_LB05)):
	sra = '%02d' %(df_LB05['RAh'].values[i])+':'+'%02d' %(df_LB05['RAm'].values[i])+':'+'%06.3f' %(df_LB05['RAs'].values[i])
	sdec = '%02d' %(df_LB05['DEd'].values[i])+':'+'%02d' %(df_LB05['DEm'].values[i])+':'+'%05.2f' %(df_LB05['DEs'].values[i])
	f.write(sra+'\t'+sdec+'\t'+'Xray-LB05-{0:02d}'.format(i+1)+'\t'+'50'+'\t'+'TARGET'+'\n')


# Background galaxies from SDSS (SDSS DR16)
df_SDSS_galaxy_cut = pd.read_pickle('./df_SDSS_galaxy_cut.pkl')
sra_sdss_galaxy = Angle(df_SDSS_galaxy_cut['ra'].values, u.deg).to_string(unit = u.hour, sep=':')
sdec_sdss_galaxy = Angle(df_SDSS_galaxy_cut['dec'].values, u.deg).to_string(unit = u.deg, sep=':')

for i in np.arange(len(df_SDSS_galaxy_cut)):
	if (float(sra_sdss_galaxy[i].split(':')[0]) < 10.0):
		sra = ['0'+sra_sdss_galaxy[i].split(':')[0], sra_sdss_galaxy[i].split(':')[1], '%06.3f' %(float(sra_sdss_galaxy[i].split(':')[2]))]
	else:
		sra = [sra_sdss_galaxy[i].split(':')[0], sra_sdss_galaxy[i].split(':')[1], '%06.3f' %(float(sra_sdss_galaxy[i].split(':')[2]))]
	sra = ':'.join(sra)	

	sdec = [sdec_sdss_galaxy[i].split(':')[0], sdec_sdss_galaxy[i].split(':')[1], '%05.2f' %(float(sdec_sdss_galaxy[i].split(':')[2]))]
	sdec = ':'.join(sdec)
	
	f.write(sra+'\t'+sdec+'\t'+'Galaxy-SDSS-{0:05d}'.format(i+1)+'\t'+'70'+'\t'+'TARGET'+'\n')


# Field stars from SDSS (SDSS DR16)
df_SDSS_star_cut = pd.read_pickle('./df_SDSS_star_cut.pkl')
sra_sdss_star = Angle(df_SDSS_star_cut['ra'].values, u.deg).to_string(unit = u.hour, sep=':')
sdec_sdss_star = Angle(df_SDSS_star_cut['dec'].values, u.deg).to_string(unit = u.deg, sep=':')

for i in np.arange(len(df_SDSS_star_cut)):
	if (float(sra_sdss_star[i].split(':')[0]) < 10.0):
		sra = ['0'+sra_sdss_star[i].split(':')[0], sra_sdss_star[i].split(':')[1], '%06.3f' %(float(sra_sdss_star[i].split(':')[2]))]
	else:
		sra = [sra_sdss_star[i].split(':')[0], sra_sdss_star[i].split(':')[1], '%06.3f' %(float(sra_sdss_star[i].split(':')[2]))]
	sra = ':'.join(sra)	

	sdec = [sdec_sdss_star[i].split(':')[0], sdec_sdss_star[i].split(':')[1], '%05.2f' %(float(sdec_sdss_star[i].split(':')[2]))]
	sdec = ':'.join(sdec)
	
	f.write(sra+'\t'+sdec+'\t'+'Star-SDSS-{0:05d}'.format(i+1)+'\t'+'99'+'\t'+'TARGET'+'\n')

# sdss_guide = ((df_SDSS_star_cut['psfMag_r'].values > 14.5) & (df_SDSS_star_cut['psfMag_r'].values < 15.5))

# for i in np.arange(np.sum(sdss_guide)):
# 	if (float(sra_sdss_star[sdss_guide][i].split(':')[0]) < 10.0):
# 		sra = ['0'+sra_sdss_star[sdss_guide][i].split(':')[0], sra_sdss_star[sdss_guide][i].split(':')[1], '%06.3f' %(float(sra_sdss_star[sdss_guide][i].split(':')[2]))]
# 	else:
# 		sra = [sra_sdss_star[sdss_guide][i].split(':')[0], sra_sdss_star[sdss_guide][i].split(':')[1], '%06.3f' %(float(sra_sdss_star[sdss_guide][i].split(':')[2]))]
# 	sra = ':'.join(sra)	

# 	sdec = [sdec_sdss_star[sdss_guide][i].split(':')[0], sdec_sdss_star[sdss_guide][i].split(':')[1], '%05.2f' %(float(sdec_sdss_star[sdss_guide][i].split(':')[2]))]
# 	sdec = ':'.join(sdec)
	
# 	f.write(sra+'\t'+sdec+'\t\t\t'+'guide'+'\n')


# guide stars
df_gsc2_cut = pd.read_pickle('./df_gsc2_cut.pkl')

sra_gsc2 = Angle(df_gsc2_cut['ra'], u.deg).to_string(unit = u.hour, sep=':')
sdec_gsc2 = Angle(df_gsc2_cut['dec'], u.deg).to_string(unit = u.degree, sep=':')

for i in np.arange(len(df_gsc2_cut)):
	if (float(sra_gsc2[i].split(':')[0]) < 10.0):
		sra = ['0'+sra_gsc2[i].split(':')[0], sra_gsc2[i].split(':')[1], '%06.3f' %(float(sra_gsc2[i].split(':')[2]))]
	else:
		sra = [sra_gsc2[i].split(':')[0], sra_gsc2[i].split(':')[1], '%06.3f' %(float(sra_gsc2[i].split(':')[2]))]
	sra = ':'.join(sra)	

	sdec = [sdec_gsc2[i].split(':')[0], sdec_gsc2[i].split(':')[1], '%05.2f' %(float(sdec_gsc2[i].split(':')[2]))]
	sdec = ':'.join(sdec)
	
	f.write(sra+'\t'+sdec+'\t\t\t'+'guide'+'\n')


f.close()