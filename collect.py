#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:22:35 2020

@author: jlee
"""


import numpy as np
import pandas as pd


# ----- GSC-2 catalog ----- #
di_gsc = 'GSC2/'

## Reading GSC-2 catalog
gsc_name = di_gsc+'GSC_240_ver1.csv'
df_gsc2 = pd.read_csv(gsc_name, comment='#', header=0, usecols=(0, 2, 3, 4, 27, 28), na_values=' ')

## V magnitude condition
vmag_cnd = ((df_gsc2['epoch'] >= 2000.0) & (df_gsc2['VMag'] > 14.0) & (df_gsc2['VMag'] < 16.0))


## Making new data frame
df_gsc2_cut = df_gsc2.loc[vmag_cnd]
df_gsc2_cut.to_pickle("./df_gsc2_cut.pkl")


# ----- M81 SNRs (Lee+15) ----- #
di_SNR = 'SNRs/'
df_L15 = pd.read_csv(di_SNR+'Lee+15_tab2.dat', comment='#', header=None, usecols=(0,5,6),
                     names=('ID','RA','Decl'))
df_L15.to_pickle("./df_L15.pkl")


# ----- M81 PNe & HII regions (Stanghellini+10) ----- #
di_PNH = 'PNe+HII/'
df_S10 = pd.read_csv(di_PNH+'Stanghellini+10_tab1.dat', header=None,
                     names=('ID','RA','Decl'))
df_S10.to_pickle("./df_S10.pkl")


# ----- M81 HII regions (Patterson+12) ----- #
df_P12 = pd.read_csv(di_PNH+'Patterson+12_tab1.dat', header=None,
                     usecols=(0,1,2), names=('ID','RA','Decl'))
df_P12.to_pickle("./df_P12.pkl")


# ----- M82 PNe (Johnson+09) ----- #
df_J09 = pd.read_csv(di_PNH+'Johnson+09_tab3.dat', header=None,
                     usecols=(0,3,4,5,6,7,8),
                     names=('ID','RAh','RAm','RAs','DEd','DEm','DEs'))
df_J09.to_pickle("./df_J09.pkl")


# ----- NGC 3077 SNRs (Leonidaki+13) ----- #
df_L13 = pd.read_csv(di_SNR+'Leonidaki+13_tab4.dat', header=None,
                     names=('ID','RA','Decl','Flag'))
df_L13.to_pickle("./df_L13.pkl")


# ----- X-ray sources (Liu & Bregman 2005) ----- #
diX = 'X-ray/'
df_LB05 = pd.read_csv(diX+'LB05_table2.dat', header=None, usecols=(0,1,2,3,4,5,6),
                      names=('ID','RAh','RAm','RAs','DEd','DEm','DEs'))
df_LB05.to_pickle("./df_LB05.pkl")


# ----- M82 star clusters (Lim+13) ----- #
diY = 'ySCs/'
df_Lim13 = pd.read_fwf(diY+'Lim+13_table2.dat', header=None, skiprows=32, usecols=(1,2,3,4,5,6,7,8),
                       names=('RAh','RAm','RAs','DEd','DEm','DEs','Vmag','e_Vmag'))

## V magnitude condition
vmag_cnd = (df_Lim13['Vmag'] < 19.0)

## Making new data frame
df_Lim13_cut = df_Lim13.loc[vmag_cnd]
df_Lim13_cut.to_pickle("./df_Lim13_cut.pkl")


# ----- SDSS field stars ----- #
di_SDSS = 'SDSS/'
df_SDSS_star = pd.read_csv(di_SDSS+'SQL_star.csv', comment='#', usecols=(0,1,3,4,5,11))

## r magnitude condition
rmag_cnd = ((df_SDSS_star['psfMag_r'] > 14.5) & \
            (df_SDSS_star['psfMag_i'] < 20.0) & \
            (df_SDSS_star['psfMag_g'] - df_SDSS_star['psfMag_i'] > 0.2) & \
            (df_SDSS_star['psfMag_g'] - df_SDSS_star['psfMag_i'] < 1.4) & \
            (df_SDSS_star['specObjID'] == 0))

## Making new data frame
df_SDSS_star_cut = df_SDSS_star.loc[rmag_cnd]
df_SDSS_star_cut.to_pickle("./df_SDSS_star_cut.pkl")


# ----- SDSS background galaxies ----- #
di_SDSS = 'SDSS/'
df_SDSS_galaxy = pd.read_csv(di_SDSS+'SQL_galaxy.csv', comment='#', usecols=(0,1,3,4,5,11))

## r magnitude condition
imag_cnd = ((df_SDSS_galaxy['fiber2Mag_i'] > 15.0) & \
            (df_SDSS_galaxy['fiber2Mag_i'] < 20.0) & \
            (df_SDSS_galaxy['specObjID'] == 0))

## Making new data frame
df_SDSS_galaxy_cut = df_SDSS_galaxy.loc[imag_cnd]
df_SDSS_galaxy_cut.to_pickle("./df_SDSS_galaxy_cut.pkl")
