# MMT_fibers


README

@ Jeong Hwan Lee

## Prerequisites
* GSC2/    # Guide stars
* PNe+HII/    # M81 PNe & HII regions
* SDSS/    # SDSS field stars & background galaxies
* SNRs/    # M81 & NGC 3077 SNRs
* X-ray/    # X-ray sources
* ySCs/    # Young stellar clusters

xfitfibs

## Python codes
```
collect.py
mk_cat.py

(run xfitfibs)
  - Load observational catalog file (open *.cat)
  - In Field Table window, 'Fit Guides' tap:
  	- Type ra, dec, start, minutes, exptime, nexp, grating, and centerwave in the rows
  	- Fill the checkbox of the rows
  	- Begin Fit
  	(If there is any problem, check start, minutes, exptime, etc.)
  - Move to 'Fit Fibers' tap:
  	- Begin Fit
  - Check in Draw Window: Window - Catalog

*** If you re-open xfitfibs, you can load the previous configuration by opening target1.cat on the Draw Window. ***

read_cfg.py
plt_spatial.py
```

