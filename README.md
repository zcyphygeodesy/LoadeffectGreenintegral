## Fortran codes for computation of residual surface load effects by Green's Integral
https://www.zcyphygeodesy.com/en/h-nd-134.html
## [Algorithm purpose]
    From the regional residual equivalent water height (EWH) variation grid (cm), compute the residual surface load effects on the geoid or height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (mE) and horizontal gravity gradient (NW, to the north and to the west, mE) by load Green's function integral.
    When computing the load effects of sea level variations, the height of the calculation point is the normal or orthometric height. When computing the load effects of surface atmosphere or land water variations, the height of the calculation point is the height relative to the Earth’s surface.
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAgtbbQuQYolb-olQIwlg44ugk.jpg)
## [Geophysical models]
    The Green function file LoadGreen.txt of the load indirect effect on all-element geodetic variations.
## [Main program for test entrance]
    LoadeffectGreenintegral.f90
    The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), height (m) relative to the Earth’s surface …
    The output parameters tdn(14): the residual surface load effects on all-element geodetic variations.
tdn(1:14) stores the residual surface load effects on 10 kinds of geodetic variations, which are the residual load effects on height anomaly tdn(1) (mm), ground gravity #tdn(2) (μGal), gravity disturbance tdn(3) (μGal), ground tilt #tdn(4:5) (SW, to the south and to the west, mas), vertical deflection tdn(6:7) (SW, to the south and to the west, mas), horizontal displacement #tdn(8:9) (EN, to the east and to the north, mm), ground radial displacement #tdn(10) (mm), ground normal or orthometric height #tdn(11) (mm), radial gravity gradient tdn(12 )(mE) and horizontal gravity gradient tdn(13:14) (NW, to the north and to the west, mE).
    The calculation point can be on the ground, low altitude, ocean or underwater space. The geodetic variations abvove marked with # are valid only when the site is fixed with the solid Earth.
## (1) Reading module of the Green function file of load indirect effect
    LGrnFunc(loadgrfl,GF)
    Input parameters: loadgrfl - The Green function file name LoadGreen.txt of the load indirect effect.
    Return parameters: GF(8000,9) - The Green functions of the load indirect effects on height anomaly (e-13), ground gravity (e-17), gravity disturbance (e-18), ground tilt (e-14), vertical deflection (e-19), horizontal displacement (e-12), ground radial displacement (e-11), radial gravity gradient (e-15) and horizontal gravity gradient (e-15). Where the integral distance of GF(i,1:9) is equal to 100i (m).
## (2) Computation module for residual surface load effects by Green's Integral
    rntGreenintegral(BLH,ewh,hd,nlat,nlon,GF,direct,indrct,GRS,dr)
    Input parameters: BLH(3) - longitude (decimal degrees), latitude (decimal degrees), height (m) of the calculation point relative to the Earth’s surface. 
    Input parameters: ewh(nlat,nlon) - the regional residual equivalent water height (EWH) variation grid (cm/hPa).
    Input parameters: dr, hd(6) - the integral radius (m) and grid specification parameters (minimum and maximum longitude, minimum and maximum latitude, longitude and latitude intervals of a cell grid).
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    When GRS(6) = -1, the atmospheric load effect is calculated. When GRS(6) > 0, the load effect of land water, sea level variation or the sum of the two is calculated.
    Return parameters: direct(10) - the residual load direct effect on the height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), radial gravity gradient (1mE) and horizontal gravity gradient (NW, to the north and to the west, mE).
    Return parameters: indrct(14)  - the residual load indirect effect on on the height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (1mE) and horizontal gravity gradient (NW, to the north and to the west, mE). Where Indrct(11) = default value for the normal or orthometric height variation.
## (3) Calculation module for the normal gravity field
    normdjn(GRS,djn); GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    RLAT_BLH(GRS, RLAT, BLH)
## (6) Other auxiliary modules
    IntpGrnF(GF,dl,vfn); PickRecord(str0, kln, rec, nn)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] PAGravf4.5 User Reference https://www.zcyphygeodesy.com/en/
    7.2.3 The geodetic numerical grid file
    8.3 Surface load effects on various geodetic variations by Green's Integral
The zip compression package includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable test file and all input and output data.
