:: Batch file for calling kdu_compress on a temporary image file out.rawl
:: This program must be run with the following parameters:
:: out.bat <height> <width> <dwttype> <Cdecomp> <layer>
:: where height, width are the dimensions of the image, dwttype is 0 if the 5x3 wavelet is desired or 1 if the 9x7 wavelet is desired, and layer specifies the depth of dyadic analysis that the image has been extracted from.

:: set up variables
@echo off
:: TODO: goto error if number of parameters is wrong

:: if dwttype==0 then select 5x3 wavelet transform, otherwise default to 9x7
if %3==0 (set kdu_string= Catk=2 Kextension:I2=SYM Kreversible:I2=no Ksteps:I2={2,0,0,0},{2,-1,0,0} Kcoeffs:I2=-0.5,-0.5,0.25,0.25) else (set kdu_string=)
set decomp=%~4
set rate0=1,0.8,0.6,0.4,0.2,0.1
set rate1=1.4,1.2,1,0.8,0.6,0.4
if %5==0 (set rate=%rate0%) else (set rate=%rate1%)
@echo on

cd tmp
kdu_compress -i out.rawl -o out.j2c Cdecomp="%decomp%,B(-:-:-)" -rate %rate% Sdims={%1,%2} Sprecision=16 Ssigned=yes Qstep=0.001 -precise %kdu_string%
kdu_expand -i out.j2c -o out1.0.rawl -layers 6
kdu_expand -i out.j2c -o out0.8.rawl -layers 5
kdu_expand -i out.j2c -o out0.6.rawl -layers 4
kdu_expand -i out.j2c -o out0.4.rawl -layers 3
kdu_expand -i out.j2c -o out0.2.rawl -layers 2
kdu_expand -i out.j2c -o out0.1.rawl -layers 1
cd ..
@goto end

:error
REM usage: out.bat <h> <w> <dwttype> <Cdecomp>, where h and w are the height and width of the full resolution image, and dwttype must be 0 or 1
:end