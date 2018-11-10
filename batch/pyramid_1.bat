@if [%4]==[] @goto error
@set /a halfH=%1/2
@set /a halfW=%2/2
@set decomp=%~4
:: if %3==0 then select 5x3 wavelet transform, otherwise default to 9x7
@if %3==0 (set kdu_string= Catk=2 Kextension:I2=SYM Kreversible:I2=no Ksteps:I2={2,0,0,0},{2,-1,0,0} Kcoeffs:I2=-0.5,-0.5,0.25,0.25) else (set kdu_string=)
cd tmp
kdu_compress -i coarse.rawl,diff0.rawl -o out.j2c Cdecomp="%decomp%,B(-:-:-)" -rate 1,0.8,0.6,0.4,0.2,0.1 Sdims={%halfH%,%halfW%},{%1,%2} Cweight:C0=2 Cweight:C1=1 Sprecision=16 Ssigned=yes Qstep=0.001 -precise %kdu_string%
kdu_expand -i out.j2c -o coarse1.0.rawl,diff0_1.0.rawl -layers 6
kdu_expand -i out.j2c -o coarse0.8.rawl,diff0_0.8.rawl -layers 5
kdu_expand -i out.j2c -o coarse0.6.rawl,diff0_0.6.rawl -layers 4
kdu_expand -i out.j2c -o coarse0.4.rawl,diff0_0.4.rawl -layers 3
kdu_expand -i out.j2c -o coarse0.2.rawl,diff0_0.2.rawl -layers 2
kdu_expand -i out.j2c -o coarse0.1.rawl,diff0_0.1.rawl -layers 1
cd ..
@goto end
:error
REM usage: pyramid_1.bat <h> <w> <dwttype> <Cdecomp>, where h and w are the height and width of the original image, and dwttype must be 0 or 1
:end