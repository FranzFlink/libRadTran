;cd /projekt_agmwend/data/EUREC4A/11_VELOX-Tools/
;.r VELOX_TB_EUREC4A_nav_smart.pro

;**********************************************************************
;Use this to perform the fdw simulations using the iMAR Navigation Data
;applies the merged bco and dropsonde data, which fits temporal closest
;**********************************************************************



pro sod2hhmmssarray,sod_input,hh,mm,ss
 sod = double(sod_input)
 hh = fix(sod/3600d)
 mmf = (sod/3600d - double(hh) ) * 60d
 mm = fix(mmf)
 ssf = ( mmf - float(mm) ) * 60.0
 ss = round(ssf)
 so = where(ss ge 60,mso)
 if mso gt 0 then begin
    ss[so]-=60
    mm[so]++
    endif
 mo=where(mm ge 60,mmo)
 if mmo gt 0 then begin
       mm[mo]-=60
       hh[mo]++
       endif
 end

function start_libradtran2,ts,lati,longi,alti,atmf,ozone,date,wlrange
  lib_file_path = '/projekt_agmwend/data/EUREC4A/11_VELOX-Tools/libradtran/'
  openw,1,lib_file_path+'VELOX_'+strtrim(string(wlrange),1)+'_TB'+date+'.inp
    printf,1,'data_files_path /opt/libradtran/2.0.4/share/libRadtran/data            # location of internal libRadtran data'
    printf,1,'rte_solver disort'
    printf,1,'mol_abs_param reptran medium'
    printf,1,'atmosphere_file /projekt_agmwend/data/EUREC4A/11_VELOX-Tools/add_data/afglt_evi.dat'
    printf,1,'radiosonde '+atmf+' H2O RH'
    printf,1,'time '+ts
    if float(lati) lt 0.0 then printf,1,'latitude S '+string(abs(float(lati)))
    if float(lati) ge 0.0 then printf,1,'latitude N '+string(abs(float(lati)))
    if float(longi) lt 0.0 then printf,1,'longitude W '+string(abs(float(longi)))
    if float(longi) ge 0.0 then printf,1,'longitude E '+string(abs(float(longi)))
    if alti lt 0.0 then printf,1,'zout 0 #'+alti
    if alti gt 0.0 then printf,1,'zout '+alti
    printf,1,'phi 0.0			# output azimuth angle'
    printf,1,'umu 1.0    	    # cosine of output zenith angle'; 1.0=0°, 0.9962=5°,  0.9848=10°, 0.9659=15°, 0.9397=20°

    if wlrange eq 0 then printf,1,'wavelength  7700 12000'       ;VELOX Broadband 7.7-12 mum
    if wlrange eq 1 then printf,1,'wavelength  8095  9195'       ;VELOX Filter 8645
    if wlrange eq 2 then printf,1,'wavelength 10350 11130'       ;VELOX Filter 10740
    if wlrange eq 3 then printf,1,'wavelength 10847 12000'       ;VELOX Filter  11660
    if wlrange eq 4 then printf,1,'wavelength 11500 12000'       ;VELOX Filter 11500
    printf,1,'source thermal'
;    printf,1,'filter_function_file /projekt_agmwend/data/EUREC4A/11_VELOX-Tools/add_data/filter_function_velox_ge_window.dat 		# to consider for the spectral transmission of the Germanium window'
;    printf,1,'source thermal /projekt_agmwend/data/EUREC4A/11_VELOX-Tools/add_data/UVSPEC_LOWTRAN_THERMAL.TRANS		# extraterr. solar flux'
;    printf,1,'wavelength_grid_file /projekt_agmwend/data/EUREC4A/11_VELOX-Tools/add_data/UVSPEC_LOWTRAN_THERMAL.TRANS	# wavelengths to be used'
    printf,1,'output_process per_nm'
    printf,1,'output_process integrate	# calculate broadband quantities'
    printf,1,'output_quantity brightness	# transfer into brightness temperatures'
    printf,1,'output_user lambda sza uu'
    printf,1,'albedo_library IGBP'
    printf,1,'brdf_rpv_type 17'
    printf,1,'sur_temperature 300       # Surface temperature'
    printf,1,'mol_modify O3 '+ozone + '  DU'
    printf,1,'aerosol_default'
    ;;printf,1,'verbose'
    ;printf,1,'quiet'
  close,1
  if wlrange eq 0 then cmd = '/opt/libradtran/2.0.4/bin/uvspec < '+lib_file_path + 'VELOX_0_TB'+date+'.inp > '+lib_file_path+'VELOX_0_TB'+date+'.out'
  if wlrange eq 1 then cmd = '/opt/libradtran/2.0.4/bin/uvspec < '+lib_file_path + 'VELOX_1_TB'+date+'.inp > '+lib_file_path+'VELOX_1_TB'+date+'.out'
  if wlrange eq 2 then cmd = '/opt/libradtran/2.0.4/bin/uvspec < '+lib_file_path + 'VELOX_2_TB'+date+'.inp > '+lib_file_path+'VELOX_2_TB'+date+'.out'
  if wlrange eq 3 then cmd = '/opt/libradtran/2.0.4/bin/uvspec < '+lib_file_path + 'VELOX_3_TB'+date+'.inp > '+lib_file_path+'VELOX_3_TB'+date+'.out'
  if wlrange eq 4 then cmd = '/opt/libradtran/2.0.4/bin/uvspec < '+lib_file_path + 'VELOX_4_TB'+date+'.inp > '+lib_file_path+'VELOX_4_TB'+date+'.out'
;  if wlrange eq 0 then cmd = 'uvspec < '+lib_file_path + 'VELOX_0_TB'+date+'.inp > '+lib_file_path+'VELOX_0_TB'+date+'.out'
;  if wlrange eq 1 then cmd = 'uvspec < '+lib_file_path + 'VELOX_1_TB'+date+'.inp > '+lib_file_path+'VELOX_1_TB'+date+'.out'
;  if wlrange eq 2 then cmd = 'uvspec < '+lib_file_path + 'VELOX_2_TB'+date+'.inp > '+lib_file_path+'VELOX_2_TB'+date+'.out'
;  if wlrange eq 3 then cmd = 'uvspec < '+lib_file_path + 'VELOX_3_TB'+date+'.inp > '+lib_file_path+'VELOX_3_TB'+date+'.out'
;  if wlrange eq 4 then cmd = 'uvspec < '+lib_file_path + 'VELOX_4_TB'+date+'.inp > '+lib_file_path+'VELOX_4_TB'+date+'.out'
  spawn,cmd,spawnresult,spawnerror
  return,spawnerror
end






close,/all




date = ['20200205a']

altitude = 10000
path1 = '/projekt_agmwend2/data_raw/EUREC4A_raw_only/06_Flights/Flight_'+date
lib_file_path = '/projekt_agmwend/data/EUREC4A/11_VELOX-Tools/libradtran/'
flightname = date
o3= 375
nav_inp=1  
navfile = file_search(path1+'/horidata/NavCommand/','Nav_GPSPos*')
openr,1,navfile
dummy=''
for i=0,12 do readf,1,dummy
data=fltarr(8,file_lines(navfile)-13)
readf,1,data
close,1

sod=data(1,*)  
lat=data(3,*)
lon=data(2,*)
alt=((data(4,*)*0.0)+Altitude)/1000. 
temp_res=5L

sod2hhmmssarray,sod,h,m,s

datestring = strmid(date,0,4) + ' ' + strmid(date,4,2) + ' ' + strmid(date,6,2) + ' '

TB_out=dblarr(5,n_elements(sod)+temp_res+1.)

j=0L
counter = 0L
jindex = 0L
check = long(sod[j])
for j=0.,n_elements(sod)-1. do begin

    if (j eq 0L) or (j eq n_elements(sod)-1) or (long(sod[j]) ge temp_res + check) then begin
        if (j eq n_elements(sod)-1) then s[j]++  ; make sure last second is fully covered

        timestring = datestring + string(h[j],format='(I2.2)') + ' ' + string(m[j],format='(I2.2)') + ' ' + string(s[j],format='(I2.2)')
        atmos_path='/projekt_agmwend/data/EUREC4A/03_Soundings/RS_for_libradtran/Merged/'+strmid(date,4,4)+'/'
        atmos_array=file_search(atmos_path,'*.dat')
        atmos_index=where(abs(strmid(atmos_array,11,5,/reverse)-sod(j)) eq min(abs(strmid(atmos_array,11,5,/reverse)-sod(j))))
        atmosfile = atmos_array(atmos_index)
        print,'SOD: '+sod(j)
        print,'using radiosonde: '+atmosfile
        wait,0.1

        waste = start_libradtran2(timestring,string(lat[j]),string(lon[j]),string(alt[j]),atmosfile,string(o3),date,'0')

        result4=0            
        wait,0.1       
        
        for zt=0,4 do begin
        ofile=lib_file_path+'VELOX_'+strtrim(string(zt),1)+'_TB'+date+'.out'
        if (j eq 0) then begin
            nwl=file_lines(ofile)
            print,file_lines(ofile)
            data0 = fltarr(3,nwl)
        endif
        openr,2,ofile
        readf,2,data0
        close,2
        TBdummy = data0[2,*] 			; brightness Temperature
        wl = reform(data0[0,*])
        TB_out(zt,j:j+temp_res) = TBdummy

        endfor

        jindex = [jindex, j]
        counter++
        check = long(sod[j])
endif
    if (long(sod[j]) gt temp_res + check)  then check = long(sod[j])
    j++
print,j,n_elements(sod)
print,'Completed '+string(j/n_elements(sod)*100.)+' % of that day.'
endfor

jindex=jindex[1:*]
print,'Data output for '+flightname

outfile = '/projekt_agmwend/data/EUREC4A/06_Flights/Flight_'+date+'/VELOX_cloud-free_TB_'+date+'_R0_ds_REPTRAN-MEDIUM_5sec_BRDF-17_zout_'+strtrim(string(altitude),1)+'m.dat'

openw,oo,outfile,/GET_LUN

    datatitle='sod         Lat(N)       Lon(E)         TB-VELOX-Full    TB-VELOX-8645    TB-VELOX-10740    TB-VELOX-11660    TB-VELOX-11500'
    printf,oo,datatitle
    for k=0L, counter-1 do begin
        j=jindex[k]
        printf,oo,sod[j],lat[j],lon[j],TB_out[0,j],TB_out[1,j],TB_out[2,j],TB_out[3,j],TB_out[4,j],format='(I5.5,2F13.6,5F15.7)'
    endfor ; k
    alt = altitude
    close,oo & free_lun,oo
  

endfor ;end loop over days

print,'Good Night and Good Luck!'

end
