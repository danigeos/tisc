#!/bin/csh

if ($6 == "") then
	echo "Syntax:    $0 projectname  lonmin lonmax  latmin latmax  resolution"
	echo "Example:   make_n_run_present_tisc_prj.csh test -10 1 36 41 10m"
	exit
endif

set prj = $1

set lonmin = $2
set lonmax = $3
set latmin = $4
set latmax = $5
set resol_model = $6

set projection = M8 #S55/35/8 #B50/30/30/45/17.5 #Poly/42/30/17.5 #B42/30/30/45/17.5
set iluminangle	= 0
set shadfac = .7

#############################
set ps = $prj.map.ps
set tmp = tmp.$prj

set Region_data	= $lonmin/$lonmax/$latmin/$latmax
set Region_plot	= $Region_data #60/20/115/40\r #60/20/125/34\r #41/20/137/39\r

echo Region_data = $Region_data
echo Region_plot = $Region_plot

set mean_lat = `awk -v a=$latmin -v b=$latmax 'BEGIN{print (b+a)/2}'`
lonlat2xy -L$mean_lat -x$lonmin -y$latmin <<END> $tmp.length.tmp
	$lonmin $latmin
	$lonmax $latmax
END
set length = `awk '{if (NR==2) print $1,$2}' $tmp.length.tmp `
echo ; echo Dimensions: $length[1] x $length[2] 


###############################################################################
echo TOPO...
set grdfile = ~/DGC_global_data/topography/ETOPO1/ETOPO1_Bed_g.grd
psbasemap -J$projection -R$Region_plot -Ba5nseW -X2 -Y14 -K -P > $ps

grdcut $grdfile -G$tmp.topobat_cont.grd.tmp -R$Region_data
grdmath $tmp.topobat_cont.grd.tmp 1000 DIV = $tmp.topobat_cont.grd.tmp 
grdsample $tmp.topobat_cont.grd.tmp -G$tmp.topobat_col.grd.tmp -I$resol_model -F
#grdinfo topobat_cont.grd.tmp

#makecpt -Cglobe -Tpaleta.cpt_Tfile.tmp > paleta.cpt.tmp
#makecpt -Ctopo > paleta.cpt.tmp
grd2cpt $tmp.topobat_cont.grd.tmp -Crelief -S-1/7/1 > $tmp.paleta.cpt.tmp

grdgradient $tmp.topobat_col.grd.tmp -G$tmp.intensity.grd.tmp -A$iluminangle -Nt$shadfac
grdimage   	$tmp.topobat_col.grd.tmp \
	-I$tmp.intensity.grd.tmp \
	-J$projection -R$Region_plot -C$tmp.paleta.cpt.tmp -O -K >> $ps
#grdcontour 	topobat_cont.grd.tmp -JM -R -C1 -Z.001 -L-1/8 -W1/0 -O -K >> $ps

#Rivers data base is shifted 1mm North and 1mm East!!:
pscoast -J$projection -R -Di -W1/0 \
	-O -K >> $ps
pscoast -J$projection -R -Di -I0/3/50/50/255 -I1/2/30/30/255 -I2/2/10/10/255 -I3/1/0/0/255 -I4/1/0/0/255 -Ii/1/0/0/255 \
	-A100/0/4 \
	-C0/0/255 \
	-Lf$lonmin/$mean_lat/$mean_lat/200+lkm \
	-Bg10 -O -K >> $ps

pstext 	-J$projection -R -G0 -O -K << END >> $ps
#	-5	43.8	40 0 0 2	M a r    C a n t \337 b r i c o
END

psscale	 -C$tmp.paleta.cpt.tmp -D4/-.4/8/.3h -L -B:"Elevation (km)": -O -K >> $ps




###############################################################################
echo PRECIPITATION...
set grdfile = ~/DGC_global_data/precipitation,evaporation/global_precipitation-cai_precip2.clim.grd

psbasemap	-J$projection -R$Region_plot -Ba5nsew -X9 -Y0 -K -O >> $ps


set min_rain = .2

cat - <<END>  $tmp.raincpt.tmp 
0		255	150	0	$min_rain 255	255	0
$min_rain 	230     255     80      0.5     21      255     170
0.5     	64      255     210     1       64      191     255
1       	106     149     255     1.5     106     149     255
1.5     	149     106     255     2       149     106     255
2       	212     43      255     3       212     43      255
B       0       0       0
F       255     255     255
N       128     128     128
END

grdcut $grdfile -G$tmp.rain.grd.tmp -R$Region_data
grdsample -fg $tmp.rain.grd.tmp -G$tmp.rain.grd.tmp -I.1
grdmath -fg $tmp.rain.grd.tmp 1000 DIV = $tmp.rain.grd.tmp

grdimage $tmp.rain.grd.tmp -fg -J$projection -R -C$tmp.raincpt.tmp -O -K >> $ps
pscoast -J$projection -R -Di -W1/0 \
	#-Lf130/60/60/1000+lkm \
	-O -K >> $ps
pscoast -J$projection -R -Di -I0/3/50/50/255 -I1/2/30/30/255 -I2/2/10/10/255 -I3/1/0/0/255 -I4/1/0/0/255 -Ii/1/0/0/255 \
	-A100/0/4 \
	-C0/0/255 \
	#-Lf45/30/10/1000+lkm \
	-Bg10 -O -K >> $ps

pstext 	-J$projection -R -G0 -O -K << END >> $ps
#	-5	43.8	40 0 0 2	M a r    C a n t \337 b r i c o
END

psscale	 -C$tmp.raincpt.tmp -D4/-.4/8/.3h -B:"Precipitation (m/yr)": -L -O >> $ps


convert $ps -trim $ps\.jpg



###############################################################################
echo Preparing topo for TISC...

set Nx = `grdinfo $tmp.topobat_col.grd.tmp | awk '{if ($(NF-1)=="nx:") {print $NF; exit;}}'`
set Ny = `grdinfo $tmp.topobat_col.grd.tmp | awk '{if ($(NF-1)=="ny:") {print $NF; exit;}}'`

echo grid size: $Nx / $Ny
echo "¡¡AUTOMATICALLY WRITTEN BY SCRIPT!!" > $tmp.prm.PRM
tisc -hp >> $tmp.prm.PRM
tisc $tmp.prm -M0 -N$Nx\/$Ny -D0/$length[1]/0/$length[2] \
	-ti0 -tf1 -td.1 -te1e-3 \
	-qhydro_model=1 -p300/0/4000 \
	-qerosed_model=6 -qerodability=1e-7 -qK_river_cap=1e3 \
	-qlost_rate=.05 \
	-f2 > $prj.PRM
rm $tmp.prm*
#awk '{if ($1==modeisost) exit}' $tmp.PRM2.tmp >$prj.PRM
grdsample $tmp.topobat_col.grd.tmp -G$tmp.topo.grd.tmp -I$resol_model
#grdinfo $tmp.topo.grd.tmp
echo "¡¡AUTOMATICALLY WRITTEN BY SCRIPT!!" > $prj.ZINI
echo mode_interp 0 >> $prj.ZINI
grd2xyz $tmp.topo.grd.tmp | awk '{print $1,$2,$3*1e3}' >> $prj.ZINI

grdsample $tmp.rain.grd.tmp -G$tmp.rain.grd.tmp -I$resol_model
#grdinfo $tmp.rain.grd.tmp
echo "¡¡AUTOMATICALLY WRITTEN BY SCRIPT!!" > $prj.RAIN
echo mode_interp 0 >> $prj.RAIN
grd2xyz $tmp.rain.grd.tmp | awk '{print $1,$2,$3*1e3}' >> $prj.RAIN


cat <<END > $prj\1.UNIT
	mode_interp 4
	time 1
	time_stop 3
	gradual 1
	2200
	950e3	390e3
	1020e3	390e3
	1020e3	460e3
	950e3	460e3
END

echo Running TISC...; echo
tisc $prj -Pc


rm -f tmp.*.tmp

