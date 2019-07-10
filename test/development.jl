#=
#rigi according to 2016 paper
y=679520.05
x=212273.44
#using BenchmarkTools 
@btime latLong($x,$y)

@code_warntype latLong(x,y)
lat,long=latLong(x,y)
latA,longA=approxLatLong(x,y)
lat_seconds=((lat-47)*60-3)*60
long_seconds=((long-8)*60-29)*60
latA_seconds=((latA-47)*60-3)*60
longA_seconds=((longA-8)*60-29)*60
-3*60


latLong(x,y)
a,b=latLong(x,y)
a,b=latLong(200_000,600_000)

a,b=latLong(-100_000,100_000)
compound(a),compound(b)


a,b=latLong(191775.030,602030.680)
compound(a),compound(b)
#we have CH1903+ ellipsoid coordinates
X,Y,Z=geoCentricCartesianCoordinates(a,b,897.361)
#X,Y,Z=geoCentricCartesianCoordinates(a,b,500)

#convert to ETRS89
Xnew=X+674.374
Ynew=Y+15.056
Znew=Z+405.346

aNew,bNew,hNew=geoCentricCartesianToEllipsoidCoordinates(Xnew,Ynew,Znew)
compound(aNew),compound(bNew)

@code_warntype geoCentricCartesianToEllipsoidCoordinates(Xnew,Ynew,Znew)
@btime geoCentricCartesianToEllipsoidCoordinates($Xnew,$Ynew,$Znew)

xx=191775.030
yy=602030.680
@btime LV03toETRS89($xx,$yy)
@code_warntype LV03toETRS89(xx,yy)

=#

#=
#Bern, Sidlerstrasse 5
#approxCHCoord(lat,long)
#http://geodesy.geo.admin.ch/reframe/lv03towgs84?easting=600000&northing=200000&altitude=550.0&format=json

lat,long = 46.95108288705888,7.438632502714563
chCoord(lat,long)

y,x=600000,200000
lat_res,long_res=approxLatLong(x,y)
lat-lat_res
long-long_res

approxLatLong(200_000,600_000)

=#

#=
    @btime approxCHCoord($lat,$long)
    @code_warntype approxCHCoord(lat,long)

    x,y=approxCHCoord(lat,long)
    approxLatLong(x,y)

    @btime approxLatLong($x,$y)
    @code_warntype approxLatLong(x,y)
=#

#=
#Rigi
#47°03'28.95659233"N 8°29'11.11127154"E
lat=47+(3*60+28.95659233)/3600
long=8+(29*60+11.11127154)/3600
#y,x=679602.034,212422.148
x,y=chCoord(lat,long)
lat_res,long_res=approxLatLong(x,y)
lat-lat_res
long-long_res
x-chCoord(lat_res,long_res)[1]
y-chCoord(lat_res,long_res)[2]
=#
#=
    @code_warntype chCoord(lat,long)
    @btime chCoord($lat,$long) #198.879ns (199.165,199.165)
=#