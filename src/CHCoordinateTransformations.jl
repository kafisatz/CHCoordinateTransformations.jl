module CHCoordinateTransformations

#Coordinate conversion between ellipsoid coordinates lat/long (degrees) and Swiss (military) coordinates (in Julia)
#see https://www.swisstopo.admin.ch/content/swisstopo-internet/de/online/calculation-services/_jcr_content/contentPar/tabs/items/dokumente_und_publik/tabPar/downloadlist/downloadItems/8_1467103085694.download/refsys_d.pdf 
#see also https://tools.wmflabs.org/geohack/geohack.php?params=46.951081_N_7.438637_E_dim:1_region:CH-BE_type:landmark&pagename=Schweizer_Landeskoordinaten&language=de&title=Fundamentalpunkt+der+Schweizer+Landeskoordinaten
#see also https://www.swisstopo.admin.ch/de/karten-daten-online/calculation-services/navref.html 
#http://geodesy.geo.admin.ch/reframe/lv95towgs84?easting=2600000&northing=1200000&altitude=550.0&format=json 

export approxCHCoord 
export LV03toETRS89 
export latLong 
export geoCentricCartesianCoordinates
export geoCentricCartesianToEllipsoidCoordinates
export approxLatLong
export chCoord
export compound

function approxCHCoord(lat,long)
    phiPrime=(lat*3600-169028.66)/10000.0
    lambdaPrime=(long*3600-26782.5)/10000.0
    x = 200147.07 + 308807.95 * phiPrime + 3745.25 * lambdaPrime^2 + 76.63 * phiPrime^2 + 119.79 * phiPrime^3 - 194.56 *lambdaPrime*lambdaPrime*phiPrime
    y = 600072.37 + 211455.93 * lambdaPrime - 10938.51 * lambdaPrime * phiPrime - 0.36 * lambdaPrime *phiPrime^2 - 44.54 * lambdaPrime^3
    return x,y
end

function approxLatLong(x,y)
    #x,y are Swiss military coordinates in meters, e.g. x=200000;y=600000
    #Y,X are in 1000 kilometers (shifted by 600_000 m and 200_000 meters respectively) 
    #output phi,labmda is in degrees (lat,long)
    Y=y/1e6 - 0.6
    X=x/1e6 - 0.2
    a1=4.72973056+X*0.7925714+X^2*0.132812+X^3*0.0255+X^4*0.0048 
    a3=-0.044270-X*0.0255-X^2*0.0096
    a5=0.00096

    p0=X*3.23864877-X^2*0.0025486-X^3*0.013245+X^4*0.000048
    p2=-0.27135379-X*0.0450442-X^2*0.007553-X^3*0.00146
    p4=0.002442+X*0.00132

    phi=(16.902866+p0+p2*Y^2+p4*Y^4)*1e4/3600
    lambda=(2.67825+a1*Y+a3*Y^3+a5*Y^5)*1e4/3600

    #within Switzerland the accuracy should be 0.00014 for labmda and 0.0004 for phi
    return phi,lambda
end


function chCoord(lat,long)
    #see: https://cms.geo.admin.ch/www.swisstopo.admin.ch/archives/cms2007/internet/swisstopo/de/home/topics/survey/sys/refsys/switzerland.parsysrelated1.24280.downloadList.42086.DownloadFile.tmp/refsysd.www.pdf
    #lat long are epxected to be latitude and longitude cooridnates
    #e.g. lat=47.05804;long=8.486419

    #input conversion to rad
    lambda=long*pi/180
    phi=lat*pi/180
    #constants
    e=Base.MathConstants.e
    alpha=1.00072913843038
    b0=(46+(54*60+27.83324844)/3600)*pi/180
    K=0.0030667323772751
    R=6378815.90365
    E=0.08169683121525584 #sqrt(0.006674372230614)
    lambda0=(7+(26*60+22.5)/3600)*pi/180
    phi0=(46+(57*60+8.66)/3600)*pi/180

    s=alpha*log(tan(pi/4+phi/2))-alpha*E/2*log((1+E*sin(phi))/(1-E*sin(phi)))+K
    b=2*(atan(e^s)-pi/4)
    I=alpha*(lambda-lambda0)
    Ibar=atan(sin(I)/(sin(b0)*tan(b)+cos(b0)*cos(I)))
    bbar=asin(cos(b0)*sin(b)-sin(b0)*cos(b)*cos(I))
    y=R*Ibar+600000
    x=R/2*log((1+sin(bbar))/(1-sin(bbar)))+200000
    return x,y
end

function latLong(x,y)
    #x,y are Swiss military coordinates in meters, e.g. x=200000;y=600000
    #Y,X are shifted by 600_000 m and 200_000 meters respectively
    #the output (phi,labmda) are CH1903+ (UTM Zone 32) coordinates in degrees (lat,long)
    Y=y-600_000
    X=x-200_000

    #constants
    e=Base.MathConstants.e
    alpha=1.00072913843038
    b0=(46+(54*60+27.83324844)/3600)*pi/180
    K=0.0030667323772751
    R=6378815.90365
    E=0.08169683121525584 #sqrt(0.006674372230614)
    lambda0=(7+(26*60+22.5)/3600)*pi/180
    phi0=(46+(57*60+8.66)/3600)*pi/180

    Ibar=Y/R
    bbar=2*(atan(exp(X/R))-pi/4)
    b=asin(cos(b0)*sin(bbar)+sin(b0)*cos(bbar)*cos(Ibar))
    I=atan(sin(Ibar)/(cos(b0)*cos(Ibar)-sin(b0)*tan(bbar)))
    lambda=lambda0+I/alpha
    #phi=BigFloat(b)
    phi=b
    counter=0
    while true
        phi0=phi
        #S=(log(tan(pi/4+BigFloat(b)/2))-BigFloat(K))/BigFloat(alpha)+BigFloat(E)*log(tan(pi/4+asin(BigFloat(E)*sin(phi)/2)))
        S=(log(tan(pi/4+(b)/2))-(K))/(alpha)+(E)*log(tan(pi/4+asin((E)*sin(phi))/2))
        #S=log(tan(pi/4+phi/2))
        phi=2*atan(exp(S))-pi/2
        counter+=1
        #@show counter
        if counter>20
            break
        end
        #if (counter>1)&&(phi==phi0)
        if (phi==phi0)
            break
        end        
    end

    #convert to degree
    lambda=lambda*180/pi
    phi=phi*180/pi
return phi,lambda
end


function compound(x::Float64)
    m1=trunc(Int,x)
    m2=trunc(Int,(x-m1)*60)
    m3=(x-m1-m2/60)*3600
    return m1,m2,m3
end


function geoCentricCartesianCoordinates(lat,long,h=0)    
    #input needs to be in degree
    phi=lat/180*pi
    lambda=long/180*pi
    #a=6378137.000 #grs80
    #b=6356752.314140 #grs80
    a=6377397.155 #bessel 1841
    b=6356078.962822 #bessel 1841
    e=sqrt(a^2-b^2)/a
    RN=a/sqrt(1-e^2*(sin(phi))^2)
    X=(RN+h)*cos(phi)*cos(lambda)
    Y=(RN+h)*cos(phi)*sin(lambda)
    Z=(RN*(1-e^2)+h)*sin(phi)
    return X,Y,Z
end

function geoCentricCartesianToEllipsoidCoordinates(x,y,z)
    a=6378137.000 #grs80
    b=6356752.314140 #grs80
    #a=6377397.155 #bessel 1841
    #b=6356078.962822 #bessel 1841
    e=sqrt(a^2-b^2)/a
    
    lambda=atan(y/x)

    phi0=atan(z/sqrt(x^2+y^2))
    phi=phi0

    RN=a/sqrt(1-e^2*(sin(phi))^2)
    h=sqrt(x^2+y^2)/cos(phi)-RN
    phi=atan((z/sqrt(x^2+y^2))/(1-(RN*e^2)/(RN+h)))

    counter=0 
    while true
        phi0=phi
        RN=a/sqrt(1-e^2*(sin(phi))^2)
        h=sqrt(x^2+y^2)/cos(phi)-RN
        phi=atan((z/sqrt(x^2+y^2))/(1-(RN*e^2)/(RN+h)))
        counter+=1
        #@show counter
        if counter>20
            break
        end
        #if (counter>1)&&(phi==phi0)
        if (phi==phi0)
            break
        end        
    end

    lambda=lambda*180/pi
    phi=phi*180/pi
    return phi,lambda,h
end

function LV03toETRS89(east,north,h=0)
    a,b=latLong(north,east)
    #compound(a),compound(b)
    #we have CH1903+ ellipsoid coordinates

    X,Y,Z=geoCentricCartesianCoordinates(a,b,h)
    #convert to ETRS89
    Xnew=X+674.374
    Ynew=Y+15.056
    Znew=Z+405.346
    
    aNew,bNew,hNew=geoCentricCartesianToEllipsoidCoordinates(Xnew,Ynew,Znew)
    return aNew,bNew,hNew
end

end # module
