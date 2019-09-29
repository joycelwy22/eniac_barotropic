
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           GetICs
%
%  Get the Initial Conditions for
%  the ENIAC forecasts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Extract the values from the NCEP
%  analysis and convert to PS grid.
%
%  The data were extracted from the
%  NCEP Re-analysis on a hemispheric
%  domail on a 2.5x2.5 degree grid.
%
%  The data have been downloaded from
%  the URL http://nomad3.ncep.noaa.gov/
%  cgi-bin/ftp2u_6p_r1.sh
%
%  Each single field was downloaded
%  separately. This was done for simplicity.
%  They are named like 1949010503.grb
%
%  Each input field is in a separate
%  file and they are read in and
%  processed individually.
%
%  There were four ENIAC forecasts:
%  The dates in question are (all in 1949):
%
%  (1) 0030 Z, January 5th, 1949  with 
%        verification analysis Jan 6
%  (2) 0030 Z, January 30th, 1949 with
%        verification analysis Jan 31
%  (3) 0030 Z, January 31st, 1949  with
%        verification analysis Feb 1
%  (4) 0030 Z, February 13th, 1949 with
%        verification analysis Feb 14
%
%  Thus, there are seven grib files. 
%  They are named like 1949010503.grb
%  They were run through degrib with the
%  script DeCode (on met.ucd.ie) to give
%  ascii files, with names such as
%      1949010503.asc 
%
%  The ascii files (1949010503.asc etc.)
%  are the INPUTS TO THIS PROGRAM. The
%  outputs on polar stereographic grids
%  are named, e.g., Case1-1949010503.z00
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear   %   Clear the memory

%%  Open a plotting window.
figure('Position',[200,100,800,400])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                   %%%
%%%   Input ascii files with global 2.5x2.5 degree    %%%
%%%   height analysis.                                %%%
%%%                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                   %%%
%%%                1949010503.asc                     %%%
%%%                1949010603.asc                     %%%
%%%                1949013003.asc                     %%%
%%%                1949013103.asc                     %%%
%%%                1949020103.asc                     %%%
%%%                1949021303.asc                     %%%
%%%                1949021403.asc                     %%%
%%%                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                     %%%
%%   Output files on ENIAC 19x16 PS grid                %%%
%%%                                                     %%%
%%% Ncase   Analysis file         Verification file     %%%
%%%  1    Case1-1949010503.z00   Case1-1949010603.z00   %%%
%%%  2    Case2-1949013003.z00   Case2-1949013103.z00   %%%
%%%  3    Case2-1949013103.z00   Case2-1949020103.z00   %%%
%%%  4    Case2-1949021303.z00   Case2-1949021403.z00   %%%
%%%                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify the Forecast Case

        Ncase = 1;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% File name tag.
Case(1,1:6) = 'Case1-';
Case(2,1:6) = 'Case2-';
Case(3,1:6) = 'Case3-';
Case(4,1:6) = 'Case4-';

NUMx = [ 19 19 19 19 ];
NUMy = [ 16 16 16 16 ];
Xpol = [ 10 10 10 10 ];
Ypol = [ 14 14 14 14 ];
DELs = [ 736E+3 736E+3 736E+3 736E+3];
Cang = [ -85 -85 -85 -85 ];

Nx = NUMx(Ncase);
Ny = NUMy(Ncase);
Xp = Xpol(Ncase);
Yp = Ypol(Ncase);
Ds = DELs(Ncase);
cenang = Cang(Ncase);

%%  Adjustment to avoid trouble at date-line
cenang = cenang + 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Read in the analysis data

YMDH1(1,1:10) = '2012102200';
YMDH1(2,1:10) = '2012102600';
YMDH1(3,1:10) = '2012102700';
YMDH1(4,1:10) = '2012102800';

YMDH2(1,1:10) = '2012102412';
YMDH2(2,1:10) = '2012102700';
YMDH2(3,1:10) = '2012102800';
YMDH2(4,1:10) = '2012102900';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Initial analysis
%% Filein1 = [YMDH1(Ncase,1:10) '.asc'];
%% 
%% %% Verification analysis
%% Filein2 = [YMDH2(Ncase,1:10) '.asc'];
%% 
%% %% Initial analysis on PS grid
%% Fileout1 = [Case(Ncase,:) YMDH1(Ncase,1:10) '.z00']
%% 
%% %% Verification analysis on PS grid
%% Fileout2 = [Case(Ncase,:) YMDH2(Ncase,1:10) '.z00']

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Grid size for Northern Hemisphere.
Nlon = 144;   %   2.5 degree NCEP grid
Nlat =  73;   %   2.5 degree NCEP grid
Dlon = 2.5;
Dlat = 2.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %%  Read in the two analyses.
%%  Zread = zeros(Nlon*Nlat,1);
%%  
%%  Zread = load(Filein1);
%%  Z1 = reshape(Zread,Nlon,Nlat);
%%  
%%  Zread = load(Filein2);
%%  Z2 = reshape(Zread,Nlon,Nlat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Read in the analyses on the global grid.
%%  Zread1 = zeros(10512,1);
if      (Ncase == 1 ) 
  Zread1 = load('2012102200.asc');
  Zread2 = load('2012102412.asc');
elseif (Ncase == 2 )
  Zread1 = load('2012102600.asc');
  Zread2 = load('2012102700.asc');
elseif (Ncase == 3 )
  Zread1 = load('2012102700.asc');
  Zread2 = load('2012102800.asc');
elseif (Ncase == 4 )
  Zread1 = load('2012102800.asc');
  Zread2 = load('2012102900.asc');
end

Z1 = reshape(Zread1,Nlon,Nlat);
Z2 = reshape(Zread2,Nlon,Nlat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Set up the lat/long grid

lon     = linspace(0,357.5,Nlon);
latcomp = linspace(0,180,Nlat);
lat      = 90 - latcomp;

[LON,LAT] = meshgrid(lon,lat);
LON = LON'; LAT = LAT';  %  Matlab!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Plot the field

contours = 4500:50:6000;
subplot(1,2,1)
contourf(LON,LAT,Z1,contours)
colorbar('horiz')
subplot(1,2,2)
contourf(LON,LAT,Z2,contours)
colorbar('horiz')
pause

rx=73:144; ry=1:37;
subplot(1,2,1)
contourf(LON(rx,ry),LAT(rx,ry),Z1(rx,ry),contours)
colorbar('horiz')
subplot(1,2,2)
contourf(LON(rx,ry),LAT(rx,ry),Z2(rx,ry),contours)
colorbar('horiz')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Switch range to [-180,+180] lon.

least  = lon(1:Nlon/2);
lwest  = lon(Nlon/2+1:Nlon)-360;
lon = [ lwest least ];

Least  = LON(1:Nlon/2,:);
Lwest  = LON(Nlon/2+1:Nlon,:)-360;
LON = [ Lwest; Least ];

Zeast = Z1(1:Nlon/2,:); Zwest = Z1(Nlon/2+1:Nlon,:);
Z1 = [ Zwest; Zeast ];
Zeast = Z2(1:Nlon/2,:); Zwest = Z2(Nlon/2+1:Nlon,:);
Z2 = [ Zwest; Zeast ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Extend the arrays Z1 and Z2 by 
%%  repeating the values on the date-line
%%  (this facilitates interpolation)

lon = [ lon(1)-Dlon  , lon , lon(Nlon)+Dlon   ];
LON = [ LON(1,:)-Dlon; LON ; LON(Nlon,:)+Dlon ];
LAT = [ LAT(Nlon,:)  ; LAT ; LAT(1,:)         ];
Z1  = [ Z1(Nlon,:)   ;  Z1 ;   Z1(1,:)        ];
Z2  = [ Z2(Nlon,:)   ;  Z2 ;   Z2(1,:)        ];
Nlon = Nlon + 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Plot the switched field

subplot(1,2,1)
contourf(LON,LAT,Z1,contours)
colorbar('horiz')
subplot(1,2,2)
contourf(LON,LAT,Z2,contours)
colorbar('horiz')
pause

rx=1:73; ry=1:37;
subplot(1,2,1)
contourf(LON(rx,ry),LAT(rx,ry),Z1(rx,ry),contours)
colorbar('horiz')
subplot(1,2,2)
contourf(LON(rx,ry),LAT(rx,ry),Z2(rx,ry),contours)
colorbar('horiz')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Transform to PS grid.
%
%  Details of polar stereographic grid
%     Nx = 19 points in x-direction
%     Ny = 16 points in y-direction
%     Ds = 736 km at North Pole.
%     Xp: x-coordinate of North Pole (varies)
%     Yp: y-coordinate of North Pole (varies)
%     Cang: Angle between Date-line
%           and positive y axis (varies).
%     Origin of grid is at North Pole.
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 2E7/pi;
r2d = 180/pi; d2r = pi/180;

Ix = 1:Nx;
Iy = 1:Ny;
for nny=Iy
for nnx=Ix
  x = (nnx-Xp)*Ds;
  y = (nny-Yp)*Ds;
  r = sqrt(x^2+y^2);
  theta = atan2(y,x);
  lambda = theta + d2r*(90+cenang);
  if(lambda>pi) lambda = lambda - 2*pi; end
  phi = 2*((pi/4)-atan(r/(2*a)));
  lamd = r2d*lambda;
  phid = r2d*phi;
  LAMXY(nnx,nny) = lamd;
  PHIXY(nnx,nny) = phid;
end
end

Z1ps = zeros(Nx,Ny);  %  Heights on PS grid.
Z2ps = zeros(Nx,Ny);  %  Heights on PS grid.
for nny=1:Ny
for nnx=1:Nx
  lonxy = LAMXY(nnx,nny);
  latxy = PHIXY(nnx,nny);

  for nl=1:Nlon-1
    if(lon(nl)<=lonxy & lon(nl+1)>lonxy )
      nlft = nl;
    end
  end
  for nl=1:Nlat-1
  if(lat(nl)>=latxy & lat(nl+1)<latxy )
      nlow = nl;
    end
  end
  
  mulon = (lonxy-lon(nlft))/(lon(nlft+1)-lon(nlft));
  mulat = (latxy-lat(nlow))/(lat(nlow+1)-lat(nlow));

    %%   Check that mulon and mulat are in range.
    if(abs(mulon)>1)
      fprintf('nlft nlow mulon %g %g %g \n', ...
               nlft,nlow,mulon);
      fprintf('%g %g %g \n',lonxy,lon(nlft),lon(nlft+1))
      return
    end
    if(abs(mulat)>1)
      fprintf('nlft nlow mulat %g %g %g \n', ...
               nlft,nlow,mulat); 
      return
    end
  
  %%  Interpolate linearly
  Z1SW = Z1(nlft  ,nlow  ); Z2SW = Z2(nlft  ,nlow  );
  Z1SE = Z1(nlft+1,nlow  ); Z2SE = Z2(nlft+1,nlow  );
  Z1NW = Z1(nlft  ,nlow+1); Z2NW = Z2(nlft  ,nlow+1);
  Z1NE = Z1(nlft+1,nlow+1); Z2NE = Z2(nlft+1,nlow+1);
  
  Z1SO = (1-mulon)*Z1SW + mulon*Z1SE;
  Z1NO = (1-mulon)*Z1NW + mulon*Z1NE;
  Z1xy = (1-mulat)*Z1SO + mulat*Z1NO;
  Z1ps(nnx,nny) = Z1xy;
  
  Z2SO = (1-mulon)*Z2SW + mulon*Z2SE;
  Z2NO = (1-mulon)*Z2NW + mulon*Z2NE;
  Z2xy = (1-mulat)*Z2SO + mulat*Z2NO;
  Z2ps(nnx,nny) = Z2xy;

  Xps(nnx,nny) = nnx;
  Yps(nnx,nny) = nny;

end
end

%%  Plot the fields on the PS grid
%%  with the lat/long axes
conlon = linspace(-180,+180,37);
conlat = linspace(0,90,10);

subplot(1,2,1); hold on
contour(Xps,Yps,Z1ps,contours)
contour(Xps,Yps,LAMXY,conlon)
contour(Xps,Yps,PHIXY,conlat)

subplot(1,2,2); hold on
contour(Xps,Yps,Z2ps,contours)
contour(Xps,Yps,LAMXY,conlon)
contour(Xps,Yps,PHIXY,conlat)

hold off; pause; clf

%%  Plot inner region
rx = 3:Nx-2; ry = 2:Ny-2;
subplot(1,2,1); hold on
contour(Xps(rx,ry),Yps(rx,ry),Z1ps(rx,ry),contours);
contour(Xps(rx,ry),Yps(rx,ry),LAMXY(rx,ry),conlon);
contour(Xps(rx,ry),Yps(rx,ry),PHIXY(rx,ry),conlat);
colorbar('horiz')
subplot(1,2,2); hold on
contour(Xps(rx,ry),Yps(rx,ry),Z2ps(rx,ry),contours);
contour(Xps(rx,ry),Yps(rx,ry),LAMXY(rx,ry),conlon);
contour(Xps(rx,ry),Yps(rx,ry),PHIXY(rx,ry),conlat);
colorbar('horiz')
pause; hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Write out the analyses on the PS grid
if      (Ncase == 1 ) 
  save Case1-2012102200.z00 -ascii Z1ps;
  save Case1-2012102412.z00 -ascii Z2ps;
elseif (Ncase == 2 )
  save Case2-2012102600.z00 -ascii Z1ps;
  save Case2-2012102700.z00 -ascii Z2ps;
elseif (Ncase == 3 )
  save Case3-2012102700.z00 -ascii Z1ps;
  save Case3-2012102800.z00 -ascii Z2ps;
elseif (Ncase == 4 )
  save Case4-2012102800.z00 -ascii Z1ps;
  save Case4-2012102900.z00 -ascii Z2ps;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
