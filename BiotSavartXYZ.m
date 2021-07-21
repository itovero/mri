%-------------------------------------------------------------------------%
%-------Cicular Coil Magnetic Field Simulation--------------------------%
%-------------------------------------------------------------------------%

clc
close all
clear all

%-------------------------------------------------------------------------%
%---Coil is in the X-Y plane and Magnetic Field is Evaluated -------------%
%-------------at every point in the Y-Z plane(X=0)------------------------%
%-------------------------------------------------------------------------%
Nz=51;  % No. of grids in Z-axis
Ny=51;  % No. of grids in Y-axis
Nx=51;  % No. of grids in X-axis
N=25;   % No of grids in the coil (X-Y plane)
Ra=6;    % Radius of the coil in the X-Y plane
I=3;    % current in the coil
u0=1;   % for simplicity, u0 is taken as 1 (permitivity)
phi=-pi/2:2*pi/(N-1):3*pi/2; % For describing a circle (coil)

Xc=Ra*cos(phi); % X-coordinates of the coil
Yc=Ra*sin(phi); % Y-coordinates of the coil

xp =-25:1:25; % X-coordinates of the plane 
yp =-25:1:25; % Y-coordinates of the plane 
zp =-25:1:25;% Z-coordinates of the plane 

X(1:Nx, 1:Ny,1:Nz)=0;
Y(1:Nx, 1:Ny,1:Nz)=0; % This array is for 1-d to 2-d conversion of coordinates
Z(1:Nx, 1:Ny,1:Nz)=0;

for i=1:Nx
    X(i, :,:)=xp(i); % all y-coordinates value in 2-d form
end
for i=1:Ny
    Y(:, i,:)=yp(i); % all y-coordinates value in 2-d form
end
for i=1:Nz
    Z(:, :,i)=zp(i);% all z-coordinates value in 2-d form
end

%-------------------------------------------------------------------------
%-variable "a" in for loop is along Y and variable "b" is along Z-axis----
%-------------------------------------------------------------------------


for a=1:Nx  
for b=1:Ny
for c=1:Nz
    
%-------------------------------------------------------------------------    
% calculate R-vector from the coil(X-Y plane)to Y-Z plane where we are 
% interested to find the magnetic field and also the dl-vector along the
% coil where current is flowing
% (note that R is the position vector pointing from coil (X-Y plane) to
% the point where we need the magnetic field (in this case Y-Z plane).)
% dl is the current element vector which will make up the coil------------
%-------------------------------------------------------------------------

for i=1:N-1
Rx(i)=-0.5*(Xc(i)+Xc(i+1));
Ry(i)=(Y(a,b,c)-(0.5*(Yc(i)+Yc(i+1))));
Rz(i)=Z(a,b,c);
dlx(i)=Xc(i+1)-Xc(i);
dly(i)=Yc(i+1)-Yc(i);
end
Rx(N)=-0.5*(Xc(N)+Xc(1));
Ry(N)=(Y(a,b,c)-(0.5*(Yc(N)+Yc(1))));
Rz(N)=Z(a,b,c);
dlx(N)=-Xc(N)+Xc(1);
dly(N)=-Yc(N)+Yc(1);

%--------------------------------------------------------------------------
% for all the elements along coil, calculate dl cross R -------------------
% dl cross R is the curl of vector dl and R--------------------------------
% XCross is X-component of the curl of dl and R, similarly I get Y and Z- 
%--------------------------------------------------------------------------

for i=1:N
Xcross(i)=dly(i).*Rz(i);
Ycross(i)=-dlx(i).*Rz(i);
Zcross(i)=(dlx(i).*Ry(i))-(dly(i).*Rx(i));
R(i)=sqrt(Rx(i).^2+Ry(i).^2+Rz(i).^2);
end

%-------------------------------------------------------------------------
% this will be the biot savarts law equation------------------------------
%--------------------------------------------------------------------------

Bx1=(I*u0./(4*pi*(R.^3))).*Xcross;
By1=(I*u0./(4*pi*(R.^3))).*Ycross;
Bz1=(I*u0./(4*pi*(R.^3))).*Zcross;
%--------------------------------------------------------------------------
% now we have  magnetic field from all current elements in the form of an
% array named Bx1,By1,Bz1, now its time to add them up to get total
% magnetic field 
%-------------------------------------------------------------------------

BX(a,b,c)=0;       % Initialize sum magnetic field to be zero first
BY(a,b,c)=0;
BZ(a,b,c)=0;

%--------------------------------------------------------------------------
% here we add all magnetic field from different current elements which 
% make up the coil
%--------------------------------------------------------------------------

for i=1:N   % loop over all current elements along coil    
    BX(a,b,c)=BX(a,b,c)+Bx1(i);
    BY(a,b,c)=BY(a,b,c)+By1(i);
    BZ(a,b,c)=BZ(a,b,c)+Bz1(i);
end

%-------------------------------------------------------------------------

end
end
end

%------------------------------------------------------------------------
%---BX is a null field and only BY and BZ are color mapped...............
%--------------------------------------------------------------------------

figure(1)
plot(Xc,Yc,'linewidth',3)
axis([-20 20 -20 20])
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('square loop co-ordinates','fontsize',14)
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(1); 
set(fh, 'color', 'white'); 
grid on

figure(2)
fig2_BZ = squeeze(BZ(1, :, :));
lim1=min(min(fig2_BZ));
lim2=max(max(fig2_BZ));
steps=(lim2-lim1)/100;
contour(zp,yp,fig2_BZ,lim1:steps:lim2)
axis([-25 25 -25 25])
xlabel('Z-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('BZ component','fontsize',14)
colorbar('location','eastoutside','fontsize',14);
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(2); 
set(fh, 'color', 'white'); 

figure(3)
fig3_BZ = squeeze(BZ(:, :, 51));
lim1=min(min(fig3_BZ));
lim2=max(max(fig3_BZ));
steps=(lim2-lim1)/100;
contour(xp,yp,fig3_BZ,lim1:steps:lim2)
axis([-25 25 -25 25])
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('BZ component','fontsize',14)
colorbar('location','eastoutside','fontsize',14);
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(2); 
set(fh, 'color', 'white'); 
%figure(3)
%lim1=min(min(BY));
%lim2=max(max(BY));
%steps=(lim2-lim1)/100;
%contour(zp,yp,BY,lim1:steps:lim2)
%axis([-25 25 -25 25])
%xlabel('Z-axis','fontsize',14)
%ylabel('Y-axis','fontsize',14)
%title('BY component','fontsize',14)
%colorbar('location','eastoutside','fontsize',14);
%h=gca; 
%get(h,'FontSize') 
%set(h,'FontSize',14)
%h = get(gca, 'ylabel');
%fh = figure(3); 
%set(fh, 'color', 'white'); 


%figure(4)
%quiver(zp,yp,BZ,BY,2)
%axis([-25 25 -25 25])
%xlabel('Z-axis','fontsize',14)
%ylabel('Y-axis','fontsize',14)
%title('B-field Vector flow','fontsize',14)
%h=gca; 
%get(h,'FontSize') 
%set(h,'FontSize',14)
%h = get(gca, 'ylabel');
%fh = figure(4); 
%set(fh, 'color', 'white'); 

%--------------------------------------------------------------------------