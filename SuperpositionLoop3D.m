clc
close all
clear all

numRows = 2;
numCols = 4;
numLoops = numRows * numCols;
for k = 1:numLoops
    x(k) = LoopClass;
    x(k).radius = 3;
    x(k).y_offset = floor((k - 1) / numCols) * 3;
    x(k).x_offset = mod((k - 1), numCols) * 3;
    x(k).I = 3; %Set all current to 3Amps
end

N=25;   % No of grids in the coil (X-Y plane)
phi=-pi/2:2*pi/(N-1):3*pi/2; % For describing a circle (coil)


%Plotting the Loop
figure(1)
axis([-20 20 -20 20])
for k = 1:numLoops
    x(k).xgraph=x(k).radius*cos(phi); % X-coordinates of the coil
    x(k).ygraph=x(k).radius*sin(phi); % Y-coordinates of the coil
    plot(x(k).xgraph + x(k).x_offset,x(k).ygraph + x(k).y_offset,'linewidth',3)
    hold on
end
hold off
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('Superimposed loops','fontsize',14)
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(1); 
set(fh, 'color', 'white'); 
grid on

Nx=51;
Nz=51;  % No. of grids in Z-axis
Ny=51;  % No. of grids in Y-axis
u0=1;   % for simplicity, u0 is taken as 1 (permitivity) 

BX(1:Nx,1:Ny,1:Nz)=0; % Initialize sum magnetic field to be zero first
BY(1:Nx,1:Ny,1:Nz)=0;
BZ(1:Nx,1:Ny,1:Nz)=0;
 
 %Initializing the Magnetic Field for all loops
 for k = 1:numLoops
     disp(k);
     
     X(1:Nx, 1:Ny,1:Nz)=0;
     Y(1:Nx, 1:Ny,1:Nz)=0; % This array is for 1-d to 2-d conversion of coordinates
     Z(1:Nx, 1:Ny,1:Nz)=0;
     
     xp =-25 + x(k).x_offset:1:25 + x(k).x_offset;
     yp =-25 + x(k).y_offset:1:25 + x(k).y_offset; % Y-coordinates of the plane 
     zp =-25:1:25;% Z-coordinates of the plane 
     
     for i=1:Nx
         X(i, :,:)=xp(i); % all y-coordinates value in 2-d form
     end
     for i=1:Ny
         Y(:, i,:)=yp(i); % all y-coordinates value in 2-d form
     end
     for i=1:Nz
         Z(:, :,i)=zp(i);% all z-coordinates value in 2-d form
     end
    
     for a = 1:Nx
        for b = 1:Ny
            for c = 1:Nz
                for i = 1:N - 1
                 x(k).Rx(i)=(X(a,b,c)-(0.5*(x(k).xgraph(i)+x(k).xgraph(i+1))));
                 x(k).Ry(i)=(Y(a,b,c)-(0.5*(x(k).ygraph(i)+x(k).ygraph(i+1))));
                 x(k).Rz(i)=Z(a,b,c);
                 x(k).dlx(i)=x(k).xgraph(i+1)-x(k).xgraph(i);
                 x(k).dly(i)=x(k).ygraph(i+1)-x(k).ygraph(i);
                end
            x(k).Rx(N)=(X(a,b,c)-(0.5*(x(k).xgraph(N)+x(k).xgraph(1))));
            x(k).Ry(N)=(Y(a,b,c)-(0.5*(x(k).ygraph(N)+x(k).ygraph(1))));
            x(k).Rz(N)=Z(a,b);
            x(k).dlx(N)=-x(k).xgraph(N)+x(k).xgraph(1);
            x(k).dly(N)=-x(k).ygraph(N)+x(k).ygraph(1);
            
            for i=1:N
            x(k).Xcross(i)=x(k).dly(i).*x(k).Rz(i);
            x(k).Ycross(i)=-x(k).dlx(i).*x(k).Rz(i);
            x(k).Zcross(i)=(x(k).dlx(i).*x(k).Ry(i))-(x(k).dly(i).*x(k).Rx(i));
            x(k).R(i)=sqrt(x(k).Rx(i).^2+x(k).Ry(i).^2+x(k).Rz(i).^2);
            end
            
            x(k).Bx1=(x(k).I*u0./(4*pi*(x(k).R.^3))).*x(k).Xcross;
            x(k).By1=(x(k).I*u0./(4*pi*(x(k).R.^3))).*x(k).Ycross;
            x(k).Bz1=(x(k).I*u0./(4*pi*(x(k).R.^3))).*x(k).Zcross;
            
            %BX(a,b,c)=0;       % Initialize sum magnetic field to be zero first
            %BY(a,b,c)=0;
            %BZ(a,b,c)=0;
            
            for i=1:N   % loop over all current elements along coil    
                BX(a,b,c)=BX(a,b,c)+x(k).Bx1(i);
                BY(a,b,c)=BY(a,b,c)+x(k).By1(i);
                BZ(a,b,c)=BZ(a,b,c)+x(k).Bz1(i);
            end
            end
        end
     end
 end 
%Plotting BZ Component of Magnetic Field
figure(2)
for k = 1:numLoops
    fig2_BZ = squeeze(BZ(1, :, :));
    lim1=min(min(BZ));
    lim2=max(max(BZ));
    steps=(lim2-lim1)/100;
    %contour(zp,yp,x(k).BZ,lim1:steps:lim2)
    imagesc(zp, yp, fig2_BZ)
    axis([-25 25 -25 25])
    hold on
end
hold off
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
for k = 1:numLoops
    fig3_BZ = squeeze(BZ(1, :, :));
    lim1=min(min(BY));
    lim2=max(max(BY));
    steps=(lim2-lim1)/100;
    %contour(zp,yp,x(k).BY,lim1:steps:lim2)
    imagesc(xp, yp, fig3_BZ)
    axis([-25 25 -25 25])
    hold on
end
hold off
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('BZ component','fontsize',14)
colorbar('location','eastoutside','fontsize',14);
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(3); 
set(fh, 'color', 'white'); 


%figure(4)
%for k = 1:numLoops
%    quiver(zp,yp,x(k).BZ,x(k).BY,2)
%    axis([-25 25 -25 25])
%    hold on
%end
%hold off
%xlabel('Z-axis','fontsize',14)
%ylabel('Y-axis','fontsize',14)
%title('B-field Vector flow','fontsize',14)
%h=gca; 
%get(h,'FontSize') 
%set(h,'FontSize',14)
%h = get(gca, 'ylabel');
%fh = figure(4); 
%set(fh, 'color', 'white'); 

%Plotting BZ Component of Magnetic Field
%figure(5)
%for k = 1:numLoops
%    lim1=min(min(x(k).BZ));
%    lim2=max(max(x(k).BZ));
%    steps=(lim2-lim1)/100;
%    contour(xp,yp,abs(x(k).BZ),lim1:steps:lim2)
%    axis([-25 25 -25 25])
%    hold on
%end
%hold off
%xlabel('X-axis','fontsize',14)
%ylabel('Y-axis','fontsize',14)
%title('BZ component','fontsize',14)
%colorbar('location','eastoutside','fontsize',14);
%h=gca; 
%get(h,'FontSize') 
%set(h,'FontSize',14)
%h = get(gca, 'ylabel');
%fh = figure(2);
%set(fh, 'color', 'white'); 
 