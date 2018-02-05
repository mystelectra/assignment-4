close all
clear all

L=200e-9;
W=100e-9;
n=10000; %change
nsteps =1000; %change

tau_mn=0.2e-12 %seconds

m0=9.109382e-31; %electron mass
mn=0.26*m0;
T=300; %Kelvin
k=physconst('Boltzman');

vth = sqrt(k*T/mn);
mfp = tau_mn*vth % meters

%inititalize particle locations
x=rand(1,n)*L;
y=rand(1,n)*W;

%initialize previous location as first location just to get the plot
%started
xp = x;
yp = y;

%initialize random velocities
vx=vth*randn(1,n)/sqrt(2);
vy=vth*randn(1,n)/sqrt(2);
speed=sqrt(vx.*vx +vy.*vy);

f1 = figure;
f2 = figure;
f3=figure;

set(0, 'CurrentFigure', f1)
subplot(2,2,1); 
histogram (vx);
title ('Distribution-Vx')
subplot(2,2,2); 
histogram (vy);
title ('Distribution-vy')
subplot(2,2,[3,4]); 
histogram (speed);
title ('Distribution-speed'); 

dt=(L/vth)/100;
av_temp = zeros(1,nsteps);

%vector of colours for particle trajectories
col=hsv(10); 
set(0, 'CurrentFigure', f2)
  for p = 1:10
        plot([x(p); xp(p)],[y(p); yp(p)],'color',col(p,:));  hold on
    end
    xlim([0 L])
    ylim([0 W])
    
%main timeloop
for aa=1:nsteps
    xp=x;
    yp=y;

    %scattering
    pscat=1-exp(-dt/tau_mn);
    scatCount= 0;
    
    for bb=1:n
        if (pscat > rand())
            vx(bb)=vth*randn()/sqrt(2);
            vy(bb)=vth*randn()/sqrt(2);
            scatCount = scatCount+1;
        end
    end

    dx=vx*dt;
    dy=vy*dt;
    
    %increment every particle over dt
    x=x+dx;
    y=y+dy;

    %xpath calc before boundary adjustment
    xpath=abs(x-xp); 
    
    %travelling restrictions (WALL)
     for a=1:n        
        %periodic boundaries at x=0 and x=L
         if (xp(a)< L && x(a)>=L)
            x(a)=x(a)-L;
            xp(a)=xp(a)-L;
         elseif (xp(a)< 0 && x(a)<0)
            x(a) = x(a)+L;
            xp(a)=xp(a)+L;
         end
         
        %specular boundaries at y=0 and y=W
        if (y(a)>=W || y(a)<=0)
            vy(a) = -vy(a);
        elseif y(a)<=0
          vy(a) = -vy(a);
        end
     end %end travelling restrictions loop
     

     %ypath calc after boundary adjustment
    ypath=abs(y-yp); 
    
    %calculate path 
    path = sqrt(xpath.*xpath + ypath.*ypath);    
    
    
    
% plot(x,y,'o');hold on
  set(0, 'CurrentFigure', f2)
    %plot trajectories 
    for p = 1:10
        plot([x(p); xp(p)],[y(p); yp(p)],'color',col(p,:));  hold on
    end
    xlim([0 L])
    ylim([0 W])
    title ('Electron Collisions with Mean Free Path')
    pause(0.01);
    
    
    velx = mean(vx.^2);
    vely = mean(vy.^2);
    v_inst=sqrt(velx+vely);
    
    Temp= v_inst*v_inst*mn/k ;
   set(0, 'CurrentFigure', f3)
    plot(aa,Temp, 'o')
    hold on
    title(sprintf('Semiconductor Temperature = %s', Temp))
     
    av_temp(aa) = Temp;
end


AverageTemperature = mean(av_temp)
meanfreepath = mean(path)
meanfreetime = meanfreepath/(mean(v_inst))
