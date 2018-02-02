close all
clear all

L=200e-9;
W=100e-9;
n=10000; %change
nsteps =1000; %change
num = n;

ang=randn(1,n)*2*pi;

m0=9.109e-31;
mn=0.26*m0;
T=300;
k=physconst('Boltzman');
tau=0.2e-9; %sec
 
vth = sqrt(2*k*T/mn)

x=rand(1,n)*L;
y=rand(1,n)*W;
xp=x;
yp=y;

vx=vth*ones(1,n).*cos(ang);
vy=vth*ones(1,n).*sin(ang);

dt=(L/vth)/100;

col=hsv(10); %vector of colours for particle trajectories

f1 = figure;
f2 = figure;
set(0, 'CurrentFigure', f2)
  for p = 1:10
        plot([x(p); xp(p)],[y(p); yp(p)],'color',col(p,:));  hold on
    end
    xlim([0 L])
    ylim([0 W])  
    
 f3 = figure;   
for i=1:nsteps

    xp=x;
    yp=y;

    dx=vx*dt;
    dy=vy*dt;

    x=x+dx;
    y=y+dy;

    xpath=abs(x-xp); %xpath calc before bounday adjustment
  
    %periodic boundaries for walls
     for a=1:n
        if (x(a)>L)
            x(a)=x(a)-L;
            xp(a)=xp(a)-L;
        elseif x(a)<0
            x(a) = x(a)+L;
            xp(a) = xp(a)+L;
        end
    %specular boundaries for ceiling and floor
        if y(a)>=W
            vy(a) = -vy(a);
         elseif y(a)<=0
            vy(a) = -vy(a);
        end
     end

    ypath=abs(y-yp); %ypath calc after boundary adjustment
    path = sqrt(xpath.*xpath + ypath.*ypath);
    
    velx = mean(abs(vx));
    vely = mean(abs(vy));
    v_inst=sqrt(velx*velx+vely*vely);
    
    Temp= v_inst*v_inst*mn/k ;

    set(0, 'CurrentFigure', f1)
    %plot trajectories 
    plot(x,y,'*');  
    xlim([0 L])
    ylim([0 W])
  title(sprintf('Monte Carlo Electron Simulation - Number of Electrons = %d', num))
  
    
    set(0, 'CurrentFigure', f2)
        %plot trajectories
    for p = 1:10
        plot([x(p); xp(p)],[y(p); yp(p)],'color',col(p,:));  hold on
    end
    title('Electron Simulation with Trajectories')
    xlim([0 L])
    ylim([0 W])
   
    
    set(0, 'CurrentFigure', f3)
    plot(i,Temp, 'o')
    hold on
    title(sprintf('Semiconductor Temperature = %s', Temp))
    
    pause(0.01);
    
end

 
meanfreepath = mean(path);%??


