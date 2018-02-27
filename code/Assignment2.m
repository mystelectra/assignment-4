%% Assignment 2
%% Question 1(a)
% The Finite Difference Method is used to solve for electrostatic potential
% of a region, with $V=V_0$ on the left side and $V=0$ on the right side of
% the region. The top and bottom of the region are not fixed, so the
% problem is essential in one dimension.

close all
clear all

nx=150; %L
ny=100; %W
BC_left=1;
BC_right=0;
BC_top=0;
BC_bottom=0;

G=sparse(nx*ny);
B=zeros(1,nx*ny);

for a=1:nx
    for b=1:ny
        n=a+(b-1)*nx;
        
        if a==1 
            %Left Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_left;
        elseif a==nx
            %Right Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_right;
        elseif b==1   
            % Top
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nym= a +(b)*nx;
            
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            B(n)=BC_top;
        elseif b==ny
            %Bottom
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
   
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1;
            B(n)=BC_bottom;
        else
            %All Central Nodes
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
            nym= a +(b)*nx;
            
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
    end
end

% figure(1)
% grid on
% spy(G)
% title("G-matrix composition")

V=G\B';

Vmap = zeros(nx,ny);
for a=1:nx
    for b = 1:ny
      n=a+(b-1)*nx;
      Vmap(a,b) = V(n);
    end
end

figure(2)
plot(Vmap)
xlabel(sprintf('Length = %d', nx))
ylabel(sprintf('Width = %d', ny)) 
title(sprintf('Electrostatic Potential: V= %d  at X=0, V=%d at X=%d',BC_left, BC_right,nx))

 
 
%% Question 1(b)
% The Finite Difference Method is used to solve for electrostatic potential
% of a region, with $V=V_0$ on the left side and $V=0$ on the right side of
% the region. The top and bottom of the region are not fixed, so the
% problem is essential in one dimension.

BC_left=1;
BC_right=1;
BC_top=0;
BC_bottom=0;

G=sparse(nx*ny);
B=zeros(1,nx*ny);

for a=1:nx
    for b=1:ny
        n=a+(b-1)*nx;
        
        if a==1 
            %Left Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_left;
        elseif a==nx
            %Right Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_right;
        elseif b==1   
            % Top
            G(:,n) = 0; 
            G(n,n) = 1;
            B(n)=BC_top;
        elseif b==ny
            %Bottom
            G(:,n) = 0; 
            G(n,n) = 1;
            B(n)=BC_bottom;
        else
            %All Central Nodes
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
            nym= a +(b)*nx;
            
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
    end
end

V=G\B';

Vmap = zeros(nx,ny);
for a=1:nx
    for b = 1:ny
      n=a+(b-1)*nx;
      Vmap(a,b) = V(n);
    end
end

figure(3)
surf(Vmap)
ylabel(sprintf('Length = %d', nx))
xlabel(sprintf('Width = %d', ny)) 
zlabel('Potential')
title(sprintf('Electrostatic Potential: V= %d  at X=0, V=%d at X=%d',BC_left, BC_right,nx))

%% Analytical Solution to electrostatic potential where $V_0 = 1$
L=nx;
W=ny;
Vx = zeros(L,W);
for n = 1:2:100
    for x = 1:L
        for y = 1:W
            Vx(x,y)=Vx(x,y)+(4/pi)*(1/n)*cosh((n*pi*(x-L/2))/W)/cosh(n*pi*L/W).*sin(n*pi*y/W);
        end
    end
   figure(4)
   mesh(Vx)
   hold on
   pause(0.01)
end
    
%% 
% Ultimately, the solution should be stopped when the solution reaches the
% desired accuracy, indicated by some pre-decided variance factor. In this
% example, I arbitrarily chose a finite number of repetitions to simulate
% the solution. 100 repetitions is more than sufficient to obtain the correct
% representation. More repetitions than this does not provide increased
% accuracy. The analytical and mesh solutions agree to the extenet that
% supports the use of the mesh solution. 

%% Question 2
% In this section we invesitgate the Electric field and current diensity of
% a region under the same conditions as above, but this time we will
% introduce a dielectric type material, with a low conductivity and see how
% it changes the potential of the region.

nx=150; %L
ny=100; %W
BC_left=1;
BC_right=0;
BC_top=0;
BC_bottom=0;

G=sparse(nx*ny);
B=zeros(1,nx*ny);

Llow = round(L/3);
Lhigh= round(L-L/3);
Wlow =round(W/3);
Whigh=round(W-W/3);

sig = ones(W,L);

% set conductivity in boxes
for a = 1:Wlow
    for b=Llow:Lhigh
        sig(a,b) = 0.01;
    end
end
for a = Whigh:W-1
    for b=Llow:Lhigh
        sig(a,b) = 0.01;
    end
end
figure(5)
mesh(sig)


for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
        
        if a==1 
            %Left Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_left;
        elseif a==ny
            %Right Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_right;
        elseif b==1   
            % Top
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nym= a +(b)*nx;
          
            ym = (sig(a,b)+sig(a,b+1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            B(n)=BC_top;
        elseif b==nx
            %Bottom
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
   
            yp = (sig(a,b)+sig(a,b-1))/2;
            xp = (sig(a-1,b)+sig(a,b))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nyp) = yp;
            B(n)=BC_bottom;
        else
            %All Central Nodes
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
            nym= a +(b)*nx;
            
            ym = (sig(a,b)+sig(a,b+1))/2;
            yp = (sig(a,b)+sig(a,b-1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            G(n,nyp) = yp;
            B(n)=0;
        end
    end
end


V=G\B';

Vmap = zeros(ny,nx);
for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
      Vmap(a,b) = V(n);
    end
end

figure(6)
mesh(Vmap)
    xlim([0 L])
    ylim([0 W])
xlabel(sprintf('Length = %d', nx))
ylabel(sprintf('Width = %d', ny)) 
title(sprintf('Electrostatic Potential with Conductivity: V= %d  at X=0, V=%d at X=%d',BC_left, BC_right,nx))

[Ex,Ey] = gradient(-Vmap);

figure(7)
quiver(Ex,Ey)
    xlim([0 L])
    ylim([0 W])
title('Electric Field')

Jx=sig.*Ex;
Jy=sig.*Ey;

figure(8)
quiver(Jx,Jy)
    xlim([0 L])
    ylim([0 W])
title('Current density')


%% Mesh Density
% In this section, we investigate how the mesh density changes the
% resulting current

nx=75; %L
for nx = 20:10:150
ny=round((2/3)*nx); %W
L=nx;
W=ny;
BC_left=1;
BC_right=0;
BC_top=0;
BC_bottom=0;

G=sparse(nx*ny);
B=zeros(1,nx*ny);

Llow = round(L/3);
Lhigh= round(L-L/3);
Wlow =round(W/3);
Whigh=round(W-W/3);

sig = ones(W,L);

% set conductivity in boxes
for a = 1:Wlow
    for b=Llow:Lhigh
        sig(a,b) = 0.01;
    end
end
for a = Whigh:W-1
    for b=Llow:Lhigh
        sig(a,b) = 0.01;
    end
end

for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
        
        if a==1 
            %Left Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_left;
        elseif a==ny
            %Right Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_right;
        elseif b==1   
            % Top
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nym= a +(b)*nx;
          
            ym = (sig(a,b)+sig(a,b+1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            B(n)=BC_top;
        elseif b==nx
            %Bottom
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
   
            yp = (sig(a,b)+sig(a,b-1))/2;
            xp = (sig(a-1,b)+sig(a,b))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nyp) = yp;
            B(n)=BC_bottom;
        else
            %All Central Nodes
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
            nym= a +(b)*nx;
            
            ym = (sig(a,b)+sig(a,b+1))/2;
            yp = (sig(a,b)+sig(a,b-1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            G(n,nyp) = yp;
            B(n)=0;
        end
    end
end


V=G\B';

Vmap = zeros(ny,nx);
for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
      Vmap(a,b) = V(n);
    end
end

[Ex,Ey] = gradient(-Vmap);

Jx=sig.*Ex;
Jy=sig.*Ey;


I=0;
for b=1:ny
    I=I+Jy(1,b)*(ny);%
end

figure(9)
plot(nx,I,'*')
hold on
ylabel('Current')
xlabel('Number of Mesh Divisions')
title('Current vs. Mesh Size')
pause(0.01);

end

%% Bottle-Neck Size
% Here we determine the effects of the size of the bottle-neck on the
% analysis.

nx=75; %L
ny=round((2/3)*nx); %W
L=nx;
W=ny;
BC_left=1;
BC_right=0;
BC_top=0;
BC_bottom=0;

G=sparse(nx*ny);
B=zeros(1,nx*ny);

for div=3:0.1:4
Llow = round(L/div);
Lhigh= round(L-L/div);
Wlow =round(W/div);
Whigh=round(W-W/div);

sig = ones(W,L);

% set conductivity in boxes
for a = 1:Wlow
    for b=Llow:Lhigh
        sig(a,b) = 0.01;
    end
end
for a = Whigh:W-1
    for b=Llow:Lhigh
        sig(a,b) = 0.01;
    end
end

for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
        
        if a==1 
            %Left Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_left;
        elseif a==ny
            %Right Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_right;
        elseif b==1   
            % Top
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nym= a +(b)*nx;
          
            ym = (sig(a,b)+sig(a,b+1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            B(n)=BC_top;
        elseif b==nx
            %Bottom
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
   
            yp = (sig(a,b)+sig(a,b-1))/2;
            xp = (sig(a-1,b)+sig(a,b))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nyp) = yp;
            B(n)=BC_bottom;
        else
            %All Central Nodes
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
            nym= a +(b)*nx;
            
            ym = (sig(a,b)+sig(a,b+1))/2;
            yp = (sig(a,b)+sig(a,b-1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            G(n,nyp) = yp;
            B(n)=0;
        end
    end
end


V=G\B';

Vmap = zeros(ny,nx);
for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
      Vmap(a,b) = V(n);
    end
end

[Ex,Ey] = gradient(-Vmap);

Jx=sig.*Ex;
Jy=sig.*Ey;


I=0;
for b=1:ny
    I=I+Jy(1,b)*(ny);%
end

figure(10)
plot(div,I,'*')
hold on
ylabel('Current')
xlabel('Proportion of Region in Bottle-Neck')
title('Current vs. Bottle-Neck')
pause(0.01);

end

%% Varying Conductivity
% In this section, we investigate how varying conductivy of materials will
% affect the resulting current in the region.

nx=75; %L
ny=round((2/3)*nx); %W
L=nx;
W=ny;
BC_left=1;
BC_right=0;
BC_top=0;
BC_bottom=0;

G=sparse(nx*ny);
B=zeros(1,nx*ny);

div = 3;
Llow = round(L/div);
Lhigh= round(L-L/div);
Wlow =round(W/div);
Whigh=round(W-W/div);



for s=0.001:0.005:2
   % set conductivity in boxes
   for a = 1:Wlow
       for b=Llow:Lhigh
           sig(a,b) = s;
       end
   end
   for a = Whigh:W-1
       for b=Llow:Lhigh
           sig(a,b) = s;
       end
   end

for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
        
        if a==1 
            %Left Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_left;
        elseif a==ny
            %Right Side
            G(n,:) = 0; 
            G(n,n) = 1;
            B(n)=BC_right;
        elseif b==1   
            % Top
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nym= a +(b)*nx;
          
            ym = (sig(a,b)+sig(a,b+1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            B(n)=BC_top;
        elseif b==nx
            %Bottom
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
   
            yp = (sig(a,b)+sig(a,b-1))/2;
            xp = (sig(a-1,b)+sig(a,b))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nyp) = yp;
            B(n)=BC_bottom;
        else
            %All Central Nodes
            nxm= a-1 +(b-1)*nx;
            nxp= a+1 +(b-1)*nx;
            nyp= a +(b-2)*nx;
            nym= a +(b)*nx;
            
            ym = (sig(a,b)+sig(a,b+1))/2;
            yp = (sig(a,b)+sig(a,b-1))/2;
            xm = (sig(a-1,b)+sig(a,b))/2;
            xp = (sig(a+1,b)+sig(a,b))/2;
          
            G(n,n) =-(yp+ym+xp+xm);
            G(n,nxm) = xm;
            G(n,nxp) = xp;
            G(n,nym) = ym;
            G(n,nyp) = yp;
            B(n)=0;
        end
    end
end


V=G\B';

Vmap = zeros(ny,nx);
for a=1:ny
    for b=1:nx
        n=a+(b-1)*nx;
      Vmap(a,b) = V(n);
    end
end

[Ex,Ey] = gradient(-Vmap);

Jx=sig.*Ex;
Jy=sig.*Ey;


I=0;
for b=1:ny
    I=I+Jy(1,b)*(ny);%
end

figure(10)
plot(s,I,'*')
hold on
ylabel('Current')
xlabel('Conductivity')
title('Current vs. Conductivity')
pause(0.01);

end



