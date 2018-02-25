%% Assignment 2
% Question 1(a)
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

 
