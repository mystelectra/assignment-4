%% Assignment 4 - Circuit modelling
%% Question 1
% The differential equations for the ciruit in the time-domain are as follows:
%
% 1) $V1 = Vin$ 
%
% 2) $V2 - L(d/dt)*I3 - V3 = 0$ 
%
% 3) $V3 - I3*R3=0$ 
%
% 4) $I3 - V1/R1 - V2*(R2+R1)/(R1*R2) -C(V1-V2) = 0$
%
% 5) $V4 - a*I3 = 0$
%
% 6) $V5*(R0+R4) -R0*V4 = 0$
%
% The same equations in the frequency domain are:
%
% 1) $V1 = Vin$
%
% 2) $V2 - j \omega L*I3  - V3 = 0$
%
% 3) $V3 - I3*R3=0$
%
% 4) $I3 - V1/R1 - V2*(R2+R1)/(R1*R2) - j \omega C*(V1-V2) = 0$
%
% 5) $V4 - a*I3 = 0$
%
% 6) $V5*(R0+R4) -R0*V4 = 0$
%

close all
clear all

R1=1;
c=0.25;
R2=2;
L=0.2;
R3=10;
a=100;
R4=0.1;
R0=1000;
Vin=1;

%V=[V1 V2 V3 I3 V4 V5]

F = [Vin 0 0 0 0 0];
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 0 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];
G=[1 0 0 0 0 0; 
    0 1 -1 0 0 0; 
    0 0 1 -R3 0 0;
    -1/R1 ((R1+R2)/(R1*R2)) 0 1 0 0; 
    0 0 0 -a 1 0;
    0 0 0 0 -R0 (R0+R4)];

% DC case
out= zeros(21,6);
in= zeros(1,21);
b=0;
for vin= -10:1:10
    b=b+1;
    F(1)=vin;
    V=G\F';
    
    out(b,:)=V;
    in(b) = vin;
end

figure(1)
plot(in,out(:,6))%plot vin vs. V5 (Vout)
title('Vout')
xlabel('Vin sweep from -10V to 10V')
ylabel('Vin')

figure(2)
plot(in,out(:,3))
title('V3')
xlabel('Vin sweep from -10V to 10V')
ylabel('Vin')

%AC case - frequency sweep
out= zeros(1001,6);
in= zeros(1,1001);
b=0;
for w= 1:1:1000
    b=b+1;
  
    V=(G+1j*w*C)\F';
    
    out(b,:)=V;
    in(b) = w;
end

figure(3)
semilogx(in,real(out(:,6)))
title('Vout')
xlabel('frequency sweep from 1Hz to 1kHz')


figure(4)
semilogx(in,real(20*log10(real(out(:,6)/out(:,1)))))
title('gain over frequency (BODE)')
xlabel('frequency sweep from 1Hz to 1kHz')

%AC case - capacitance sweep (cs)
n=1000;
out= zeros(n,6);
gain= zeros(1,n);
b=0;
w=pi;
 cs=c+0.05*randn(1,n); %std=0.05; c=0.25;
for a=1:n
   c=cs(a);
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 0 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];
    
    V=(G+1j*w*C)\F';
     out(a,:)=V;
    gain(a) = (abs(V(6))/abs(V(1)));
  
end


figure(5)
histogram(gain,100)
title('gain as a function of capacitance')
xlabel('gain')

%% Question 2- Transient Circuit Simulation
% Using the above formulation, we can do a transient analysis of the
% circuit. 
% a) This circuit functions as a low pass filter. 
%
% b) For a low pass filter, signals at frequencies below some cut-off frequency 
% are sent to the output whereas signals at frequencies above the cut-off 
% frequency are attenuated.
%
% The finite difference equation used for the numerical solution of
% transient analysis is $V_i = (F + C*V_{(i-1)}/dt)(C/dt +G)^{(-1)}$


% AC case - transient sweep
% Step
n=1000;
out= zeros(n,6);
in= zeros(1,n);
Vin=0;

F = [Vin 0 0 0 0 0];
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 0 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];

dt = 1e-3;

Vt=zeros(n,6); %present voltage
Vt(1,:) = G\F' ; %V(t=0)

% C = C*0;
H = (C/dt + G);
Hi = inv(H);

t(1) = 0;

for k= 2:1:n %ms
    t(k) = (k-1)*1e-3;
    
    if k<30
        Vin=0;
        F = [Vin 0 0 0 0 0];
        Vt(k,:) = Hi*(F' + C/dt*Vt(k-1,:)');
    elseif k>=30
        Vin=1;
        F = [Vin 0 0 0 0 0];
        Vt(k,:) = Hi*(F' + C/dt*Vt(k-1,:)');
    end
    out(k,:)=Vt(k,:);
    %plot(k,out(k,1),'+');hold on; plot(k,out(k,6),'*');
    in(k) = k; %ms
end

Fs=1./t;

figure(6)
% Vin and Vout for step
plot(in,out(:,6))
hold on;
plot(in,out(:,1))
title('Vout and Vin as a Step')
xlabel('Time(ms)')
hold off
legend('Vout','Vin')

Fin = fft(out(:,1));
Fout = fft(out(:,6));


figure(7)
semilogy(abs(fftshift(Fin)))
hold on
 semilogy(abs(fftshift(Fout)))
title('Fourier Transform - Frequency Response of Step')
 xlabel('f (Hz)')
 ylabel('|P1(f)|')
legend('Vin','Vout')

%-----------------
%Sine Function
n=1000;
out= zeros(n,6);
in= zeros(1,n);

dt = 1e-3;

Vt=zeros(n,6); %present voltage
f=1/0.03; %Hz
t(1) = 0;

Vin=sin(2*pi*f*t(1));

F = [Vin 0 0 0 0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*(F' );

t(1) = 0;

for  k= 2:1:n %ms
    t(k) = (k-1)*1e-3;
    
    Vin=sin(2*pi*f*t(k));
    F = [Vin 0 0 0 0 0];
    
        Vt(k,:) = Hi*(F' + C/dt*Vt(k-1,:)');
   
    out(k,:)=Vt(k,:);
    in(k) = k; %ms
   
end
Fs=1./t;


figure(8)
plot(in,out(:,6))
hold on;
plot(in,out(:,1))
title('Vout and Vin as a Sine Function')
xlabel('Time(ms)')
hold off
legend('Vout','Vin')

Fin = fft(out(:,1));
Fout = fft(out(:,6));


figure(9)
semilogy(abs(fftshift(Fin)))
hold on
 semilogy(abs(fftshift(Fout)))
title('Frequency Response of Sine Function')
 xlabel('f (Hz)')
 ylabel('|P1(f)|')
legend('Vin','Vout')


%-----------------
% Gaussian Pulse
n=1000;
out= zeros(n,6);
in= zeros(1,n);

delay = 0.06;%s
sd = 0.03; %s
dt = 1e-3;

Vt=zeros(n,6); %present voltage
t(1) = 0;

Vin=exp(-((t(1)-delay)/sd) *((t(1)-delay)/sd)/2);

F = [Vin 0 0 0 0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*(F' );

for  k= 1:1:n %ms
    t(k) = (k-1)*1e-3;
    
    Vin=exp(-((t(k)-delay)/sd) *((t(k)-delay)/sd)/2);
    F = [Vin 0 0 0 0 0];
    if k>1
       Vt(k,:) = Hi*(F' + C/dt*Vt(k-1,:)');
    end
    out(k,:)=Vt(k,:);
    in(k) = k; %ms
   
end

Fs=1./t;

figure(10)
plot(t,out(:,6))
hold on;
plot(t,out(:,1))
title('Vout and Vin as a Gaussian Pulse')
xlabel('Time(ms)')
hold off
legend('Vout','Vin')

Fin = fft(out(:,1));
Fout = fft(out(:,6));


figure(11)
semilogy(abs(fftshift(Fin)))
hold on
semilogy(abs(fftshift(Fout)))
title('Frequency Response of Gaussian Pulse')
 xlabel('f (Hz)')
 ylabel('|P1(f)|')
legend('Vin','Vout')

%% Increased Time Step
% With the increased time step, we can see that the circuit changes from
% being an underdamped system to an overdamped system. The output is more
% sluggish and not as responsive to changes in the input.

n=100;
out= zeros(n,6);
in= zeros(1,n);
t= zeros(1,n);

delay = 0.06;%s
sd = 0.03; %s
dt = 1e-3;

Vt=zeros(n,6); %present voltage
t(1) = 0;

Vin=exp(-((t(1)-delay)/sd) *((t(1)-delay)/sd)/2);

F = [Vin 0 0 0 0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*(F' );
b=0;
for  k= 11:10:(n*10-1) %ms
    b=b+1;
    t(b) = (k-1)*1e-3;
    
    Vin=exp(-((t(b)-delay)/sd) *((t(b)-delay)/sd)/2);
    F = [Vin 0 0 0 0 0];
    if k>11
       Vt(b,:) = Hi*(F' + C/dt*Vt(b-1,:)');
    end
    out(b,:)=Vt(b,:);
    in(b) = k; %ms
   
end

Fs=1./t;

figure(12)
plot(t,out(:,6),'*')
hold on;
plot(t,out(:,1),'*')
title('Vin as a Gaussian Pulse with Increased Time Step')
xlabel('Time(ms)')
hold off
legend('Vout','Vin')

%% Question 3
% Circuit with Noise
% The differential equations for the ciruit in the time-domain are as follows:
% 1) $V1 = Vin$
%
% 2) $V2 - L(d/dt)*I3 - V3 = 0$
%
% 3) $I3 - V3/R3 - Cn*(d/dt)*V3 = In$
%
% 4) $I3 - V1/R1 + V2*(R2+R1)/(R1*R2) -C*(d/dt)*(V1-V2) = 0$
%
% 5) $V4 - a*I3 = 0$
%
% 6) $V5*(R0+R4) - R0*V4 = 0$
%
% The same equations in the frequency domain are:
%
% 1) $V1 = Vin$
%
% 2) $V2 - j \omega L*I3 - V3 = 0$
%
% 3) $I3 - V3/R3 - j \omega Cn*(d/dt)*V3 = In$
%
% 4) $I3 - V1/R1 + V2*(R2+R1)/(R1*R2) - j \omega C*(d/dt)*(V1-V2) = 0$
%
% 5) $V4 - a*I3 = 0$
%
% 6) $V5*(R0+R4) - R0*V4 = 0$


clear all

R1=1;
c=0.25;
R2=2;
L=0.2;
R3=10;
a=100;
R4=0.1;
R0=1000;

Cn=0.00001;
In=0.001*randn();

%V=[V1 V2 V3 I3 V4 V5]
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 -Cn 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];
G=[1 0 0 0 0 0; 
    0 1 -1 0 0 0; 
    0 0 -1/R3 1 0 0;
    -1/R1 ((R1+R2)/(R1*R2)) 0 1 0 0; 
    0 0 0 -a 1 0;
    0 0 0 0 -R0 (R0+R4)];


% Gaussian Pulse
n=1000;
out= zeros(n,6);
in= zeros(1,n);
noise= zeros(1,n);

delay = 0.06;%s
sd = 0.03; %s
dt = 1e-3;

Vt=zeros(n,6); %present voltage
t(1) = 0;

Vin=exp(-((t(1)-delay)/sd) *((t(1)-delay)/sd)/2);

F = [Vin 0 In 0  0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*F';

for  k= 1:1:n %ms
    t(k) = (k-1)*1e-3;
    
    In=0.001*randn();
    Vin=exp(-((t(k)-delay)/sd) *((t(k)-delay)/sd)/2);
    F = [Vin 0 In 0  0 0];
    if k>1
       Vt(k,:) = Hi*(F' + C/dt*Vt(k-1,:)');
    end
    out(k,:)=Vt(k,:);
    in(k) = k; %ms
    noise(k) = In;
   
end

Fs=1./t;

figure(13)
plot(t,out(:,6))
hold on;
plot(t,out(:,1))
title('Vout and Vin as a Gaussian Pulse - Cn = 0.00001')
xlabel('Time(s)')
hold off
legend('Vout','Vin')

Fin = fft(out(:,1));
Fout = fft(out(:,6));


figure(14)
semilogy(abs(fftshift(Fin)))
hold on
 semilogy(abs(fftshift(Fout)))
title('Frequency Response of Gaussian Pulse - Cn = 0.00001')
 xlabel('f (Hz)')
 ylabel('|P1(f)|')
legend('Vin','Vout')
hold off

% figure (15)
% hist(noise,100)
% title('Distribution of Noise')

%% Different Capacitance Values
% By increasing the value of Cn, the noise is proportionally overwhelmed
% and the signal is smoothed out. The bankdwidth of the signal does not
% appear to be affected.

Cn=0.0001;
In=0.001*randn();

%V=[V1 V2 V3 I3 V4 V5]
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 -Cn 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];


n=1000;
out= zeros(n,6);
in= zeros(1,n);
noise= zeros(1,n);

delay = 0.06;%s
sd = 0.03; %s
dt = 1e-3;

Vt=zeros(n,6); %present voltage
t(1) = 0;

Vin=exp(-((t(1)-delay)/sd) *((t(1)-delay)/sd)/2);

F = [Vin 0 In 0  0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*F';

for  k= 1:1:n %ms
    t(k) = (k-1)*1e-3;
    
    In=0.001*randn();
    Vin=exp(-((t(k)-delay)/sd) *((t(k)-delay)/sd)/2);
    F = [Vin 0 In 0  0 0];
    if k>1
       Vt(k,:) = Hi*(F' + C/dt*Vt(k-1,:)');
    end
    out(k,:)=Vt(k,:);
    in(k) = k; %ms
    noise(k) = In;
   
end

Fs=1./t;

figure(16)
plot(t,out(:,6))
hold on;
plot(t,out(:,1))
title('Vout and Vin as a Gaussian Pulse - Cn = 0.0001')
xlabel('Time(s)')
hold off
legend('Vout','Vin')

Fin = fft(out(:,1));
Fout = fft(out(:,6));


figure(17)
semilogy(abs(fftshift(Fin)))
hold on
 semilogy(abs(fftshift(Fout)))
title('Frequency Response of Gaussian Pulse - Cn = 0.0001')
 xlabel('f (Hz)')
 ylabel('|P1(f)|')
legend('Vin','Vout')
hold off
%-------

Cn=0.001;
In=0.001*randn();

%V=[V1 V2 V3 I3 V4 V5]
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 -Cn 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];

n=1000;
out= zeros(n,6);
in= zeros(1,n);
noise= zeros(1,n);

delay = 0.06;%s
sd = 0.03; %s
dt = 1e-3;

Vt=zeros(n,6); %present voltage
t(1) = 0;

Vin=exp(-((t(1)-delay)/sd) *((t(1)-delay)/sd)/2);

F = [Vin 0 In 0  0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*F';

for  k= 1:1:n %ms
    t(k) = (k-1)*1e-3;
    
    In=0.001*randn();
    Vin=exp(-((t(k)-delay)/sd) *((t(k)-delay)/sd)/2);
    F = [Vin 0 In 0  0 0];
    if k>1
       Vt(k,:) = Hi*(F' + C/dt*Vt(k-1,:)');
    end
    out(k,:)=Vt(k,:);
    in(k) = k; %ms
    noise(k) = In;
   
end

Fs=1./t;

figure(18)
% Vin and Vout for gaussian pulse
plot(t,out(:,6))
hold on;
plot(t,out(:,1))
title('Vout and Vin as a Gaussian Pulse - Cn = 0.001')
xlabel('Time(s)')
hold off
legend('Vout','Vin')

Fin = fft(out(:,1));
Fout = fft(out(:,6));


figure(19)
semilogy(abs(fftshift(Fin)))
hold on
 semilogy(abs(fftshift(Fout)))
title('Frequency Response of Gaussian Pulse - Cn =0.001')
 xlabel('f (Hz)')
 ylabel('|P1(f)|')
legend('Vin','Vout')
hold off

% Increased Time Step
%
% Using the original capacitance with an increased time step, the same
% effect occurs where the output becomes overdamped, but in this case, it
% is affected by some noise. The effect of the noise is diminished with the
% increased time step, but the relationship between Vin and Vout is also
% significantly altered.

Cn=0.00001;
In=0.001*randn();

%V=[V1 V2 V3 I3 V4 V5]
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 -Cn 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];


n=100;
out= zeros(n,6);
in= zeros(1,n);
t= zeros(1,n);

delay = 0.06;%s
sd = 0.03; %s
dt = 1e-3;

Vt=zeros(n,6); %present voltage
t(1) = 0;

Vin=exp(-((t(1)-delay)/sd) *((t(1)-delay)/sd)/2);

F = [Vin 0 0 0 0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*(F' );
b=0;
for  k= 11:10:(n*10-1) %ms
    b=b+1;
    t(b) = (k-1)*1e-3;
     
    In=0.001*randn();
    Vin=exp(-((t(b)-delay)/sd) *((t(b)-delay)/sd)/2);
    F = [Vin 0 In 0  0 0];
    if b>1
       Vt(b,:) = Hi*(F' + C/dt*Vt(b-1,:)');
    end
    out(b,:)=Vt(b,:);
    in(b) = k; %ms
    noise(b) = In;
   
end
   

Fs=1./t;

figure(20)
plot(t,out(:,6),'*')
hold on;
plot(t,out(:,1),'*')
title('Vin as a Gaussian Pulse with Increased Time Step (x10)')
xlabel('Time(s)')
hold off
legend('Vout','Vin')

Cn=0.00001;
In=0.001*randn();

%V=[V1 V2 V3 I3 V4 V5]
C=[0 0 0 0 0 0; 
    0 0 0 -L 0 0; 
    0 0 -Cn 0 0 0;
    -c c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];


n=100;
out= zeros(n,6);
in= zeros(1,n);
t= zeros(1,n);

delay = 0.06;%s
sd = 0.03; %s
dt = 1e-3;

Vt=zeros(n,6); %present voltage
t(1) = 0;

Vin=exp(-((t(1)-delay)/sd) *((t(1)-delay)/sd)/2);

F = [Vin 0 0 0 0 0];

H = (C/dt + G);
Hi = inv(H);
Vt(1,:) = Hi*(F' );
b=0;
for  k= 11:20:(n*10-1) %ms
    b=b+1;
    t(b) = (k-1)*1e-3;
     
    In=0.001*randn();
    Vin=exp(-((t(b)-delay)/sd) *((t(b)-delay)/sd)/2);
    F = [Vin 0 In 0  0 0];
    if b>1
       Vt(b,:) = Hi*(F' + C/dt*Vt(b-1,:)');
    end
    out(b,:)=Vt(b,:);
    in(b) = k; %ms
    noise(b) = In;
   
end
   

Fs=1./t;

figure(21)
plot(t,out(:,6),'*')
hold on;
plot(t,out(:,1),'*')
title('Vin as a Gaussian Pulse with Increased Time Step (x20)')
xlabel('Time(s)')
hold off
legend('Vout','Vin')

%% Question 4 - Non-Linearity
% There are likely several ways to achieve modeling a parameter  with
% polynomial-form dependence. One method would require the use of an additional, or several additional matrices with another source
% value. A matrix of this type would have many zeros in it, and strategically placed
% coefficients to create the $I^2$ and $I^3$ factors. These would be
% multiplied by the G matrix to create a new composite matrix that would
% fit in the equation $C*dV/dt +GV = F$.
% Another method would be the implementation of a B matrix $C*dV/dt +GV +
% B(V)= F$ which is typically used for non-linear situations.
% Lastly, it may be possible to do something with the wonderful polynomial
% functions that Matlab has at its disposal, along the lines of polyfit,
% but any solution of this type would almost certainly also require the use
% of a matrix to be applicable to the type of analysis in the rest of this
% simulation.
