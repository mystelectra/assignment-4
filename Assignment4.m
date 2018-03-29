%% Assignment 4 - Circuit modelling
%% Question 1
% The differential equations for the ciruit in the time-domain are as follows:
% 1) V1 = Vin
% 2) V2 + L(d/dt)*I3 - V3 = 0
% 3) V3 - I3*R3=0
% 4) I3 - V1/R1 - V2*(R2+R1)/(R1*R2) = 0
% 5) V4 - a*I3 = 0
% 6) V5*(R0+R4) -R0*V4 = 0
%
% The same equations in the frequency domain are:
% 1) V1 = Vin
% 2) V2 + L(d/dt)*I3 - V3 = 0
% 3) V3 - I3*R3=0
% 4) I3 - V1/R1 - V2*(R2+R1)/(R1*R2) -C*(d/dt)*(V1-V2)= 0
% 5) V4 - a*I3 = 0
% 6) V5*(R0+R4) -R0*V4 = 0
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
    0 0 0 L 0 0; 
    0 0 0 0 0 0;
    c -c 0 0 0 0; 
    0 0 0 0 0 0;
    0 0 0 0 0 0];
G=[1 0 0 0 0 0; 
    0 1 -1 0 0 0; 
    0 0 1/R3 -1 0 0;
    -1/R1 ((R1+R2)/(R1*R2)) 0 1 0 0; 
    0 0 0 -a 1 0;
    0 0 0 0 -R0 (R0+R4)];

% DC case
out= zeros(21,6);
in= zeros(1,21);
it=0;
for vin= -10:1:10
    it=it+1;
    F(1)=vin;
    V=G\F';
    
    out(it,:)=V;
    in(it) = vin;
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
it=0;
for w= 1:1:1000
    it=it+1;
  
    V=(G+1j*w*C)\F';
    
    out(it,:)=V;
    in(it) = w;
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
it=0;
w=pi;
 cs=c+0.05*randn(1,n); %std=0.05; c=0.25;
for a=1:n
   c=cs(a);
  C=[0 0 0 0 0 0; 
    0 0 0 L 0 0; 
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
% b) For a low pass filter, signals at frequencies below some cut-off frequency 
% are sent to the output whereas signals at frequencies above the cut-off 
% frequency are attenuated.
%
% The finite difference equation used for the numerical solution of
% transient analysis is $V_i = (F + C*V_(i-1)/dt)(C/dt +G)^(-1)$


%AC case - transient sweep  - Step
n=1000;
out= zeros(n,6);
in= zeros(1,n);

Vin=1;

F = [Vin 0 0 0 0 0];

Vt=zeros(n,6); %present voltage
Vt_1 = zeros(n+1,6); %previous voltage
Vin=0;
Vt(1,:) = inv(C +G)*F' ; %V(t=0)

for t= 2:1:n %ms
    
    if t<3
        Vin=0;
        F = [Vin 0 0 0 0 0];
        Vt(t,:) = inv(C +G)*(F' + C*(Vt(t-1,:)'));
    elseif t>=3
        Vin=1;
        F = [Vin 0 0 0 0 0];
        Vt(t,:) = inv(C +G)*(F' + C*(Vt(t-1,:)'));
    end
    out(t,:)=Vt(t,:);
    in(t) = t; %ms
   
end

figure(6)
plot(in,out(:,6))
hold on;
plot(in,out(:,1))
title('Vout')
xlabel('Time(ms)')
hold off
legend('Vout','Vin')

Fin = fft(out(:,1));
Fout = fft(out(:,6));

figure(7)
plot(real(Fin),real(Fout))
% hold on
% plot(in,real(Fout))
title('Fourier Transform - frequency response ')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
legend('Vin','Vout')


%AC case - transient sweep - Sine function
n=1000;
out= zeros(n,6);
in= zeros(1,n);

Vin=1;

F = [Vin 0 0 0 0 0];

Vt=zeros(n,6); %present voltage
%Vt_1 = zeros(n+1,6); %previous voltage
Vt(1,:) = inv(C +G)*F' ; %V(t=0)
f=1/0.03; %Hz

for t= 2:1:n %ms
    
    Vin=sin(2*pi*f*t);
    F = [Vin 0 0 0 0 0];
    
    Vt(t,:) = inv(C +G)*(F' + C*(Vt(t-1,:)'));
   
    out(t,:)=Vt(t,:);
    in(t) = t; %ms
   
end
Fs = ones(1,n)*1/n;

figure(8)
plot(in,out(:,6))
hold on;
plot(in,out(:,1))
title('Vout')
xlabel('Time(ms)')
legend('Vout','Vin')
hold off

%need to create frequency vector

Fin = fft(out(:,1));
Fout = fft(out(:,6));

figure(9)
plot(Fs,real(Fin))
hold on
plot(Fs,real(Fout))
title('Fourier Transform - frequency response ')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('Vin','Vout')



