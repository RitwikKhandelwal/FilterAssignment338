%Butterworth Analog LPF parameters
Wc = 1.08;              %cut-off frequency
N = 8;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/16) + 1i*Wc*sin(pi/2 + pi/16);
p2 = Wc*cos(pi/2 + pi/16) - 1i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + 1i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - 1i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - 1i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + 1i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - 1i*Wc*sin(pi/2 + pi/16+3*pi/8);

%Band Edge speifications
fs1 = 48.7e3;
fp1 = 44.7e3;
fp2 = 72.7e3;
fs2 = 68.7e3;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260e3;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog BPF
[Nl, Dl] = numden(analog_lpf(s));                   %numerical simplification
Nl = sym2poly(expand(Nl));                          
Dl = sym2poly(expand(Dl));                          %collect coeffs into matrix form
k = Dl(1);    
Dl = Dl/k;
Nl = Nl/k;

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;

%BANDSTOP ANALOG
[H_transfer,freq] = freqs(ns,ds,10e3);
figure,plot(freq,abs(H_transfer));grid;xlim([0 2]);
ylabel("Magnitude response");title("Bandstop Filter FT (Analog domain)");
figure,plot(freq,angle(H_transfer));grid;xlim([0 2]);
ylabel("Phase response");title("Bandstop Filter FT (Analog domain)");

%LOWPASS ANALOG
[H_transfer,freq] = freqs(Nl,Dl,10e3);
figure,plot(freq,abs(H_transfer));grid;xlim([0 4]);
ylabel("Magnitude response");title("Lowpass Filter FT (Analog domain)");
figure,plot(freq,angle(H_transfer));grid;xlim([0 4]);
ylabel("Phase response");title("Lowpass Filter FT (Analog doomain)");

%BANDSTOP DIGITAL
freq = -pi:3.1415e-4:pi;
[H_transfer] = freqz(nz,dz,freq, 2*pi);
figure,plot(freq,abs(H_transfer));grid;xlim([-pi pi]);
ylabel("Magnitude response");title("Bandstop Filter DTFT");
xticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2', '3\pi/4','\pi'});
figure,plot(freq,angle(H_transfer));grid;xlim([-pi pi]);
ylabel("Phase response");title("Bandstop Filter DTFT");
xticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2', '3\pi/4','\pi'});