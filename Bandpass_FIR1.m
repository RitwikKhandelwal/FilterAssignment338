f_samp = 330e3;

%Band Edge speifications
fs1 = 49.5e3;
fp1 =  53.5e3;
fp2 = 73.5e3;
fs2 = 77.5e3;

Ws1 = fs1*2*pi/f_samp;
Ws2  = fs2*2*pi/f_samp;
Wp1 = fp1*2*pi/f_samp;
Wp2  = fp2*2*pi/f_samp;
Wc1 = (Ws1+Wp1)/2;
Wc2 = (Ws2+Wp2)/2;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil(1+(A-7.95) / (2.285*0.025*pi));           %empirical formula for N_min

%Window length for Kaiser Window
n= N_min+4;

alpha = beta/ n;

%Ideal bandpass impulse response of length "n"
b = (n-1)/2;
num = [0:1:(n-1)];
m = num - b + eps;

bp_ideal = 2 .* cos((Wc2+Wc1).*m./2) .* sin((Wc2-Wc1).*m./2) ./ (pi*m);

%Kaiser Window of length "N" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;

%magnitude response
freq = 0:3.1415e-4:pi;

H_transfer = freqz(FIR_BandPass, 1,length(freq));   
figure,plot(freq,abs(H_transfer));grid;
xticks([Ws1 Wp1 Wp2 Ws2]);xticklabels({'Ws1', 'Wp1', 'Wp2', 'Ws2'});
ylabel("Magnitude Response");title("BandPass Filter DTFT");

disp(FIR_BandPass);

%Phase response
figure,plot(freq,angle(H_transfer));grid;
xticks([Ws1 Wp1 Wp2 Ws2]);xticklabels({'Ws1', 'Wp1', 'Wp2', 'Ws2'});
ylabel("Phase Response");title("BandPass Filter DTFT");