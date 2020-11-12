f_samp = 260e3;

%Band Edge speifications
fs1 = 48.7e3;
fp1 = 44.7e3;
fp2 = 72.7e3;
fs2 = 68.7e3;

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

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil(1+(A-7.95) / (2.285*0.03*pi));       %empirical formula for N_min
%disp(N_min);
%Window length for Kaiser Window
n=N_min + 2;

alpha = beta/N;

%Ideal bandpass impulse response of length "n"
b = (n-1)/2;
num = [0:1:(n-1)];
m = num - b + eps;

bs_ideal = (sin(Wc1.*m)./(pi*m))+(2.*cos((pi+Wc2).*m./2).*(sin((pi-Wc2).*m./2)./(pi*m))) ;

%Kaiser Window of length "N" with shape paramter beta calculated above
kaiser_win = (kaiser(n,alpha))';

FIR_BandStop = bs_ideal.*kaiser_win;

%magnitude response
freq = 0:3.1415e-4:pi;

H_transfer = freqz(FIR_BandStop, 1,length(freq));   
figure,plot(freq,abs(H_transfer));grid;
xticks([Wp1 Ws1 Ws2 Wp2]);xticklabels({'Wp1', 'Ws1', 'Ws2', 'Wp2'});
ylabel("Magnitude Response");title("BandStop Filter DTFT");

disp(FIR_BandStop);

%Phase response
figure,plot(freq,angle(H_transfer));grid;
xticks([Wp1 Ws1 Ws2 Wp2]);xticklabels({'Wp1', 'Ws1', 'Ws2', 'Wp2'});
ylabel("Phase Response");title("BandPass Filter DTFT");