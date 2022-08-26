close all;
clear all;
clc

M = 4;
bitsPerSym = log2(M);
NumofSymble = 1000;

BaudRate            = 1e6;                        
SymbolDuration      = 1/BaudRate;                   
numSamplesPerSymbol = 2000;
Fs                  = BaudRate*numSamplesPerSymbol; 
Ts                  = 1/Fs;                         

x = randi([0 M-1],NumofSymble,1);
x_bitstream = de2bi(x,4);


%%%%%%%%%% PAM Modulation %%%%%%%%%%%%%%%%%%%%
symbols = pammod(x,M,0,'bin');
symbols = symbols;
TestInput = 1:M-1;
Constellation = pammod(TestInput,M,0,'bin');
% figure();
% plot(Constellation,'.');

%%%%%%%%%% parameters for filter %%%%%%%%%%%%%%%%%%%%%%
FilterSampleFactor  = numSamplesPerSymbol;          
L                   = 6;                           
a                   = 0.5;                         
b                   = 0;                                     
t                   = -SymbolDuration*L/2:Ts:SymbolDuration*L/2;
FilterCosin         = rcosdesign(a,L,FilterSampleFactor,'sqrt');            %check 'rcosdesign' again
% impz(FilterCosin);

%up sampling
SimbolsUpSample = UpSample(symbols,numSamplesPerSymbol);        %adding (numSamplesPerSymbol) zeros

% figure();
%  plot(SimbolsUpSample,'.');
 

TxSignal = conv(SimbolsUpSample,FilterCosin);
%TxSignal_Aligned = TxSignal((length(FilterCosin)-1)/2+1:end-(length(FilterCosin)-1)/2);

TxSignal = TxSignal* 25;


%% Drive LED (Low pass filter)
FrequencyPortionRemainAfterLED = (100e6)/(Fs/2);
d  = fdesign.lowpass('N,F3dB', 50, FrequencyPortionRemainAfterLED);
flatLowpass = design(d, 'maxflat','SystemObject',true);
%fvtool(flatLowpass);
Tx_passingLED = filter(flatLowpass.Numerator,1,TxSignal);  
%Tx_passingLED = Tx;


load('IR_FL_55Point.mat');
y = conv(TxSignal,impulse_response_FL);
%y = Tx_passingLED;

figure();
for SNR = 5:5:20
rxSig(:,SNR/5) = awgn(y,SNR);
%rxSig = y;

% subplot(2,2,SNR/5);
% plot(rxSig(:,SNR/5),'.');
% title(sprintf("Constellation with noise, SNR = %d", SNR));

% %%%%%%%%%% Matched Filter %%%%%%%%%%%%%%%%%%%%
Signal_NoISI(:,SNR/5) = conv(rxSig(:,SNR/5),FilterCosin);

%Signal_NoISI_Aligned(:,SNR/5) = Signal_NoISI((length(FilterCosin)-1)/2+1:length(Signal_NoISI(:,SNR/5))-(length(FilterCosin)-1)/2);

%Signal_downsampled(:,SNR/5) = DownSample(Signal_NoISI_Aligned(:,SNR/5),numSamplesPerSymbol);
Signal_downsampled_temp(:,SNR/5) = DownSample(Signal_NoISI(:,SNR/5),numSamplesPerSymbol);
Signal_downsampled(:,SNR/5) = Signal_downsampled_temp(7:end-7,SNR/5);

temp(:,SNR/5) = Signal_downsampled(:,SNR/5)*0;

subplot(2,2,SNR/5);
plot(Signal_downsampled(:,SNR/5), temp(:,SNR/5),'.');
title(sprintf("Constellation with noise, SNR = %d", SNR));

%Signal_downsampled_Aligned = Signal_downsampled(SNR/5,(length(FilterCosin)-1)/2+1:end-(length(FilterCosin)-1)/2);

symbols_demodulated(:,SNR/5) = pamdemod(Signal_downsampled(17:end-16,SNR/5),M,0,'bin');

% subplot(2,2,SNR/5);
% plot(Signal_downsampled(17:end-16,SNR/5),0);
% title(sprintf("Constellation Diagram Before Demodulation, SNR = %d", SNR));

end
% figure();
% plot([5:5:20],BER);

% title(sprintf("BER when SNR = %d", SNR));
for i = 1:4
    bitstream_output = de2bi(symbols_demodulated(:,i),4);
    %ErrorNumber_temp = find(bitstream_output-x_bitstream);
    ErrorNumber = size(find(bitstream_output-x_bitstream),1);
    BER(i)= ErrorNumber/(NumofSymble*bitsPerSym);

%BER(SNR/5) = ErrorNumber/(NumofSymble*bitsPerSym);
end
figure();
plot([5:5:20],BER);
title('Bit Error Rate for different SNR');
xlabel('SNR');
ylabel('BER');