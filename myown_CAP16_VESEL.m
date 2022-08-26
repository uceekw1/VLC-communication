
close all;
clear all;
clc

% System Parameters
M = 16;
NumOfBitsPerSymbol = log2(M);
NumOfSymbols_InputData = 1000000;
NumofBits_InputData = NumOfSymbols_InputData*NumOfBitsPerSymbol;
SymbolRate = 2.5e9;                 % Targeted Symbol Rate
NumOfSamplesPerSymbol = 200;        % Adding Zeros Before RRC Filter (5 or 8 gives good result;28;40)

SymbolDuration = 1/SymbolRate;
SamplingFrequency = SymbolRate*NumOfSamplesPerSymbol;
SamplingPeriod = 1/SamplingFrequency;

BER = zeros(1,20);

% QAM16 Related Parameters
Row_Qam16 = [-3, -1, 1, 3];
A = repmat(Row_Qam16, 4, 1);
B = flipud(A');
Constellation_QAM_Temp = A+i*B;
Constellation_QAM = Constellation_QAM_Temp(:); 


% RRC Filter Design
FilterSampleFactor  = NumOfSamplesPerSymbol;        
L                   = 6;                           % Cutting Length of Filter
a                   = 0.5;                         % Roll Off Factor 0.95
t                   = -SymbolDuration*L/2:SamplingPeriod:SymbolDuration*L/2;
FilterCosin         = rcosdesign(a,L,FilterSampleFactor,'sqrt'); 

% CAP Filter is based on RRC filter with an extra item (coswt and sinwt)
Fc  = (1+a)/(2*SymbolDuration);
%Fc  = (1+a)/(2*SamplingPeriod);
coswt = cos(2*pi*Fc*t);
sinwt = sin(2*pi*Fc*t);
InFilter = FilterCosin.*coswt;
QnFilter = FilterCosin.*sinwt; 

% Input Data Generation
 BitStream = randi([0 1],NumofBits_InputData,1);
 %load('RandomBitStream.mat');


% QAM4 Mapping
%[I_SymbolsTx, Q_SymbolsTx] = QPSKEncoder(BitStreamOne, BitStreamTwo); 
%BitStream_Padded = [repmat('0', 1, ceil(numel(BitStream)/NumOfBitsPerSymbol)*NumOfBitsPerSymbol-numel(BitStream)) BitStream];                                                     
BitStream_Column = reshape(BitStream, NumOfBitsPerSymbol, []).';
for i = 1:size(BitStream_Column,1)
    temp = sum(BitStream_Column(i,:) .* 2.^(NumOfBitsPerSymbol-1:-1:0));
    BitStream_Column_Int(i,:) = temp;
end
BitStream_Column_Int_ForQAM = BitStream_Column_Int + 1;
Symbols = Constellation_QAM(BitStream_Column_Int_ForQAM);

% Upsampling 
SymbolsOverSampled = zeros(size(1:NumOfSymbols_InputData*NumOfSamplesPerSymbol));
SymbolsOverSampled(1:NumOfSamplesPerSymbol:NumOfSymbols_InputData*NumOfSamplesPerSymbol) = Symbols;

% I/Q separation
%[SymbolStreamOne, SymbolStreamTwo] = IQ_Separation(SymbolsOverSampled,NumOfSamplesPerSymbol);
Real_SymbolsOverSampled = real(SymbolsOverSampled);
Imag_SymbolsOverSampled = imag(SymbolsOverSampled);

% Filter I&Q Components 
% SymbolStreamOne_Filtered = conv(InFilter,SymbolStreamOne);
% SymbolStreamTwo_Filtered = conv(QnFilter,SymbolStreamTwo);
Real_SymbolsOverSampled_Filtered = conv(InFilter,Real_SymbolsOverSampled);
Imag_SymbolsOverSampled_Filtered = conv(QnFilter,Imag_SymbolsOverSampled);

% Adding I/Q Components Together. why subtract here?
% TxSignal_temp = SymbolStreamOne_Filtered-SymbolStreamTwo_Filtered;
% TxSignal_temp2 = TxSignal_temp((length(InFilter)-1)/2+1:end-(length(InFilter)-1)/2);
TxSignal_temp = Real_SymbolsOverSampled_Filtered - Imag_SymbolsOverSampled_Filtered;
TxSignal_temp2 = TxSignal_temp((length(InFilter)-1)/2+1:end-(length(InFilter)-1)/2);

Tx = TxSignal_temp2*2;     % Amplify 25 times!

% figure();
% plot(Tx);

% %% Signal Tx's Spectrum
% F_TxSignal = fft(Tx./length(Tx));
% FSingle_TxSignal = F_TxSignal(1:length(F_TxSignal)/2);
% SamplingNum = length(F_TxSignal);
% for i = 1:length(FSingle_TxSignal)
%     F_xlabel(i) = (i-1)/(SamplingNum*SamplingPeriod);
% end
% figure(1)
% plot(F_xlabel/1e6, 20*log10(abs(FSingle_TxSignal)),'r');grid on;
% % xlim([0 500]);ylim([-140 -40]);
% xlabel('Frequency(MHz)');ylabel('Power(dBm)');

%% Drive Laser (Low pass filter)
%FrequencyPortionRemainAfterLED = 0.5;   %LED cut off frequency set to 100MHz %100e6/SamplingFrequency;
FrequencyPortionRemainAfterLaser = (35e9)/(SamplingFrequency/2);
%d  = fdesign.lowpass('N,F3dB', 50, 0.5);
d  = fdesign.lowpass('N,F3dB', 50, FrequencyPortionRemainAfterLaser);
flatLowpass = design(d, 'maxflat','SystemObject',true);
%fvtool(flatLowpass);
Tx_passingLaser = filter(flatLowpass.Numerator,1,Tx);  
%Tx_passingLaser = Tx;

%% Passing Channel in room
% load('IR_FL_55Point.mat');
% Tx_passedChannel_temp = conv(impulse_response_FL,Tx_passingLaser);
% %Tx_passedChannel = Tx_passedChannel_temp((length(impulse_response_FL)-1)/2+1:end-(length(impulse_response_FL)-1)/2);
% %Tx_passedChannel = Tx_passedChannel_temp;
 Tx_passedChannel = Tx_passingLaser;

%% 经过AWGN信道
%for SNR = 1:20
SNR = 20; 
Rx = awgn(Tx_passedChannel,SNR,'measured');
%Rx = awgn(Tx_passingLED,SNR,'measured');
%Rx = awgn(Tx,SNR,'measured');
% Rx = Tx_passedChannel;

%% Receiver Filter
RxSig_I = conv(InFilter,Rx);
RxSig_Q = conv(QnFilter,Rx);

%% Addjust sampling phase shift
% RxSig_I = RxSig_I((length(InFilter)-1)/2+1:end-(length(InFilter)-1)/2);
% RxSig_Q = RxSig_Q((length(QnFilter)-1)/2+1:end-(length(QnFilter)-1)/2);
RxSig_I = RxSig_I((length(InFilter)-1)/2+28:end+27-(length(InFilter)-1)/2);   %good sampling phase shift
RxSig_Q = RxSig_Q((length(QnFilter)-1)/2+28:end+27-(length(QnFilter)-1)/2);

%Rx_Filtered = IQ_Combination(RxSig_I,RxSig_Q,NumOfSamplesPerSymbol);
Rx_Filtered = complex(RxSig_I,RxSig_Q);

%下采样
%RxTempDownSamp_temp = DownSample(RxTemp,numSamplesPerSymbol);
RxSampOut_Unaligned = Rx_Filtered(1:NumOfSamplesPerSymbol:end);
RxSampOut = RxSampOut_Unaligned(4:end-6); % Manually checked and aligned.
figure();
plot(RxSampOut,'.');
%axis([-5,5,-5,5]);
%title(sprintf("Downsampled Signal Contellation with SNR = %d",SNR));
title("Downsampled Signal Contellation with Symbol Rate 2.5 GS/s");


% %% Equalisation
% Equaliser = comm.LinearEqualizer("Algorithm","LMS",'NumTaps',3,'StepSize',0.005,'Constellation',Constellation_QAM,'ReferenceTap',1);
% %   Equaliser = comm.DecisionFeedbackEqualizer('Algorithm','LMS','NumForwardTaps',3, ...
% %        'NumFeedbackTaps',1,'ForgettingFactor',0.95,'Constellation',Constellation_QAM,'InputSamplesPerSymbol',1,'ReferenceTap',2);
% MessageReceivedLen = length(RxSampOut);
% TrainingMessageLen = floor(MessageReceivedLen/5);
% 
% TrainingMessage = Symbols(1:TrainingMessageLen);%取其中1/5训练均衡器抽头系数
% MessageReceived = RxSampOut';     
% %MessageReceived = Rx_Filtered';   
% y = Equaliser(MessageReceived, TrainingMessage);
% 
% figure();
% plot(y(TrainingMessageLen+1:TrainingMessageLen+5000),'.');
% title("Equalized Signal Contellation Alignment ");
% % % axis([-5,5,-5,5])

%Equaliser.release();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %%%%%%%%% BER without EQ %%%%%%%%%%%%%%%
% RealRx = real(RxSampOut);
% ImagRx = imag(RxSampOut);
% 
% [RealRxBits,ImagRxBits] = QPSKDecoder(RealRx,ImagRx);
% 
% RxBits = ParallelToStream(RealRxBits,ImagRxBits);
% 
% ErrorNumber = size(find(BitStream(41:length(RxBits))-RxBits(1:end-40)),2);  %%%%为什么星座图看着基本没问题了，但是ber却差这么多？？？滤波器函数的输出截取不同！ 同样的，均衡器的trainingSignal的对齐也要调整才行！！！最好把滤波器的系数拿出来手动卷积再对齐！！！
% BER = ErrorNumber/length(RxBits(1:end-40));
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     sn_block = repmat(RxSampOut,16,1);
%     msglen = length(RxSampOut);
%     konst_block = repmat(Constellation_QAM,1,msglen);
%     distance = abs(sn_block-konst_block);
%     [dmin,ind_2] = min(distance);
%     qam_det = Constellation_QAM(ind_2);
%     qamlen = length(qam_det);

%     number_of_errors = sum(qam(1:qamlen-2) ~= qam_det(1:end));  %%I mannually aligned this!!!wrong code!
%     sep_simu(i) = number_of_errors/qamlen;

Bitstream_Out_temp = qamdemod(RxSampOut,M,'bin','OutputType','bit');
Bitstream_Out = reshape(Bitstream_Out_temp,[],1);
number_of_errors = sum(BitStream(5:end-4) ~= Bitstream_Out);
%  BER(SNR) = number_of_errors/length(Bitstream_Out);
% end
% figure();
% plot(BER);
% title("Bit Error Rate with different SNR");

