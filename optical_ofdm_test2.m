clc;
clear all;
%% Global parameters and init

M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
N_sc = 64;            % 子载波数量
N_cp = 16;            % 循环前缀
maxBitErrors = 100;    % 比特错误最大值
maxNumBits = 1e5;      % 最大对比比特数
PdBm = 0.5;
symbrate = 10;          % symbol rate [Gbaud].
tx.rolloff = 0.01;      % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'qpsk';        % modulation format
Dtot = 0.3; % residual dispersion [ps/nm] at link end
lam=1550;
normx=1.05;
% Transmission fiber
ft.length     = 100e3;       % length [m]
ft.lambda     = 1550;       % wavelength [nm] of fiber parameters
ft.alphadB    = 0.2;        % attenuation [dB/km]
ft.disp       = 17;         % dispersion [ps/nm/km] @ ft.lambda
ft.slope      = 0;          % slope [ps/nm^2/km] @ ft.lambda
ft.n2         = 0;          % nonlinear index [m^2/W]
ft.aeff       = 80;         % effective area [um^2]

% compensating fiber
fc = ft;                    % same parameters but:
fc.length     = 1e3;        % fixed length [m]
fc.alphadB    = 0.2;        % attenuation [dB/km]
fc.disp  = (Dtot*1e3 - ft.disp*ft.length)/fc.length; % [ps/nm/km] to get Dtot at end link
fc.slope      = 0;          % slope [ps/nm^2/km] @ fc.lambda

% Optical amplifier
amp.gain = (ft.length*ft.alphadB + fc.length*fc.alphadB)*1e-3; % gain [dB]

%% Rx parameters
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.oftype = 'gauss';        % optical filter type
rx.obw = Inf;               % optical filter bandwidth normalized to symbrate
rx.eftype = 'rc';       % electrical filter type
rx.ebw = 0.5;               % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;
rx.type = 'bin';            % binary pattern
%% 设置QPSK调制解调
qpskMod = comm.QPSKModulator('BitInput',true);
qpskDemod = comm.QPSKDemodulator('BitOutput',true);

%% 设置OFDM参数
ofdmMod = comm.OFDMModulator('FFTLength',N_sc,'CyclicPrefixLength',N_cp,...
    'PilotInputPort',true,...
    'PilotCarrierIndices',[12;26;40;54],...
    'InsertDCNull',true,...
    'Windowing',true,...
    'WindowLength',8);
ofdmDemod = comm.OFDMDemodulator('FFTLength',N_sc,'CyclicPrefixLength',N_cp,...
    'PilotOutputPort',true,...
    'RemoveDCCarrier',true);
ofdmDims = info(ofdmMod);
disp(ofdmMod);
showResourceMapping(ofdmMod);
numDC = ofdmDims.DataInputSize(1);
frameSize = [k*numDC 1];
SNRdb = (0:10)';

avgber = zeros(length(SNRdb),3);
berget = zeros(1,3);

Nsamp = frameSize(1);        % overall number of samples
fs = symbrate*32;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% 计算误码率
errorRate = comm.ErrorRate('ResetInputPort',true);
snrlin = SNRdb + 10*log10(k) + 10*log10(numDC/N_sc);
Plin = 10.^(PdBm/10);   % [mW]
E_laser = lasersource(Plin,lam,struct('pol','single'));
rng('default');

for m = 1:length(SNRdb)
    while berget(2) <= maxBitErrors && berget(3) <= maxNumBits
        bits_ipt = randi([0,1],frameSize);% 生成原始数据
        plt_ipt=complex(rand(ofdmDims.PilotInputSize),rand(ofdmDims.PilotInputSize));
        qpsk_data = qpskMod(bits_ipt);% 映射
        after_modu = ofdmMod(qpsk_data,plt_ipt);% OFDM
        after_modu = [after_modu;zeros(16,1)];
        
        after_iq = iqmodulator(E_laser,after_modu,struct('norm',normx));
        %通过光纤信道
        after_fiber=fiber(after_iq,ft);
        after_compensating=fiber(after_fiber,fc);
        after_amp=ampliflat(after_compensating,amp);
        myrx=after_amp;
        sigma2 = (2.5e-4)*10^(-SNRdb(m)/10)*Plin*fs/symbrate;
        noise = sqrt(sigma2/2)*(randn(Nsamp,1)+1i*randn(Nsamp,1));
        myrx.field = after_amp.field + noise; 
        %相干检测
        after_channel = rxfrontend(myrx,lam,symbrate,rx);
        after_channel = after_channel(1:80);
        [qpsk_opt,plt_opt] = ofdmDemod(after_channel);% OFDM解调
        bits_opt = qpskDemod(qpsk_opt);% 解映射
        berget = errorRate(bits_ipt,bits_opt,0);% 计算误码率
    end

    avgber(m,:) = berget;
    berget = errorRate(bits_ipt,bits_opt,1);
end

scatterplot(qpsk_data);
title('input');
scatterplot(qpsk_opt);
title('output');

bertheory = berawgn(SNRdb,'psk',M,'nondiff');
%% plot
figure;
semilogy(SNRdb,avgber(:,1),'o');
hold on;
semilogy(SNRdb,bertheory);
legend('Sim','Theory');
xlabel('snr (dB)');
ylabel('Bit Error Rate');
title('coherent optical ofdm');
grid on;
hold off;
