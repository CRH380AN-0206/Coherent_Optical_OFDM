clc;
clear all;
%% Global parameters and init

M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
N_sc = 1024;            % 子载波数量
N_cp = N_sc/4;            % 循环前缀
maxBitErrors = 2500;    % 比特错误最大值
maxNumBits = 1e5;      % 最大对比比特数
PdBm = 0.5;
symbrate = 10;          % symbol rate [Gbaud].
tx.rolloff = 0.01;      % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'qpsk';        % modulation format
Dtot = 0.3; % 残余色散 [ps/nm]连接末尾
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
rx.eftype = 'rc';           % electrical filter type
rx.ebw = 0.5;               % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;
rx.type = 'bin';            % binary pattern
%% 设置QPSK调制解调
qpskMod = comm.QPSKModulator('BitInput',true);
qpskDemod = comm.QPSKDemodulator('BitOutput',true);

%% 设置OFDM参数
guard_band=[10;10];
plt_num=32;
plt_idc=round(linspace(guard_band(1)+7,N_sc-guard_band(2)-8,plt_num))';
ofdmMod = comm.OFDMModulator('FFTLength',N_sc,'CyclicPrefixLength',N_cp,...
    'NumGuardBandCarriers',guard_band,...
    'PilotInputPort',true,...
    'PilotCarrierIndices',plt_idc,...
    'InsertDCNull',true,...
    'Windowing',true,...
    'WindowLength',8);
ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDims = info(ofdmMod);
disp(ofdmMod);
% figure(1)
% showResourceMapping(ofdmMod);
constell = comm.ConstellationDiagram('NumInputPorts',2);
numDC = ofdmDims.DataInputSize(1);
frameSize = [k*numDC 1];
SNRdb = (0:10)';
without_gb=numDC+ofdmDims.PilotInputSize(1);
val1=N_sc/without_gb;
val2=32/val1;

avgber = zeros(length(SNRdb),3);
berget = zeros(1,3);

Nsamp = ofdmDims.OutputSize(1);        % overall number of samples
fs = symbrate*32;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% 计算误码率
errorRate = comm.ErrorRate;

snrlin = SNRdb + 10*log10(k) + 10*log10(numDC/N_sc);
Plin = 10.^(PdBm/10);   % [mW]
E_laser = lasersource(Plin,lam,struct('pol','single'));
rng('default');

for m = 1:length(SNRdb)
    reset(errorRate);
    berget = zeros(1,3);
    for idx=1:1000
        bits_ipt = randi([0,1],frameSize);% 生成原始数据
        plt_ipt=complex(rand(ofdmDims.PilotInputSize),rand(ofdmDims.PilotInputSize));
        qpsk_data = qpskMod(bits_ipt);% 映射
        after_modu = ofdmMod(qpsk_data,plt_ipt);% OFDM
        tx_pw=var(after_modu);
        after_iq = iqmodulator(E_laser,after_modu,struct('norm',normx));
        %通过光纤信道
        after_fiber=fiber(after_iq,ft);
        after_compensating=fiber(after_fiber,fc);
        after_amp=ampliflat(after_compensating,amp);
        myrx=after_amp;
        trans_pw=mean(abs(after_amp.field).^2)/(val2*2);
        sigma2 = trans_pw*10^(-SNRdb(m)/10)*Plin*fs/symbrate;
        noise = sqrt(sigma2/2)*complex(randn(Nsamp,1),randn(Nsamp,1));
        myrx.field = after_amp.field + noise;
        %相干检测
        after_channel = rxfrontend(myrx,lam,symbrate,rx);
        [qpsk_opt,plt_opt] = ofdmDemod(after_channel);% OFDM解调
        chanfr_dp=plt_opt./plt_ipt;
        chanfr_int=interp1(plt_idc,chanfr_dp,guard_band(1)+1:N_sc-guard_band(2),'linear');
        chanfr_int([plt_idc;(N_sc/2+1)]-guard_band(1))=[];
        after_chanfr=qpsk_opt./chanfr_int.';
        bits_opt = qpskDemod(qpsk_opt);% 解映射
        berget = errorRate(bits_ipt,bits_opt);% 计算误码率
    end

    avgber(m,:) = berget;
%     berget = errorRate(bits_ipt,bits_opt,1);
end

bertheory = berawgn(SNRdb,'psk',M,'nondiff');
%% plot
figure(2)
pspectrum(after_modu);
figure(3)
pspectrum(after_channel);
figure(4)
semilogy(SNRdb,avgber(:,1),'o');
hold on;
semilogy(SNRdb,bertheory);
legend('Sim','Theory');
xlabel('snr (dB)');
ylabel('Bit Error Rate');
title('coherent optical ofdm');
grid on;
hold off;
constell(qpsk_data,qpsk_opt);
chanfr_int2=interp1(plt_idc,chanfr_dp,guard_band(1)+1:N_sc-guard_band(2),'pchip');
figure(5)
plot(plt_idc,chanfr_dp,'bo',...
    1+guard_band(1):length(chanfr_int2)+guard_band(2),chanfr_int2,'b-',...
    [guard_band(1) guard_band(1)],[0 1],'r-',...
    [N_sc-guard_band(2) N_sc-guard_band(2)],[0 1],'r-');
xlim([1 N_sc]);
grid on;
ylabel('Channel Frequency Response'); xlabel('Subcarrier');