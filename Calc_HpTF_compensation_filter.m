% CALCULATE HPTF COMPENSATION FILTERS
%
% This script calculates headphone transfer function (HpTF) compensation filters 
% with high-shelf regularisation in the frequency domain. Confer the following
% papers for details on this method:
% 
%  - Schärer, Z. and Lindau, A. (2009): Evaluation of Equalization Methods for
%    Binaural Synthesis. Proc. of the 126th AES Convention
%
%  - Norcross, S. G, Soulodre G. A. and Lavoie, M. C. (2004): Distortion
%    Audibility in Inverse Filtering. Proc. of the 117th AES Convention
% 
%  - Kirkeby, O., Nelson, P. A., Hamada, H. and Orduna-Bustamante, F. (1998):
%    Fast Deconvolution of Multichannel Systems using Regularization. IEEE Trans.
%    on Speech and Audio Proc. 6(2)
% 
% The HpTFs are loaded from the folder "measurements" and the filters can be
% saved in the folder "compensation_filters".

%********************************************************************************
% Copyright (c) 2016 Vera Erbes                                                 *
%                                                                               *
% Permission is hereby granted, free of charge, to any person obtaining a copy  *
% of this software and associated documentation files (the "Software"), to deal *
% in the Software without restriction, including without limitation the rights  *
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
% copies of the Software, and to permit persons to whom the Software is         *
% furnished to do so, subject to the following conditions:                      *
%                                                                               *
% The above copyright notice and this permission notice shall be included in    *
% all copies or substantial portions of the Software.                           *
%                                                                               *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     *
% THE SOFTWARE.                                                                 *
%********************************************************************************

close all; clear; clc

%% Choose parameters
clen = 2048; %target length of compensation filter
wlen = 1024; %length of window for raw IR measurements
%regularisation
fc = 6000; %upper corner frequency for regularisation high-shelf
beta = .4; %regularisation effort
%1: same filter for left/right ear averaged, 0: different filter for left/right
one_filter = 0;

%script preferences
save_filter = 0; %1: save the resulting filters, 0: do not
all_plots = 0; %1: show all plots for intermediate steps, 0: do not

%% Choose headphones
%several headphones can be simultaneously selected to be averaged
HD600           = 1;
HED415N_1       = 0;
HED415N_3       = 0;
HED415N_4       = 0;
K271            = 0;
K601_001918     = 0;
K601_001935     = 0;
K601_14539      = 0;
K612_PRO_000744 = 0;

%% Load headphone measurements
%paths to stored measurements
meas_path = 'measurements/'; %main path
hp_path = {}; %collect measurement paths for each headphone

headphones = {}; %collect headphone names
N = 10; %number of measurements per channel

if HD600
    hp_path{end+1} = [meas_path 'HD600/2014-08-12_18-31-56_pass%02d.mat'];
    %select measurements of left and right channel
    idxl = [1 2 3 4 5 6 7 8 9 10]; 
    idxr = [1 2   4 5 6 7 8 9 10];
    headphones{end+1} = 'HD600';
end

if HED415N_1
    hp_path{end+1} = [meas_path 'HED415N_1/2014-08-12_19-34-00_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1 2 3 4 5 6 7 8 9 10];
        idxr_tmp = [1 2 3 4 5 6 7 8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)]; 
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'HED415N_1';
end

if HED415N_3
    hp_path{end+1} = [meas_path 'HED415N_3/2014-08-12_19-53-45_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1   3 4 5 6   8 9 10];
        idxr_tmp = [1   3 4 5 6   8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)];
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'HED415N_3';
end

if HED415N_4
    hp_path{end+1} = [meas_path 'HED415N_4/2014-08-12_19-59-18_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1   3 4 5 6 7 8 9 10];
        idxr_tmp = [1 2 3 4 5 6 7 8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)];
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'HED415N_4';
end

if K271
    hp_path{end+1} = [meas_path 'K271/2014-08-12_18-20-09_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1 2 3 4 5 6 7   9   ];
        idxr_tmp = [1 2 3 4     7 8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)];
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'K271';
end

if K601_001918
    hp_path{end+1} = [meas_path 'K601_001918/2014-08-12_19-04-08_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1 2 3 4 5 6 7 8 9 10];
        idxr_tmp = [1 2 3 4 5 6 7 8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)];
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'K601_001918';
end

if K601_001935
    hp_path{end+1} = [meas_path 'K601_001935/2014-08-12_19-09-44_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1 2 3 4 5 6 7 8 9 10];
        idxr_tmp = [1 2 3 4 5 6 7 8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)];
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'K601_001935';
end

if K601_14539
    hp_path{end+1} = [meas_path 'K601_14539/2014-08-12_18-43-25_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1 2 3 4 5 6 7 8 9 10];
        idxr_tmp = [1 2 3 4 5 6 7 8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)];
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'K601_14539';
end

if K612_PRO_000744
    hp_path{end+1} = [meas_path 'K612_PRO_000744/2014-08-12_19-21-42_pass%02d.mat'];
    %select measurements of left and right channel
        idxl_tmp = [1 2 3 4 5 6 7 8 9 10];
        idxr_tmp = [1 2 3 4 5 6 7 8 9 10];
    if exist('idxl') %concatenate indices when several headphones are selected
        idxl = [idxl idxl_tmp+length(idxl)];
        idxr = [idxr idxr_tmp+length(idxr)];
    else
        idxl = idxl_tmp;
        idxr = idxr_tmp;
    end
    headphones{end+1} = 'K612_PRO_000744';
end

%load all
idx = 0; %init
for n = 1:length(hp_path)
    for m = 1:N
        idx = idx+1;
        load(sprintf(hp_path{n},m))
        hl_raw(:,idx) = data.ir(:,1);
        hr_raw(:,idx) = data.ir(:,2);
    end
end
N = size(hl_raw,2); %new number of measurements (for more than one headphone)
fs = data.fs; %sampling frequency

%% Window and select IRs
%normalise all IRs (for each headphone separately)
for n = 0:N/10-1
    idx = (1:10)+n*10;
    maxvalue(n+1) = max(max(abs([hl_raw(:,idx); hr_raw(:,idx)])));
    hl_raw(:,idx) = hl_raw(:,idx)./maxvalue(n+1);
    hr_raw(:,idx) = hr_raw(:,idx)./maxvalue(n+1);
end

%window IRs and truncate to target length
win = blackmanharris(wlen);
win(1:wlen/2) = 1;
win = [win; zeros(clen-wlen,1)];
hl_win = zeros(clen,N); %preallocation
hr_win = zeros(clen,N); %preallocation
for n = 1:N
    hl_win(:,n) = hl_raw(1:clen,n).*win;
    hr_win(:,n) = hr_raw(1:clen,n).*win;
end

f = (0:clen-1)/clen*fs; %freq. vector for shortened IRs

if all_plots
    %Plot all unprocessed HpIRs and window
    figure
    subplot(1,2,1)
        for n = 1:N
            plot(db(abs(hl_raw(:,n)))), hold on, grid on
        end
        plot(db(abs(win)),'k')
        xlabel('samples'), ylabel('amplitude in dB')
        axis([0 length(hl_raw) -200 0])
        title('Raw impulse responses (left) and window')
    subplot(1,2,2)
        for n = 1:N
            plot(db(abs(hr_raw(:,n)))), hold on, grid on
        end
        plot(db(abs(win)),'k')
        xlabel('samples'), ylabel('amplitude in dB')
        axis([0 length(hl_raw) -200 0])
        title('Raw impulse responses (right) and window')
  
    %Plot all unprocessed and windowed HpTFs
    fraw = (0:size(hl_raw,1)-1)/size(hl_raw,1)*fs; %freq. vector for raw IRs
    figure
    subplot(1,2,1)
        for n = 1:N
            ph1 = semilogx(fraw,db(abs(fft(hl_raw(:,n)))),'k');
            hold on, grid on
        end
        for n = 1:N
            ph2 = semilogx(f,db(abs(fft(hl_win(:,n)))),'r');
            hold on, grid on
        end
        xlabel('frequency in Hz'), ylabel('amplitude in dB')
        axis([10 fs/2 -40 20])
        legend([ph1 ph2],'unprocessed','windowed','Location','SW')
        title('Unprocessed and windowed HpTFs (left)')
    subplot(1,2,2)
        for n = 1:N
            ph1 = semilogx(fraw,db(abs(fft(hr_raw(:,n)))),'k');
            hold on, grid on
        end
        for n = 1:N
            ph2 = semilogx(f,db(abs(fft(hr_win(:,n)))),'r');
            hold on, grid on
        end
        xlabel('frequency in Hz'), ylabel('amplitude in dB')
        axis([10 fs/2 -40 20])
        legend([ph1 ph2],'unprocessed','windowed','Location','SW')
        title('Unprocessed and windowed HpTFs (right)')
end

%pick only selected windowed IRs
hl = hl_win(:,idxl);
hr = hr_win(:,idxr);

if all_plots
    %Plot all and selected HpTFs in different colours
    figure
    subplot(1,2,1)
        for n = 1:N
            ph1(n) = semilogx(f,db(abs(fft(hl_win(:,n)))));
            hold on, grid on
        end 
        for n = 1:length(idxl)
            ph2 = semilogx(f,db(abs(fft(hl(:,n)))),'k');
            hold on, grid on
        end
        xlabel('frequency in Hz'), ylabel('amplitude in dB')
        axis([10 fs/2 -40 20])
        if N == 10
            legend([ph1 ph2],'1','2','3','4','5','6','7','8','9','10',...
                'selected','Location','S')
        else
            legend(ph2,'selected','Location','S')
        end
        title('All and selected HpTFs (left)')
    subplot(1,2,2)
        for n = 1:N
            ph1(n) = semilogx(f,db(abs(fft(hr_win(:,n)))));
            hold on, grid on
        end 
        for n = 1:length(idxr)
            ph2 = semilogx(f,db(abs(fft(hr(:,n)))),'k');
            hold on, grid on
        end
        xlabel('frequency in Hz'), ylabel('amplitude in dB')
        axis([10 fs/2 -40 20])
        if N == 10
            legend([ph1 ph2],'1','2','3','4','5','6','7','8','9','10',...
                'selected','Location','S')
        else
            legend(ph2,'selected','Location','S')
        end
        title('All and selected HpTFs (right)')
end

%% Calculate complex mean of selected HpTFs
%if only one filter for both ears: concatenate IRs of left and right channel
if one_filter
    hl = [hl hr];
    hr = hl;
    
    if all_plots
        %Plot selected HpTFs of left and right channel together
        figure
            for n = 1:size(hl,2)
                semilogx(f,db(abs(fft(hl(:,n)))),'k')
                hold on, grid on
            end
            xlabel('frequency in Hz'), ylabel('amplitude in dB')
            axis([10 fs/2 -40 20])
            title('Selected HpTFs of left and right channel together')
    end
end

%normalise selected IRs
maxvalue(end+1) = max(abs([hl(:); hr(:)]));
hl = hl./maxvalue(end);
hr = hr./maxvalue(end);

%FFT of selected IRs
Hl = fft(hl,[],1);
Hr = fft(hr,[],1);

%complex mean
Hml = sum(Hl,2)/size(Hl,2);
Hmr = sum(Hr,2)/size(Hr,2);

%Plot processed and mean transfer function
figure
subplot(1,2,1)
    for n = 1:size(Hl,2)
        ph1 = semilogx(f,db(abs(Hl(:,n))),'r');
        hold on, grid on
    end
    ph2 = semilogx(f,db(abs(Hml)),'LineWidth',1,'Color','k');
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([10 fs/2 -40 20])
    legend([ph1 ph2],'selected HpTFs','mean','Location','SW')
    title('Chosen and mean HpTFs (left)')
subplot(1,2,2)
    for n = 1:size(Hr,2)
        ph1 = semilogx(f,db(abs(Hr(:,n))),'r');
        hold on, grid on
    end
    ph2 = semilogx(f,db(abs(Hmr)),'LineWidth',1,'Color','k');
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([10 fs/2 -40 20])
    legend([ph1 ph2],'selected HpTFs','mean','Location','SW')
    title('Chosen and mean HpTFs (right)')

%% Design target bandpass
fl = 20; %lower corner frequency
fh = 20000; %upper corner frequency
stopatt = 60; %stopband attenuation

%calculate beta for kaiser window to satisfy desired stopband attenuation
beta_kaiser = .1102*(stopatt-8.7);
w = kaiser(clen+1,beta_kaiser);
b = fir1(clen,[fl/(fs/2) fh/(fs/2)],w);
b = b(2:end);
D = fft(b);
D = D.';

if all_plots
    %Plot target bandpass impulse and frequency response
    figure
    subplot(1,2,1)
        plot(0:length(D)-1,ifft(D))
        grid
        xlabel('samples'), ylabel('amplitude')
        set(gca,'XLim',[0 length(D)-1])
        title('Impulse response of target bandpass')
    subplot(1,2,2)
        semilogx(f,db(abs(D)))
        grid
        xlabel('frequency in Hz'), ylabel('amplitude in dB')
        axis([10 fs/2 -40 20])
        title('Frequency response of target bandpass')
end

%% Design regularisation filter (high-shelf)
freg = [0 2000 fc fs/2]/(fs/2); %specified frequencies
G = [-20 -20 0 0]; %gain in dB at specified frequencies
g = 10.^(G/20); %linear gain at specified frequencies
b = fir2(50,freg,g);
b = [b'; zeros(clen-length(b),1)];
B = fft(b);

if all_plots
    %Plot regularisation high-shelf
    figure
        semilogx(f,db(abs(B))), hold on
        grid
        xlabel('frequency in Hz'), ylabel('amplitude in dB')
        axis([0 fs/2 -40 20])
        title('Regularisation filter')
end

%% Calculate inverse filter in frequency domain
Hcl = D.*conj(Hml)./(Hml.*conj(Hml)+beta*B.*conj(B));
Hcr = D.*conj(Hmr)./(Hmr.*conj(Hmr)+beta*B.*conj(B));

hcl = ifft(Hcl);
hcr = ifft(Hcr);

%Plot inverse filter
figure
subplot(1,2,1)
    plot(0:clen-1,hcl), hold on
    plot(0:clen-1,hcr)
    grid
    xlabel('samples'), ylabel('amplitude')
    set(gca,'XLim',[0 clen-1])
    legend('left','right','Location','SW')
    title('Impulse responses of compensation filters')
subplot(1,2,2)
    semilogx(f,db(abs(Hcl))), hold on
    semilogx(f,db(abs(Hcr)))
    grid
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([0 fs/2 -40 20])
    legend('left','right','Location','SW')
    title('Frequency responses of compensation filters')

%% Calculate compensation results
%compensation result for mean HpTF
hml_comp = conv(hcl,ifft(Hml));
hmr_comp = conv(hcr,ifft(Hmr));

%compensation result for all selected single HpTFs
hl_comp = zeros(2*clen-1,length(idxl)); %preallocation
hr_comp = zeros(2*clen-1,length(idxr)); %preallocation
for n = 1:size(hl,2)
    hl_comp(:,n) = conv(hcl,hl(:,n));
end
for n = 1:size(hr,2)
    hr_comp(:,n) = conv(hcr,hr(:,n));
end

%Plot compensation results for mean HpTF
fcomp = (0:length(hl_comp)-1)/length(hl_comp)*fs; %freq. vector for comp. results
tcomp = (0:length(hl_comp)-1)/fs*1000; %time vector for comp. results in ms
figure
subplot(1,2,1)
    semilogx(fcomp,db(abs(fft(ifft(D),length(fcomp)))),'k'), hold on
    semilogx(fcomp,db(abs(fft(hml_comp))))   
    semilogx(fcomp,db(abs(fft(hmr_comp))))
    grid
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([0 fs/2 -40 20])
    legend('target bandpass','left','right','Location','SW')
    title('Frequency responses of compensated mean HpTF')
subplot(1,2,2)
    plot(tcomp,db(hml_comp)), hold on
    plot(tcomp,db(hmr_comp))
    grid
    xlabel('time in ms'), ylabel('amplitude in dB')
    set(gca,'XLim',[0 tcomp(end)])
    legend('left','right','Location','SW')
    title('Impulse responses of compensated mean HpTF')

%Plot compensation results for all selected HpTFs
figure
subplot(1,2,1)
    semilogx(fcomp,db(abs(fft(ifft(D),length(fcomp)))),'LineWidth',1,'Color','k')
    hold on
    semilogx(fcomp,db(abs(fft(hl_comp))))
    grid
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([0 fs/2 -40 20])
    legend('target bandpass','Location','SW')
    title('Frequency responses of compensated selected HpTFs (left)')
subplot(1,2,2)
    semilogx(fcomp,db(abs(fft(ifft(D),length(fcomp)))),'LineWidth',1,'Color','k')
    hold on
    semilogx(fcomp,db(abs(fft(hr_comp))))
    grid
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([0 fs/2 -40 20])
    legend('target bandpass','Location','SW')
    title('Frequency responses of compensated selected HpTFs (right)')
    
%% Save compensation filters
if save_filter
    %create filename depending on selected headphones
    selected_hp = ''; %init
    for n = 1:length(headphones)
        selected_hp = [selected_hp '_' headphones{n}];
    end
    if one_filter
        one_filter_string = '1Filter';
    else
        one_filter_string = 'LRFilter';
    end
    
    %one array for left and right channel
    hc = [hcl hcr]; 
    
    wavwrite(hc,fs,24,['compensation_filters/hpComp' selected_hp '_' ...
        one_filter_string '.wav']);
end