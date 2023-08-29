% WFTHR() - Concatenate waveforms based on threshold and epoch length
%   Usage
%       [wf,idx_wf,thr] = wfthr(X,thr,L,pr)
% 
%   Inputs 
%       X = voltage recording (time samples)
%       thr = amplitude threshold [default: Quiroga et al 2004]
%                                 [default: define sign by highest number of samples]
%       L = epoch length [default: 32 samples]
%       pr = "refractory period" [default: L]
% 
%   Output 
%       wf = concatenated waveforms (obs,samples)
%       idx_wf = index of waveforms (sample)
% 
% Referências:
%   
% Autor: Danilo Benette Marques, 2022

function [wf,idx_wf,thr] = wfthr(X,thr,L,pr)

if ~isvector(X)
    error('X must be a vector (samples)')
end

X = shiftdim(X);

if nargin==1 
    thr = []; L = []; pr = [];
elseif nargin==2
    L = []; pr = [];
elseif nargin==3
    pr = [];
end

if isempty(thr)
    %get highest N default thr (+/-)
    thr = 5*median(abs(X)/0.6745);
    
    thrU = thr; NthrU = numel(find(X>thrU));
    thrL = -thr; NthrL = numel(find(X<thrL));
    
    if NthrU>NthrL
        thr = thrU;
    elseif NthrU<NthrL
        thr = thrL;
    end
    %display thr
    disp(['Waveform detection based on threshold  = ' num2str(thr)])
end
if isempty(L)
    L = 32; 
end
if isempty(pr)
    pr = L;
end

%mV threshold
if thr<0
    idx_thr = find(X<thr);
elseif thr>0
    idx_thr = find(X>thr);
elseif thr==0
    idx_thr = find(X); 
    %pensar em fazer segmentação baseada em L
    idx_wf = (1:length(X))';
    wf = X; 
    return 
end

% Threshold-based spike detection considering pr
    % turns continuous data above thr to 0/1
    idx_thr01 = zeros(size(X));
    idx_thr01(idx_thr) = 1;
    % find start of continuous data above thr
    diff01 = diff([0 ; idx_thr01]);
    idx_thr1 = find(diff01==1);
    % get starts between 'pr' interval limits
    diff1 = diff(idx_thr1);
    idx_thr1_pr = [1 ; find(diff1>pr)+1];
    
idx_wf = idx_thr1(idx_thr1_pr);

% Concatenate waveforms
preL = round(L/4);

    idx_wf(idx_wf-preL+L>length(X))=[]; %exclude if last idx_wf(iwf)-preL+L greater than end
    idx_wf(idx_wf-preL+1<preL+1)=[]; %exclude if first idx_wf(iwf)-preL+1 lesser than start

wf = zeros(numel(idx_wf),L);
for iwf = 1:numel(idx_wf)
wf(iwf,:) = X((idx_wf(iwf)-preL+1):(idx_wf(iwf)-preL+L),:)';
end

% figure,plot(wf','color','k')
%     ylabel('Amp.')
%     xlabel('Time (samples)')


end

