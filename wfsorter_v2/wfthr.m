% WFTHR() - Concatenate waveforms based on threshold and epoch length
%   Usage
%       [wf,idx_wf,thr] = wfthr(X,thr,L,pr,align)
% 
%   Inputs 
%       X = voltage recording (time samples)
%       thr = amplitude threshold [default: Quiroga et al 2004]
%                                 [default: define sign by highest number of samples]
%       L = epoch length [default: 32 samples]
%       pr = "refractory period" [default: L]
%       align = method to align waveforms. {'thr','mod','max','min'}
%                                 [default: thr]
% 
%   Output 
%       wf = concatenated waveforms (obs,samples)
%       idx_wf = index of waveforms (sample)
% 
% Referências:
%   
% Autor: Danilo Benette Marques, 2022

function [wf,idx_wf,thr] = wfthr(X,thr,L,pr,align)

%Error if X is not vector
if ~isvector(X)
    error('X must be a vector (samples)')
end

% Return if thr is set to 0
if nargin>=2
    if thr==0
        %pensar em fazer segmentação baseada em L
        idx_wf = (1:length(X))';
        wf = X; 
        return 
    end
end

%Define waveform alignment method
if nargin<5 | isempty(align)
    align = 'thr';
end

X = shiftdim(X);

if nargin==1 
    thr = []; L = []; pr = [];
elseif nargin==2
    L = []; pr = [];
elseif nargin==3
    pr = [];
end

%Define default parameters
if isempty(thr)
    autothr = 1;
    thr = 5*median(abs(X)/0.6745);
    disp(['Waveform detection based on threshold  = ' num2str(thr)])
else
    autothr = 0;
end
if isempty(L)
    L = 32; 
end
if isempty(pr)
    pr = L;
end

%Define samples around detection
preL = round(L/4); %before thr. to concatenate
posL = 3*round(L/4); %after thr. to detect spike maxima

%Threshold-based spike detection considering pr
    % finds continuous data above thr
    if autothr %automatic threshold
        idx_thr = find(abs(X)>thr);
    else %manual threshold
        if thr>0
            idx_thr = find(X>thr);
        elseif thr<0
            idx_thr = find(X<thr);
        end
    end
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

% Align waveforms
    idx_wf_align = zeros(numel(idx_wf),1);
    for iwf = 1:numel(idx_wf)
        if idx_wf(iwf)-preL+1>0 && idx_wf(iwf)+posL<=length(X) %if within signal length
            switch align
                case 'mod' %find maximum absolute value around thr detection
                    [~,idx_wf_align(iwf)] = max(abs(X(idx_wf(iwf)-preL+1:idx_wf(iwf)+posL))); 
                case 'max' %find maximum around thr detection
                    [~,idx_wf_align(iwf)] = max(X(idx_wf(iwf)-preL+1:idx_wf(iwf)+posL)); 
                case 'min' %find minimum around thr detection
                    [~,idx_wf_align(iwf)] = min(X(idx_wf(iwf)-preL+1:idx_wf(iwf)+posL)); 
            end
        end
    end
    
    if ~strcmp(align,'thr')
        idx_wf_align = idx_wf_align - preL+1; %align to found indices
    end
    
idx_wf = idx_wf + idx_wf_align; %adds indices

% Concatenate waveforms
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

