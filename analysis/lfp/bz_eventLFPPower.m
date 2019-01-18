function [ LFPPower ] = bz_eventLFPPower (lfp, events, varargin)

% [ LFPPower ] = bz_eventLFPpower (lfp, events, varargin)
% Calculates event-triggered (i.e. SWRs) MUA power map from a linear array of LFPs
%
% INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    events          events timestamps (in sec)
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       channels    vector with channels to inlcude. (default 1:64)
%       twin        time window around events to calculate average. Default [0.1 0.1]
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       doDetrend   Default false.
%       plotCSD     true/false. Default true.
%       plotLFP     true/false. Default true.
%       saveName    Default eventCSD.
%    =========================================================================
%
% OUTPUT:
%    eventCSD            a buzcode structure with fields eventcsd.data,
%                                                        eventcsd.timestamps
%                                                        eventcsd.samplingRate
%                                                        eventcsd.channels
%                                                        eventcsd.params
%   
%
% Munoz W - 190112
%
% TODO: Make so it can read from the binary lfp file in chunks instead of from a lfp.mat
%
%% Parse inputs

p = inputParser;
addParameter(p,'channels',[1:64],@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_sm',0,@isnumeric);
addParameter(p,'temp_sm',0,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);
addParameter(p,'cwin',[]);
addParameter(p,'lfp',[]);
addParameter(p,'freqs',[1 128],@isnumeric);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'saveName','eventCSD')

parse(p,varargin{:});
lfp = p.Results.lfp;
freqs = p.Results.freqs;
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;
cwin = p.Results.cwin;
saveMat = p.Results.saveMat;
saveName = p.Results.saveName;

%Defaults
if ~exist('basePath','var')
    basePath = pwd;
end
[baseFolder,baseName] = fileparts(basePath);

%% File Management
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');

savefile = fullfile(basePath,[baseName,'.',saveName,'.lfp.mat']);

%% Load the LFP
if isempty(lfp)
    lfp = bz_GetLFP('all','noPrompts',true);
end

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

twin = p.Results.twin*samplingRate;
events = round(events*samplingRate);

%% Filtering data

for i = 1:length(channels)
    data(:,i) = abs(hilbert(bz_Filter(data(:,i),'passband',freqs,'filter','fir1','order',4)));
%     [~,MUAPower] = FiltNPhase(double(data(:,i)),freqs,samplingRate);
%     data(:,i) = MUAPower;
end
%% Compute event-triggered MUA average

[~,chanidx] = ismember(channels,lfp.channels);

pow_temp = nan(twin(1)+twin(2)+1,length(channels),length(events));
for e = 1:length(events)
    if events(e)-twin(1) > 0 && events(e)+twin(2) < size(data,1)
        pow_temp(:,:,e) = data(events(e)-twin(1):events(e)+twin(2),chanidx);
    else
    end
end

pow_avg = nanmean(pow_temp,3)*-1;


% generate output structure
eventLFPPower.Powerdata = pow_avg;
eventLFPPower.timestamps = -twin(1):twin(2);
eventLFPPower.samplingRate = samplingRate;
eventLFPPower.channels = channels;
eventLFPPower.params.freqs = freqs;
eventLFPPower.params.spat_sm = spat_sm;
eventLFPPower.params.temp_sm = temp_sm;
eventLFPPower.params.detrend = doDetrend;

if saveMat
    save(savefile,'eventLFPPower')
end

%% Plot
taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
cmax = max(max(pow_avg));

if ~isempty(cwin)
    cmax = cwin;
end

if plotLFP
    
    figure;
    subplot(1,2,1);
    imagesc(taxis,1:size(pow_avg,2),pow_avg');hold on;
    colormap jet; caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (ms)');ylabel('channel');title('LFP Power');
    plot([0 0],[1 size(pow_avg,2)],'--k');hold on;
      
end

NiceSave(saveName,figfolder,baseName)

end
