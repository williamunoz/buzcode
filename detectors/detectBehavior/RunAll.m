
%% Spontaneous database and analysis

%180515_WT_EM1M3\Spont3 crashed on line 400 GetPupil... should run in debugging mode!
% spathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180515_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180515_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180515_WT_EM1M3\Spont3';...

%180602_KO_EM1M3\Spont3 failed on eventCSD analysis...
% spathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180515_WT_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180530_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180530_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180530_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180530_KO_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180601_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180601_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180601_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180601_KO_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180602_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180602_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180602_KO_EM1M3\Spont3';...

%DONE
% spathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171209_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171209_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171209_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180211_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180211_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180211_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180212_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180212_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180212_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180213_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180213_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180213_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180602_KO_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180603_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180603_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180603_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180605_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180605_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180605_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180605_WT_EM1M3\Spont4'];

%DONE!
% spathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180607_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180607_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180608_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180608_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180608_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180704_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180704_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180704_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180705_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180705_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180705_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180708_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180708_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180708_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180708_WT_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180709_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180709_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180709_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180709_KO_EM1M3\Spont4';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180710_WT_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180710_WT_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180710_WT_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180711_KO_EM1M3\Spont1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180711_KO_EM1M3\Spont2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180711_KO_EM1M3\Spont3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180711_KO_EM1M3\Spont4'];

%%
for i = 14:size(spathdatabase,1)
    
    basePath = spathdatabase(i,:);
    cd(basePath);
    
    %Basic analysis
    
%     GetPupilDilation; close all;
%     bz_LFPfromDat; close all;
%     GetWhiskFromEMG; close all;
%     DetectSlowWaves(basePath,'noSpikes',true,'sensitivity',0.4,'NREMInts',[0 Inf],'noPrompts',true); close all;
    
    %eventCSD analysis
    
    [baseFolder,baseName] = fileparts(basePath);
    
    savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
    load(savfile,'sessionInfo');
    badchannels = sessionInfo.badchannels;
    usechannels = sessionInfo.AnatGrps.Channels;
    usechannels(ismember(usechannels,badchannels))=[];
    
    lfp = bz_GetLFP('all','noPrompts',true);
    
    savfile = fullfile(basePath,[baseName,'.EMGwhisk.behavior.mat']);
    load(savfile,'EMGwhisk');
    
    savfile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);
    load(savfile,'SlowWaves');
    
    bz_eventCSD(lfp,SlowWaves.ints.UP(:,1),'channels',usechannels,'twin',[0.25 0.25],'cwin',40,'spat_sm',11,'saveName','SlowWavesCSD');
    bz_eventCSD(lfp,EMGwhisk.ints.Wh(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',11,'saveName','EMGwhiskCSD');
    close all;
    
    %LFP power analysis
%     bz_eventLFPPower(lfp,SlowWaves.ints.UP(:,1),'freqs',[100 200],'channels',usechannels,'twin',[0.25 0.25],'cwin',200,'saveName','SlowWavesLFPPower_100_200');
%     bz_eventLFPPower(lfp,SlowWaves.ints.UP(:,1),'freqs',[300 400],'channels',usechannels,'twin',[0.25 0.25],'cwin',200,'saveName','SlowWavesLFPPower_300_400');
%     
%     bz_eventLFPPower(lfp,EMGwhisk.ints.Wh(:,1),'freqs',[100 200],'channels',usechannels,'twin',[0.75 0.75],'cwin',200,'saveName','EMGwhiskLFPPower_100_200');
%     bz_eventLFPPower(lfp,EMGwhisk.ints.Wh(:,1),'freqs',[300 400],'channels',usechannels,'twin',[0.75 0.75],'cwin',200,'saveName','EMGwhiskLFPPower_300_400');
%     
%     
    %HF analysis
    
    %     [baseFolder,baseName] = fileparts(basePath);
    %     figfolder = fullfile(basePath,'DetectionFigures');
    %
    %     savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
    %     load(savfile,'sessionInfo');
    %     badchannels = sessionInfo.badchannels;
    %     usechannels = sessionInfo.AnatGrps.Channels;
    %     usechannels(ismember(usechannels,badchannels))=[];
    %
    %     lfp = bz_GetLFP('all','noPrompts',true);
    %
    %     logfilter = [100 200];
    %     higfilter = [300 400];
    %
    %     logpower = zeros(length(usechannels),1);
    %     higpower = zeros(length(usechannels),1);
    %
    %     [~,chanidx] = ismember(usechannels,lfp.channels);
    %     for ii = 1:length(usechannels)
    %         [hpData,MUAPower] = FiltNPhase(double(lfp.data(:,chanidx(ii))),logfilter,lfp.samplingRate);
    %         logpower(ii) = mean(MUAPower);
    %         [hpData,MUAPower] = FiltNPhase(double(lfp.data(:,chanidx(ii))),higfilter,lfp.samplingRate);
    %         higpower(ii) = mean(MUAPower);
    %     end
    %
    %     logpower = logpower-logpower(1);
    %     higpower = higpower-higpower(1);
    %
    %     figure;
    %     subplot(1,2,1)
    %     barh(logpower./max(logpower),'b');
    %     axis tight
    %     title('Gamma 100-200 Hz');
    %     set(gca,'YDir','reverse');
    %     ylabel('channel');xlabel('norm. power');
    %
    %     subplot(1,2,2)
    %     barh(higpower./max(higpower),'r');
    %     axis tight
    %     title('Gamma 300-400 Hz');
    %     set(gca,'YDir','reverse');
    %     ylabel('channel');xlabel('norm. power');
    %
    %     NiceSave('lohiMUADepth',figfolder,baseName)
    %
    %     close all;
end

clear;

%If re-running failed detection w/ manually specified thresholds

%basePath = pwd;

%Wh
%GetWhiskFromEMG(basePath,'Whthreshold',0.75,'NWhthreshold',0.75);

%% Touch database and analysis

%180212_KO_EM1M3\Touch1 failed on EMGwhisk... line 255
% tpathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180212_KO_EM1M3\Touch1'];

%180515_WT_EM1M3\Touch1 failed on EMGwhisk... line 255
%     tpathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180213_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180515_WT_EM1M3\Touch1';...
    
%180515_WT_EM1M3\Touch2 failed on EMGwhisk... line 255
% 'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180515_WT_EM1M3\Touch2';...
 
%DONE
% tpathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180530_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180530_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180601_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180601_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180602_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180602_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180603_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180603_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180605_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180605_WT_EM1M3\Touch2'];

%DONE
% tpathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180607_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180607_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180608_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180608_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180704_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180704_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180705_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180705_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180709_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180709_KO_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180710_WT_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180710_WT_EM1M3\Touch2';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180711_KO_EM1M3\Touch1';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180711_KO_EM1M3\Touch2'];

%%
for i = 1:size(tpathdatabase,1)
    
    basePath = tpathdatabase(i,:);
    cd(basePath);
    
%     bz_LFPfromDat; close all;
%     GetWhiskFromEMG(basePath,'PulseChannel',2,'EMGChannel',3); close all;
%     GetWhiskFromEMG(basePath,'PulseChannel',2,'EMGChannel',1,'EMGanalysis',false); close all;
%     DetectSlowWaves(basePath,'noSpikes',true,'sensitivity',0.4,'NREMInts',[0 Inf],'noPrompts',true); close all;
    
    %eventCSD analysis
            [baseFolder,baseName] = fileparts(basePath);
    
            savfile = fullfile(basePath,[baseName,'.sessionInfo.mat']);
            load(savfile,'sessionInfo');
            badchannels = sessionInfo.badchannels;
            usechannels = sessionInfo.AnatGrps.Channels;
            usechannels(ismember(usechannels,badchannels))=[];
    
            lfp = bz_GetLFP('all','noPrompts',true);
    
            savfile = fullfile(basePath,[baseName,'.EMGwhisk.behavior.mat']);
            load(savfile,'EMGwhisk');
    
            savfile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);
            load(savfile,'SlowWaves');
    
            savfile = fullfile(basePath,[baseName,'.Piezotouch.behavior.mat']);
            load(savfile,'Piezotouch');
    
            bz_eventCSD(lfp,SlowWaves.ints.UP(:,1),'channels',usechannels,'twin',[0.25 0.25],'cwin',40,'spat_sm',11,'saveName','SlowWavesCSD');
            bz_eventCSD(lfp,EMGwhisk.ints.Wh(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',11,'saveName','EMGwhiskCSD');
            bz_eventCSD(lfp,Piezotouch.ints.Touch(:,1),'channels',usechannels,'twin',[0.75 0.75],'spat_sm',11,'saveName','PiezotouchCSD');
    
            close all;
end

clear;

%If re-running failed detection w/ manually specified thresholds

%basePath = pwd;

%Wh
%GetWhiskFromEMG(basePath,'Whthreshold',1,'NWhthreshold',1,'PulseChannel',2,'EMGChannel',3);
%Touch
%GetWhiskFromEMG(basePath,'Whthreshold',0.3,'NWhthreshold',0.3,'PulseChannel',2,'EMGChannel',1,'EMGanalysis',false);

%% Main database and concatenation

%Don't run until solving issue with 180515 Spont3
% mpathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180515_WT_EM1M3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180525_KO_EM1M3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180601_KO_EM1M3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180602_KO_EM1M3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180603_WT_EM1M3';...
%     'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180605_WT_EM1M3'];

%DONE
%'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180607_WT_EM1M3';...
    
mpathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180608_KO_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180703_WT_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180704_KO_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180705_KO_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180706_WT_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180708_WT_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180709_KO_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180710_WT_EM1M3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180711_KO_EM1M3'];

%This one recorded at different depth Spont 1-3 vs Spont4-Touch2
%... so no concatenation?
%'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180530_KO_EM1M3';...
%%
for i = 1:size(mpathdatabase,1)

    basePath = mpathdatabase(i,:);
    cd(basePath);

    bz_ConcatenateDats(basePath,0,0);
    
    bz_LFPfromDat; close all;

end
