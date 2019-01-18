%% Loading pupil dilation data

if ~exist('basePath','var')
    basePath = pwd;
end
[baseFolder,baseName] = fileparts(basePath);

savefile = fullfile(basePath,[baseName,'.pupildiameter.behavior.mat']);
load(savefile,'pupildilation');

%% Interpolating data associated w/ blinking and below-size-threshold pupil

%Blinks
unstableframes = pupildilation.unstabledetection.unstableframes;
x = diff(unstableframes)==1;
f = find([false,x]~=[x,false]);

puparea_raw = pupildilation.puparea_pxl;

for i = 1:2:length(f);
    
    tempmns = [round(nanmean(puparea_raw(unstableframes(f(i))-150:unstableframes(f(i))-1))) round(nanmean(puparea_raw(unstableframes(f(i+1))+1:unstableframes(f(i+1))+150)))];
    tempvals = spline([unstableframes(f(i))-1 unstableframes(f(i+1))+1],tempmns,[unstableframes(f(i)):1:unstableframes(f(i+1))]);
    puparea_raw(unstableframes(f(i)):unstableframes(f(i+1))) = tempvals;
    
    figure(1); hold on;
    plot([unstableframes(f(i))-100:1:unstableframes(f(i))-1],[puparea_raw(unstableframes(f(i))-100:unstableframes(f(i))-1)],'k');
    plot([unstableframes(f(i)):1:unstableframes(f(i+1))],tempvals,'r');
    plot([unstableframes(f(i+1))+1:1:unstableframes(f(i+1))+100],[puparea_raw(unstableframes(f(i+1))+1:unstableframes(f(i+1))+100)],'k');
    
end

%Below-size-threshold pupil
smallpidx = isnan(puparea_raw);
smallpidx = find(smallpidx');

for i = 1:length(smallpidx)
    
    puparea_raw(smallpidx(i)) = round(nanmean(puparea_raw(smallpidx(i)-150:smallpidx(i))));
    
end

%Re-smoothening
sf_eff = pupildilation.samplingRate;
smoothwin_s = 0.5;  %s
smoothwin_frames = round(smoothwin_s.*sf_eff); %Calculate window in frames
puparea = smooth(puparea_raw,smoothwin_frames,'moving'); %,'rloess');
puparea = puparea./nanmedian(puparea);

%% Excluding data related to squinting (with loss of pupil visibility) and general detection failures.

bad = [3500 4000 6100 6150 46000 47500 49000 52500 79000 79100 79550 79650 84000 84700]; %to be specified after evaluating pupil detection video;
badidx = [];
for i = 1:2:length(bad)
    puparea_raw(bad(i):bad(i+1)) = NaN;
    puparea(bad(i):bad(i+1)) = NaN;
    badidx = cat(2,badidx,[bad(i):1:bad(i+1)]);
end

%% Saving to behavior structure

pupildilation.data_corr = puparea;
pupildilation.puparea_pxl_corr = puparea_raw;
pupildilation.smallpidx = smallpidx;
pupildilation.badidx = badidx;

save(savefile,'pupildilation');

%% Saving corrected eye video w/ RGB tags for blinks, non-detectable pupil, and manual bad epoch exclusion, respectively.

pupilmaxsize = 5; %for plotting purposes

figfolder = fullfile(basePath,'DetectionFigures');
savevid = fullfile(figfolder,[baseName,'.pupilvid.corrected.avi']);

vidName = fullfile(basePath,[baseName,'.avi']);
pupilvidobj = VideoReader(vidName);
NumberOfFrames = pupilvidobj.NumberOfFrames;

SAVEVID = true;
savevidfr = 10;
if SAVEVID
    if ~exist(figfolder,'dir')
        mkdir(figfolder)
    end
    pupdiamVid = VideoWriter(savevid);
    pupdiamVid.FrameRate = 1./(0.015.*savevidfr);
    open(pupdiamVid);
end

%Loading eyeline
savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
load(savfile,'preprocess')
eyeline = preprocess.eyeline;

pupilbox = pupildilation.pupilbox;
X = pupildilation.pupilxy(:,1);
Y = pupildilation.pupilxy(:,2);

pupfig = figure;
colormap('gray');
for ff = 1:NumberOfFrames
    
    vidframe=read(pupilvidobj,ff);
    if size(vidframe,3)==3
        vidframe = rgb2gray(vidframe); %greyscale the frame
    end
    vidframe_orig = vidframe;
     
    if SAVEVID && mod(ff,savevidfr)==0
        subplot(1,2,1)
        %yrange = [0 max(puparea)]; max normlization
        windur = 3000; %frames
        earlypoint = max(1,ff-windur);
        %plot(earlypoint:ff,(puparea(earlypoint:ff)-yrange(1))./diff(yrange),'k')
        plot(earlypoint:ff,puparea(earlypoint:ff)./nanmedian(puparea),'k')
        hold on
        %plot(ff,puparea(ff),'ro')
        plot(unstableframes,zeros(size(unstableframes)),'r.','markersize',10);
        axis square;
        plot(smallpidx,zeros(size(smallpidx)),'g.','markersize',10);
        plot(badidx,zeros(size(badidx)),'b.','markersize',10);
        hold off
        xlim([earlypoint ff])
        ylabel({'Pup. Area', '(med.^-^1)'})
        %if yrange(1)~=yrange(2); ylim([min(puparea) max(puparea)]); end
        ylim([0 pupilmaxsize])
    end
    
     if SAVEVID && mod(ff,savevidfr)==0;
        subplot(1,2,2)
        imagesc(vidframe_orig); axis square;
        hold on
        plot(eyeline(:,1),eyeline(:,2),'b:','linewidth',0.5)
        rectangle('Position',pupilbox(ff,:),'EdgeColor',[1 0 0],...
            'Curvature', [1,1],'LineWidth',0.5)
        plot(X(ff),Y(ff),'g+','markersize',3)
        set(gca,'ytick',[]);set(gca,'xtick',[]);
        hold off
     end
    
    if SAVEVID && mod(ff,savevidfr)==0
        imgFrame = getframe(gcf);
        writeVideo(pupdiamVid,imgFrame.cdata);
    end
end

if SAVEVID; close(pupdiamVid); end
