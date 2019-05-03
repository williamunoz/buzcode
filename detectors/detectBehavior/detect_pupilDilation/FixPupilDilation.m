%% Loading pupil dilation data
basePath = pwd;
[baseFolder,baseName] = fileparts(basePath);

figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.pupildiameter.behavior.mat']);
savfile = fullfile(basePath,[baseName,'.pupilcorrection.behavior.mat']);
load(savefile,'pupildilation');

% uncomment to re-run 
% load(savfile,'pupilcorrection');
% pupildilation.puparea_pxl = pupilcorrection.puparea_pxl_raw;

%% Interpolating data associated w/ blinking, below-size-threshold pupil
puparea_raw = pupildilation.puparea_pxl;

%Below-size-threshold pupil
smallidx = isnan(puparea_raw);
smallidx = find(smallidx');

for i = 1:length(smallidx)
    if smallidx(i)-20 > 0
        puparea_raw(smallidx(i)) = round(nanmean(puparea_raw(smallidx(i)-20:smallidx(i))));
    end
end

%Blinks
unstableframes = pupildilation.unstabledetection.unstableframes;
x = diff(unstableframes)==1;
f = find([false,x]~=[x,false]);

blinkidx = [];
for i = 1:2:length(f);
    if unstableframes(f(i))-20 > 0 & unstableframes(f(i+1))+20 < unstableframes(end)
        tempmns = [round(nanmean(puparea_raw(unstableframes(f(i))-20:unstableframes(f(i))-1))) round(nanmean(puparea_raw(unstableframes(f(i+1))+1:unstableframes(f(i+1))+20)))];
        tempvals = spline([unstableframes(f(i))-1 unstableframes(f(i+1))+1],tempmns,[unstableframes(f(i)):1:unstableframes(f(i+1))]);
        puparea_raw(unstableframes(f(i)):unstableframes(f(i+1))) = tempvals;
        blinkidx = cat(2,blinkidx,[unstableframes(f(i)):1:unstableframes(f(i+1))]);
        
        figure(1); hold on;
        plot([unstableframes(f(i))-20:1:unstableframes(f(i))-1],[puparea_raw(unstableframes(f(i))-20:unstableframes(f(i))-1)],'k');
        plot([unstableframes(f(i)):1:unstableframes(f(i+1))],tempvals,'r');
        plot([unstableframes(f(i+1))+1:1:unstableframes(f(i+1))+20],[puparea_raw(unstableframes(f(i+1))+1:unstableframes(f(i+1))+20)],'k');
    end
end

%% If needing to eliminate whisker associated noise
% inFs = 50;
% lopass = 2;
% ratio =lopass/(inFs/2) ;
% 
% puparea_raw = iosr.dsp.sincFilter(puparea_raw,ratio);

%% Interpolating data associated w/ bad detection
figure; plot(puparea_raw,'k');

dpuparea = abs(diff(puparea_raw));
jumpthresh = mean(dpuparea)+1.*std(dpuparea);
y = 1;
while ~isempty(y)
    plot(dpuparea,'k'); hold on;
    plot(get(gca,'xlim'),[1 1].*jumpthresh,'r--')
    
    title({'Click to adjust threshold',...
        'or press RETURN for current threshold'})
    [~,y] = ginput(1);
    hold off
    if isempty(y); break; end
    jumpthresh = y;
end

jumpwindow = 15;
jumpframes = [];
for i = 1:length(dpuparea)
    if dpuparea(i) > jumpthresh
        jumpframes = unique([jumpframes,[i-jumpwindow:i+jumpwindow]]);
        jumpframes(jumpframes<1 | jumpframes>length(puparea_raw))=[];
        puparea_raw(jumpframes)=nan;
    end
end

badidx = isnan(puparea_raw);
badidx = find(badidx');
for i = 1:length(badidx)
    if badidx(i)-1 > 0
        puparea_raw(badidx(i)) = round(nanmean(puparea_raw(badidx(i)-1:badidx(i))));
    end
end

% f = [5892 5931 11490 11520 20960 20990 23770 23810 24990 25060 32520 32560 41010 41220 45530 45570 48280 48390 57550 57630 71590 71670 74200 74470 74960 75000];
% badidx = [];
% for i = 1:2:length(f);
%
%     tempmns = [round(nanmean(puparea_raw(f(i)-20:f(i)-1))) round(nanmean(puparea_raw(f(i+1)+1:f(i+1)+20)))];
%     tempvals = spline([f(i)-1 f(i+1)+1],tempmns,[f(i):1:f(i+1)]);
%     puparea_raw(f(i):f(i+1)) = tempvals;
%     badidx = cat(2,badidx,[f(i):1:f(i+1)]);
%
%     figure(1); hold on;
%     plot([f(i)-100:1:f(i)-1],[puparea_raw(f(i)-100:f(i)-1)],'k');
%     plot([f(i):1:f(i+1)],tempvals,'b');
%     plot([f(i+1)+1:1:f(i+1)+100],[puparea_raw(f(i+1)+1:f(i+1)+100)],'k');
% end

figure; plot(pupildilation.puparea_pxl,'k'); hold on; plot(puparea_raw,'b');

%% Excluding data related to squinting (with loss of pupil visibility) and general detection failures.
% figure; plot(puparea_raw,'k');

trashp = [43000 length(puparea_raw)];
%trashp = [10500 11100 16000 19500 25500 length(puparea_raw)];
%trashp = []; %to be specified after evaluating pupil detection video;
trashidx = [];
for i = 1:2:length(trashp)
    puparea_raw(trashp(i):trashp(i+1)) = NaN;
    trashidx = cat(2,trashidx,[trashp(i):1:trashp(i+1)]);
end

%% Re-smoothening and re-normalizing
sf_eff = pupildilation.samplingRate;
smoothwin_s = 0.1;  %s
smoothwin_frames = round(smoothwin_s.*sf_eff); %Calculate window in frames
puparea = smooth(puparea_raw,smoothwin_frames,'moving'); %,'rloess');
puparea = puparea./nanmedian(puparea);

%Saving to behavior structure

pupildilation.data = puparea;
pupilcorrection.puparea_pxl_raw = pupildilation.puparea_pxl;
pupildilation.puparea_pxl = puparea_raw;

pupilcorrection.correctiondate = today('datetime');
pupilcorrection.correctionparams.blinkidx = blinkidx;
pupilcorrection.correctionparams.smallidx = smallidx;
pupilcorrection.correctionparams.badidx = badidx;
pupilcorrection.correctionparams.excidx = trashidx;

save(savefile,'pupildilation');
save(savfile,'pupilcorrection');

close all;

%% Pupil correction figure
figure;
subplot(2,1,2);
plot(pupildilation.timestamps,pupilcorrection.puparea_pxl_raw./nanmedian(pupilcorrection.puparea_pxl_raw),'k'); hold on;
plot(pupildilation.timestamps,puparea,'b');
hold on;
plot(pupildilation.timestamps(unstableframes),zeros(size(unstableframes)),'r.','markersize',10)
plot(pupildilation.timestamps(smallidx),zeros(size(smallidx)),'m.','markersize',10)
plot(pupildilation.timestamps(badidx),zeros(size(badidx)),'b.','markersize',10)
plot(pupildilation.timestamps(trashidx),zeros(size(trashidx)),'g.','markersize',10)
legend({'raw','corrected'},'location','northwest');
hold on;
set(gca,'xticklabel',[])
ylabel({'Pupil Diameter','(medd^-^1)'})
xlim(pupildilation.t_pulse([1 end]))

subplot(2,2,1);
[counts,centers] = hist(pupilcorrection.puparea_pxl_raw./nanmedian(pupilcorrection.puparea_pxl_raw),linspace(0,4,60));
bar(centers,counts,'FaceColor',[1 1 1],'EdgeColor','k')
hold on;
[counts,centers] = hist(puparea,linspace(0,4,60));
bar(centers,counts,'FaceColor',[1 1 1],'EdgeColor','b')
hold on;
legend({'raw','corrected'},'location','northeast');
xlim([0 5])
xlabel('Pupil Diameter')

subplot(2,2,2)
noncorr = isnan(puparea); noncorr = find(noncorr');
noncorr = length(noncorr)+1;
totalcorr = length(smallidx)+length(unstableframes)+length(badidx)+1;
if length(puparea)-totalcorr-noncorr < 0
    x =[totalcorr noncorr length(puparea)-noncorr];
else
    x =[totalcorr noncorr length(puparea)-totalcorr-noncorr];
end
labels = {'Corrected','Excluded','Good'};
pie(x,[1 0 0]);
colormap([0 0 1;      %// red
    .5 .5 .5;
    1 1 1])
legend(labels,'Location','southoutside','Orientation','horizontal');
%title(ax2,'Pupil detection outcome');

NiceSave('PupilCorrection',figfolder,baseName);

%% Concatenating & renormalizing
basePath = pwd;
bz_ConcatenateBehavior('pupildiameter',basePath);

[baseFolder,baseName] = fileparts(basePath);
savefile = fullfile(basePath,[baseName,'.pupildiameter.behavior.mat']);
load(savefile,'pupildiameter');

pupildiameter.data = pupildiameter.puparea_pxl./nanmedian(pupildiameter.puparea_pxl);

save(savefile,'pupildiameter');

%% OLD: Saving corrected eye video w/ RGB tags for blinks, non-detectable pupil, and manual bad epoch exclusion, respectively.

% pupilmaxsize = 5; %for plotting purposes
%
% figfolder = fullfile(basePath,'DetectionFigures');
% savevid = fullfile(figfolder,[baseName,'.pupilvid.corrected.avi']);
%
% vidName = fullfile(basePath,[baseName,'.avi']);
% pupilvidobj = VideoReader(vidName);
% NumberOfFrames = pupilvidobj.NumberOfFrames;
%
% SAVEVID = true;
% savevidfr = 10;
% if SAVEVID
%     if ~exist(figfolder,'dir')
%         mkdir(figfolder)
%     end
%     pupdiamVid = VideoWriter(savevid);
%     pupdiamVid.FrameRate = 1./(0.015.*savevidfr);
%     open(pupdiamVid);
% end
%
% %Loading eyeline
% savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
% load(savfile,'preprocess')
% eyeline = preprocess.eyeline;
%
% pupilbox = pupildilation.pupilbox;
% X = pupildilation.pupilxy(:,1);
% Y = pupildilation.pupilxy(:,2);
%
% pupfig = figure;
% colormap('gray');
% for ff = 1:NumberOfFrames
%
%     vidframe=read(pupilvidobj,ff);
%     if size(vidframe,3)==3
%         vidframe = rgb2gray(vidframe); %greyscale the frame
%     end
%     vidframe_orig = vidframe;
%
%     if SAVEVID && mod(ff,savevidfr)==0
%         subplot(1,2,1)
%         %yrange = [0 max(puparea)]; max normlization
%         windur = 3000; %frames
%         earlypoint = max(1,ff-windur);
%         %plot(earlypoint:ff,(puparea(earlypoint:ff)-yrange(1))./diff(yrange),'k')
%         plot(earlypoint:ff,puparea(earlypoint:ff)./nanmedian(puparea),'k')
%         hold on
%         %plot(ff,puparea(ff),'ro')
%         plot(unstableframes,zeros(size(unstableframes)),'r.','markersize',10);
%         axis square;
%         plot(smallpidx,zeros(size(smallpidx)),'g.','markersize',10);
%         plot(badidx,zeros(size(badidx)),'b.','markersize',10);
%         hold off
%         xlim([earlypoint ff])
%         ylabel({'Pup. Area', '(med.^-^1)'})
%         %if yrange(1)~=yrange(2); ylim([min(puparea) max(puparea)]); end
%         ylim([0 pupilmaxsize])
%     end
%
%      if SAVEVID && mod(ff,savevidfr)==0;
%         subplot(1,2,2)
%         imagesc(vidframe_orig); axis square;
%         hold on
%         plot(eyeline(:,1),eyeline(:,2),'b:','linewidth',0.5)
%         rectangle('Position',pupilbox(ff,:),'EdgeColor',[1 0 0],...
%             'Curvature', [1,1],'LineWidth',0.5)
%         plot(X(ff),Y(ff),'g+','markersize',3)
%         set(gca,'ytick',[]);set(gca,'xtick',[]);
%         hold off
%      end
%
%     if SAVEVID && mod(ff,savevidfr)==0
%         imgFrame = getframe(gcf);
%         writeVideo(pupdiamVid,imgFrame.cdata);
%     end
% end
%
% if SAVEVID; close(pupdiamVid); end
