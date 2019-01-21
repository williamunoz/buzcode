%% Loading pupil dilation data

basePath = pwd;

[baseFolder,baseName] = fileparts(basePath);

figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.pupildiameter.behavior.mat']);
savfile = fullfile(basePath,[baseName,'.pupilcorrection.behavior.mat']);
load(savefile,'pupildilation');

%% Interpolating data associated w/ blinking, below-size-threshold pupil

%Blinks
unstableframes = pupildilation.unstabledetection.unstableframes;
x = diff(unstableframes)==1;
f = find([false,x]~=[x,false]);

puparea_raw = pupildilation.puparea_pxl;
blinkidx = [];
for i = 1:2:length(f);
    
    tempmns = [round(nanmean(puparea_raw(unstableframes(f(i))-20:unstableframes(f(i))-1))) round(nanmean(puparea_raw(unstableframes(f(i+1))+1:unstableframes(f(i+1))+20)))];
    tempvals = spline([unstableframes(f(i))-1 unstableframes(f(i+1))+1],tempmns,[unstableframes(f(i)):1:unstableframes(f(i+1))]);
    puparea_raw(unstableframes(f(i)):unstableframes(f(i+1))) = tempvals;
    blinkidx = cat(2,blinkidx,[unstableframes(f(i)):1:unstableframes(f(i+1))]);
    
    figure(1); hold on;
    plot([unstableframes(f(i))-100:1:unstableframes(f(i))-1],[puparea_raw(unstableframes(f(i))-100:unstableframes(f(i))-1)],'k');
    plot([unstableframes(f(i)):1:unstableframes(f(i+1))],tempvals,'r');
    plot([unstableframes(f(i+1))+1:1:unstableframes(f(i+1))+100],[puparea_raw(unstableframes(f(i+1))+1:unstableframes(f(i+1))+100)],'k');
    
end

%Below-size-threshold pupil
smallidx = isnan(puparea_raw);
smallidx = find(smallidx');

for i = 1:length(smallidx)
    
    puparea_raw(smallidx(i)) = round(nanmean(puparea_raw(smallidx(i)-20:smallidx(i))));
    
end

%% Interpolating data associated w/ bad detection

% figure; plot(puparea_raw,'k');

f = [23320 23360 56510 56730 56910 56990 57530 57640 57690 57880 69510 69600 69600 69650 102900 103100];
badidx = [];    
for i = 1:2:length(f);
    
    tempmns = [round(nanmean(puparea_raw(f(i)-20:f(i)-1))) round(nanmean(puparea_raw(f(i+1)+1:f(i+1)+20)))];
    tempvals = spline([f(i)-1 f(i+1)+1],tempmns,[f(i):1:f(i+1)]);
    puparea_raw(f(i):f(i+1)) = tempvals;
    badidx = cat(2,badidx,[f(i):1:f(i+1)]);
    
    figure(1); hold on;
    plot([f(i)-100:1:f(i)-1],[puparea_raw(f(i)-100:f(i)-1)],'k');
    plot([f(i):1:f(i+1)],tempvals,'b');
    plot([f(i+1)+1:1:f(i+1)+100],[puparea_raw(f(i+1)+1:f(i+1)+100)],'k');
end

%% Excluding data related to squinting (with loss of pupil visibility) and general detection failures.

trashp = []; %to be specified after evaluating pupil detection video;
trashidx = [];
% for i = 1:2:length(trashp)
%     puparea(trashp(i):trashp(i+1)) = NaN;
%     trashidx = cat(2,trashidx,[trashp(i):1:trashp(i+1)]);
% end

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
plot(pupildilation.timestamps,pupildilation.puparea_pxl_raw./nanmedian(pupildilation.puparea_pxl_raw),'k'); hold on; 
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
[counts,centers] = hist(pupildilation.puparea_pxl_raw./nanmedian(pupildilation.puparea_pxl_raw),linspace(0,4,60));
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
noncorr = length(noncorr)+length(trashidx)+1;
totalcorr = length(smallidx)+length(unstableframes)+length(badidx)+1;
x =[totalcorr noncorr length(puparea)-totalcorr-noncorr];
labels = {'Corrected','Excluded','Good'};
pie(x,[1 0 0]);
colormap([0 0 1;      %// red
    .5 .5 .5;
    1 1 1])
legend(labels,'Location','southoutside','Orientation','horizontal');
%title(ax2,'Pupil detection outcome');

NiceSave('PupilCorrection',figfolder,baseName);

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
