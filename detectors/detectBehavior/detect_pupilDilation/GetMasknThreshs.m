spathdatabase = ['R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171006_WT_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171206_WT_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171209_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171209_WT_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171209_WT_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\171211_KO_EM1M3\Spont4';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180208_KO_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180209_WT_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180210_KO_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180211_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180211_WT_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180211_WT_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180212_KO_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180212_KO_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180212_KO_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180213_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180213_WT_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180213_WT_EM1M3\Spont3';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Spont1';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Spont2';...
    'R:\rudylab\archive\William\181205_ChrmKO_Layers_Analysis\180214_WT_EM1M3\Spont3'];

%% Obtaining mean frame and minimum intensity projections for eye/pupil masks

for i = 1:size(spathdatabase,1)
    
    basePath = spathdatabase(i,:);
    cd(basePath);
    
    [baseFolder,baseName] = fileparts(basePath);
    
    vidName = fullfile(basePath,[baseName,'.avi']);
    pupilvidobj = VideoReader(vidName);
    NumberOfFrames = pupilvidobj.NumberOfFrames;
    
    %Obtaning mean frame
    tempframe1 = double(rgb2gray(read(pupilvidobj,1)));
    for ff = 2:NumberOfFrames
        tempframe1 = tempframe1 + double(rgb2gray(read(pupilvidobj,ff)));
    end
    meanframe = tempframe1./NumberOfFrames;
    meanframe = uint8(meanframe);
    %     figure(1); colormap('gray'); subplot(2,1,1); imagesc(meanframe); axis tight; title('average frame');
    preprocess.meanframe = meanframe;
    
    %Obtaning min frame
    tempframe2 = double(rgb2gray(read(pupilvidobj,1)));
    for ff = 2:NumberOfFrames
        
        tempf = double(rgb2gray(read(pupilvidobj,ff)));
        %         tempidx = find(tempf < 1);
        %         tempf(tempidx) = NaN;
        tempframe2 = cat(3,tempframe2,tempf);
        %         tempframe2 = nanmin(tempframe2,[],3);
        tempframe2 = min(tempframe2,[],3);
    end
    minframe = uint8(tempframe2);
    %     figure(1); colormap('jet'); subplot(2,1,2); imagesc(minframe); axis tight; title('minimum intensity projection');
    preprocess.minframe = minframe;
    
    savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
    save(savfile,'preprocess');
    
    close all;
end

%% Tracing eye/pupil masks

basePath = pwd;
cd(basePath);

%Loading meanframe and min int projection
[baseFolder,baseName] = fileparts(basePath);
savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
load(savfile,'preprocess')

meanframe = preprocess.meanframe;
minframe = preprocess.minframe;

%Tracing eye
figure; colormap('jet');
subplot(2,2,1);
imagesc(meanframe);
title({'Trace the Eye','Hit "ESC" for next frame'})
set(gca,'ytick',[]);set(gca,'xtick',[]);
h_eye = imfreehand;

eyeline = getPosition(h_eye);
eyemask = createMask(h_eye);
noneyeval = 255;

meanframe(~eyemask) = noneyeval;

%Tracing pupil
subplot(2,2,3);
imagesc(meanframe)
title('Trace the Pupil')
h = imfreehand;
pupilmask = createMask(h);

%Get the pixel values for pupil, eye, and not-pupil eye
pupilpixels = double(meanframe(pupilmask))./255;
irispixels  = double(meanframe(eyemask & ~pupilmask))./255;
eyepixels = double(meanframe(eyemask))./255;

intensitybins = linspace(0,1,40);
pupilhist = hist(pupilpixels,intensitybins);
irishist = hist(irispixels,intensitybins);
eyehist = hist(eyepixels,intensitybins);

%Starting intensity threshold: 2.5std above mean pupil
intensitythresh = mean(pupilpixels)+1.*std(pupilpixels);
x = 1;
while ~isempty(x)
    tentativepupil=2.*~im2bw(meanframe,intensitythresh);
    tentativenotpupil=im2bw(meanframe,intensitythresh);
    tentativemap = eyemask+tentativepupil+tentativenotpupil;
    subplot(2,2,3)
    imagesc(tentativemap)
    subplot(2,2,2)
    bar(intensitybins,eyehist)
    hold on
    plot(intensitybins,irishist,'g','linewidth',2)
    plot(intensitybins,pupilhist,'r','linewidth',2)
    plot([1 1].*intensitythresh,get(gca,'ylim'),'r--')
    xlim([0 1])
    legend('Whole Eye','Not Pupil','Pupil',...
        'location','northwest')
    
    %Show eye with over/under pixels.
    %allow user to select threshold and show new over/under
    title({'Click to adjust threshold',...
        'or press RETURN for current threshold'})
    [x,~] = ginput(1);
    hold off
    if isempty(x); break; end
    intensitythresh = x;
end

subplot(2,2,4);
imagesc(minframe); hold on;
plot(eyeline(:,1),eyeline(:,2),'r:','linewidth',2)
title({'Trace the full pupil area','Hit "ESC" for next frame'})
set(gca,'ytick',[]);set(gca,'xtick',[]);
set(gca,'YDir','reverse');
axis tight;
h = imfreehand;
pupilmask = createMask(h);

%% Saving updated preprocess.mat

preprocess.eyeline = eyeline;
preprocess.eyemask = eyemask;
preprocess.pupilmask = pupilmask;
preprocess.intensitythresh = intensitythresh;

savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
save(savfile,'preprocess');

clear; close all;

%% Obtaining mean, median and std for eye area 

for i = 1:size(spathdatabase,1)
    
    basePath = spathdatabase(i,:);
    cd(basePath);
    
    [baseFolder,baseName] = fileparts(basePath);
    
    vidName = fullfile(basePath,[baseName,'.avi']);
    pupilvidobj = VideoReader(vidName);
    NumberOfFrames = pupilvidobj.NumberOfFrames;
    
    savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
    load(savfile,'preprocess')
    
    noneyeval = 255;
    eyemask = preprocess.eyemask;
    
    meaneyepixel = zeros(1,NumberOfFrames);
    stdeyepixel = zeros(1,NumberOfFrames);
    for ff = 1:NumberOfFrames 
        vidframe=read(pupilvidobj,ff);
        vidframe(~eyemask) = noneyeval;
        
        %Detect blinking very large/small avgerage/std eye pixel intensity
        %Mean/std of the eye pixels for unstable detection
        meaneyepixel(ff) = median(double(vidframe(eyemask))./noneyeval);
        stdeyepixel(ff) = 1.4826.*mad(double(vidframe(eyemask))./noneyeval); 
    end
    
    %Saving updated preprocess.mat
    preprocess.meaneyepixel = meaneyepixel;
    preprocess.stdeyepixel = stdeyepixel;
    
    savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
    save(savfile,'preprocess');
end

%% Selecting frame instability theresholds

for i = 1:size(spathdatabase,1)
    
    basePath = spathdatabase(i,:);
    cd(basePath);
    
    [baseFolder,baseName] = fileparts(basePath);
     
    savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
    load(savfile,'preprocess');
    
    meaneyepixel = preprocess.meaneyepixel;
    stdeyepixel = preprocess.stdeyepixel;
    
    unstablemeanthresh = mean(meaneyepixel)+1.*std(meaneyepixel);
    y = 1;
    while ~isempty(y)
        plot(meaneyepixel,'k'); hold on;
        plot(get(gca,'xlim'),[1 1].*unstablemeanthresh,'r--')
        
        title({'Click to adjust threshold',...
            'or press RETURN for current threshold'})
        [~,y] = ginput(1);
        hold off
        if isempty(y); break; end
        unstablemeanthresh = y;
    end
    
    unstablestdthresh = mean(stdeyepixel)+1.*std(stdeyepixel);
    y = 1;
    while ~isempty(y)
        plot(stdeyepixel,'k'); hold on;
        plot(get(gca,'xlim'),[1 1].*unstablestdthresh,'r--')
        
        title({'Click to adjust threshold',...
            'or press RETURN for current threshold'})
        [~,y] = ginput(1);
        hold off
        if isempty(y); break; end
        unstablestdthresh = y;
    end
    
    preprocess.unstablemeanthresh = unstablemeanthresh;
    preprocess.unstablestdthresh = unstablestdthresh;
    
    savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
    save(savfile,'preprocess');
    
    close all;
end

%% MISC_OLD WAY: Loading video

% if ~exist('basePath','var')
%     basePath = pwd;
% end
% [baseFolder,baseName] = fileparts(basePath);
% 
% vidName = fullfile(basePath,[baseName,'.avi']);
% pupilvidobj = VideoReader(vidName);
% NumberOfFrames = pupilvidobj.NumberOfFrames;

%% Obtaining mean frame and minimum intensity projection for masks

% %Obtaning mean frame
% tempframe1 = double(rgb2gray(read(pupilvidobj,1)));
% for ff = 2:NumberOfFrames
%     tempframe1 = tempframe1 + double(rgb2gray(read(pupilvidobj,ff)));
% end
% meanframe = tempframe1./NumberOfFrames;
% meanframe = uint8(meanframe);
% preprocess.meanframe = meanframe;
% 
% %Obtaning min frame
% tempframe2 = double(rgb2gray(read(pupilvidobj,1)));
% for ff = 2:NumberOfFrames
%     tempframe2 = cat(3,tempframe2,double(rgb2gray(read(pupilvidobj,ff))));
%     tempframe2 = min(tempframe2,[],3);
% end
% minframe = uint8(tempframe2);
% % figure; imagesc(minframe);

%% Obtaining eye mask and tentative pupil mask for intensity threshold

% pupfig = figure;
% colormap('gray');
% subplot(2,2,1);
% imagesc(meanframe);
% % subplot(2,2,2);
% % imagesc(eyemask);
% title({'Trace the Eye','Hit "ESC" for next frame'})
% set(gca,'ytick',[]);set(gca,'xtick',[]);
% h_eye = imfreehand;
% 
% eyeline = getPosition(h_eye);
% eyemask = createMask(h_eye);
% noneyeval = 255;
% 
% meanframe(~eyemask) = noneyeval;
% 
% %Show the Masked frame
% subplot(2,2,3)
% imagesc(meanframe)
% title('Trace the Pupil')
% %Trace the pupil
% h = imfreehand;
% pupilmask = createMask(h);
% 
% %Get the pixel values for pupil, eye, and not-pupil eye
% pupilpixels = double(meanframe(pupilmask))./255;
% irispixels  = double(meanframe(eyemask & ~pupilmask))./255;
% eyepixels = double(meanframe(eyemask))./255;
% 
% intensitybins = linspace(0,1,40);
% pupilhist = hist(pupilpixels,intensitybins);
% irishist = hist(irispixels,intensitybins);
% eyehist = hist(eyepixels,intensitybins);
% 
% %Starting intensity threshold: 2.5std above mean pupil
% intensitythresh = mean(pupilpixels)+1.*std(pupilpixels);
% x = 1;
% while ~isempty(x)
%     tentativepupil=2.*~im2bw(meanframe,intensitythresh);
%     tentativenotpupil=im2bw(meanframe,intensitythresh);
%     tentativemap = eyemask+tentativepupil+tentativenotpupil;
%     subplot(2,2,3)
%     imagesc(tentativemap)
%     subplot(2,2,2)
%     bar(intensitybins,eyehist)
%     hold on
%     plot(intensitybins,irishist,'g','linewidth',2)
%     plot(intensitybins,pupilhist,'r','linewidth',2)
%     plot([1 1].*intensitythresh,get(gca,'ylim'),'r--')
%     xlim([0 1])
%     legend('Whole Eye','Not Pupil','Pupil',...
%         'location','northwest')
%     
%     %Show eye with over/under pixels.
%     %allow user to select threshold and show new over/under
%     title({'Click to adjust threshold',...
%         'or press RETURN for current threshold'})
%     [x,~] = ginput(1);
%     hold off
%     if isempty(x); break; end
%     intensitythresh = x;
% end

%% Obtaning pupil mask

% % tempidx = find(minframe < 0.0000000000001);
% % newframe = meanframe;
% % newframe(tempidx) = noneyeval;
% % figure; imagesc(newframe);
% 
% figure; imagesc(minframe);
% 
% title('Trace the Pupil')
% %Trace the pupil
% h = imfreehand;
% pupilmask = createMask(h);

%% Obtaining mean/std eye area for thresholds

% meaneyepixel = zeros(1,NumberOfFrames);
% stdeyepixel = zeros(1,NumberOfFrames);
% for ff = 1:NumberOfFrames
%     
%     vidframe=read(pupilvidobj,ff);
%     vidframe(~eyemask) = noneyeval;
%     
%     %Detect blinking very large/small avgerage/std eye pixel intensity
%     %Mean/std of the eye pixels for unstable detection
%     meaneyepixel(ff) = median(double(vidframe(eyemask))./255);
%     stdeyepixel(ff) = 1.4826.*mad(double(vidframe(eyemask))./255);
%     
% end
% 
% unstablemeanthresh = mean(meaneyepixel)+1.*std(meaneyepixel);
% y = 1;
% while ~isempty(y)
%     plot(meaneyepixel,'k'); hold on;
%     plot(get(gca,'xlim'),[1 1].*unstablemeanthresh,'r--')
%     
%     title({'Click to adjust threshold',...
%         'or press RETURN for current threshold'})
%     [~,y] = ginput(1);
%     hold off
%     if isempty(y); break; end
%     unstablemeanthresh = y;
% end
% 
% unstablestdthresh = mean(stdeyepixel)+1.*std(stdeyepixel);
% y = 1;
% while ~isempty(y)
%     plot(stdeyepixel,'k'); hold on;
%     plot(get(gca,'xlim'),[1 1].*unstablestdthresh,'r--')
%     
%     title({'Click to adjust threshold',...
%         'or press RETURN for current threshold'})
%     [~,y] = ginput(1);
%     hold off
%     if isempty(y); break; end
%     unstablestdthresh = y;
% end

%% Saving preprocessing structure

% preprocess.minframe = minframe;
% preprocess.eyeline = eyeline;
% preprocess.eyemask = eyemask;
% preprocess.pupilmask = pupilmask;
% % preprocess.shadowmask = shadowmask;
% preprocess.intensitythresh = intensitythresh;
% preprocess.meaneyepixel = meaneyepixel;
% preprocess.stdeyepixel = stdeyepixel;
% preprocess.unstablemeanthresh = unstablemeanthresh;
% preprocess.unstablestdthresh = unstablestdthresh;
% 
% savfile = fullfile(basePath,[baseName,'.preprocess.mat']);
% save(savfile,'preprocess');
% 
% clear; close all;
