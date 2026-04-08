function [h,meanBead] = psfeval (filename, showPlot, zoomIn, redThreshold, showMicroscopeEvalPlot, showZProjection, beadSize, extractMetaDataFromPath)
%%psfeval2('test.czi','on')
%%PSF Eval determines a beads size in a multibeads measurement. It detects
%%automaticly all beads, produces a mean bead and performs a pca on it for
%%getting the main axes if it is tilted in any direction. The output h is a
%%list with multiple axis handles for information about it. showPlot should
%%be 1 for'on' and 0 for 'off'
%%same goes for zoomIn, which simply zooms in the contourplots for a 3
%%sigma boundary
%%defaults: showPlot = 1, zoomIn = 0
%%h: handle to figure collection
%%1-6: projections and fits
%%7: compact version for microscope evaluation
%%8: z-projection with used beads and filtered beads
%%mean Bead: mean bead
%%Comments: -Using CZI-Files is recommended but not necessary
%% psfeval('200917.czi',1,1, 1.7, 1, 1, 0, 1);
%%by Florian Vollrath 2015

%% check for input variables
disp('Start PSF Evalution For: ');
disp(filename);
if nargin < 8
    extractMetaDataFromPath = 0;
end

if nargin < 7
    beadSize = 0;
end

if nargin < 6
    showZProjection = 0;
end

if nargin < 5
    showMicroscopeEvalPlot = 0;
end

if nargin < 4
    redThreshold = 1.5;
end

if nargin < 3
    zoomIn = 1;
end

if nargin < 2
    showPlot = 1;
end

if nargin == 0
    error('You have to specify at least a filename');
end

%initialize waitbar for indicating the status of evaluation
wbar = waitbar(0,'please wait...');
%input variables
binningfactor=1;
brightnessTolerance=0.25;
%Read Image
waitbar (0,wbar,'Loading Image file into memory...');
[data, metadata]=imageIORead(filename);
metadata.scaleSize = metadata.scaleSize*10^6;
%get datatype of Image
data_info=whos('data');
datatype=data_info.class;

%datatype='uint32';

%when necessary information is missing in metadata ask the user for it
if isnan(metadata.wavelengthExc)
    if extractMetaDataFromPath
        if contains(filename, '488')
            metadata.wavelengthExc = 488;
        elseif contains(filename, '561')
            metadata.wavelengthExc = 561;
        elseif contains(filename, '405')
            metadata.wavelengthExc = 405;
        elseif contains(filename, '633')
            metadata.wavelengthExc = 633;
        else
            metadata.wavelengthExc = 488;
        end
        disp(['Detected ' num2str(metadata.wavelengthExc) ' as exc wavelength from path!']);
    else
        answer = inputdlg('Excitation Wavelength','Missing Variable');
        metadata.wavelengthExc = str2double(answer{1});
        if isnan(metadata.wavelengthExc)
            %set some standard value
            metadata.wavelengthExc = 488;
        end
    end
end

if isnan(metadata.refractiveIndex)
    if extractMetaDataFromPath
        if contains(filename, 'Air')
            metadata.refractiveIndex = 1;
        elseif contains(filename, 'Oil')
            metadata.refractiveIndex = 1.518;
        elseif contains(filename, 'Water')
            metadata.refractiveIndex = 1.33;
        elseif contains(filename, 'Multi')
            metadata.refractiveIndex = 1.518; %assume oil..
        else
            metadata.refractiveIndex = 1;
        end
    else
        disp(['Detected refractive index: ' num2str(metadata.refractiveIndex)]);
        answer = inputdlg('Refractive Index','Missing Variable');
        metadata.refractiveIndex = str2double(answer{1});
        if isnan(metadata.refractiveIndex)
            metadata.refractiveIndex = 1;
        end
    end
end


if isnan(metadata.NA)
    if extractMetaDataFromPath
        objectiveString = strsplit(filename, filesep);
        objectiveString = objectiveString{end-3};
        NASplit = strsplit(objectiveString, '_');
        metadata.NA = str2double(NASplit{2});
    else
        answer = inputdlg('Numerical Aperture','Missing Variable');
        metadata.NA = str2double(answer{1});
        if isnan(metadata.NA)
            metadata.NA = 1;
        end
    end
end
if metadata.refractiveIndex.^2 < metadata.NA
    warning('RI^2 is smaller than NA in meta data. This is not possible. Setting RI to 1.518 (Immersion Oil)');
    metadata.refractiveIndex = 1.518;
end

if any([isnan(metadata.scaleSize(1)), isnan(metadata.scaleSize(2)), isnan(metadata.scaleSize(3))])
    answer = inputdlg({'X Scale','Y Scale', 'Z Scale'},'Missing Variable');
    metadata.scaleSize(1) = str2double(answer{1});
    metadata.scaleSize(2) = str2double(answer{2});
    metadata.scaleSize(3) = str2double(answer{3});
end

%Leica files are given in µm, adjust for this case
if any(metadata.scaleSize > 10000)
    metadata.scaleSize = metadata.scaleSize/1000000;
end

%determine optimal beadsize:
%resolution in lateral and axial dimension
res_lat=0.61*metadata.wavelengthExc/metadata.NA;
res_ax=0.88*metadata.wavelengthExc/(metadata.refractiveIndex-sqrt(metadata.refractiveIndex.^2-metadata.NA.^2));
%sz in pixels
szFac = 3;
sz=round([res_lat/metadata.scaleSize(1)/1000, res_lat/metadata.scaleSize(2)/1000, res_ax/metadata.scaleSize(3)/1000]*szFac);
sz(3) = sz(3) / szFac * 2;

%if z-stack is not broad enough we have to truncate the z-size
if ((sz(3)*2)+1)>=size(data,3)
    sz(3)=round((size(data,3)-1)/2*0.66);
end
%check for airy processing/deconvolution
[~,fn,~] = fileparts(filename);
if any(strfind(fn, 'Airy')) || any(strfind(fn, 'Airy'))
    %then make sz (size of cut out bead) smaller
    sz = [sz(1), sz(2), round(sz(3)/8*5)];
    disp('Detected Airy Processed file, make assumed beads smaller');
end

%for performance issues we are using the same datatype as the data is in
%and not double, so we need to cast sz to this datatype
sz=cast(sz,datatype);
%%
%get mean pixel intensity as offset value and set it constant in fit..
offsetValue = mean(data(:));
data = data-offsetValue;
%% pre processing of data
%binning
if binningfactor~=1
    data_b=double(binning(data,binningfactor));
else
    data_b=double(data);
end
%gaussian filter
data_b = imgaussfilt(data_b, 2);

%threshold imagestack
%[~, thr] = threshold(data_b,'triangle');
%make thresholded binary data representation
%data_b_thr = data_b>thr;

%expects values between 0 and 1 -> normalize
minInt = min(data_b(:));
maxInt = max(data_b(:));
data_b = (data_b-minInt)/(maxInt-minInt);
%subtract mean intensity and set everything below 0 to 0
data_b = data_b - mean(data_b(:));
data_b(data_b<0)=0;
data_b = max(data_b, [], 3);

data_b_thr = imbinarize(data_b, 'global');


%close data again to avoid holes in beads
data_b_thr=imclose(data_b_thr,strel('disk',cast(sz(1)/20,'double')));
%perform opening for deleting open pixels
data_b_thr=imopen(data_b_thr,strel('disk',cast(sz(1)/20,'double')));



%label the binary array
%data_l = double(label(data_b_thr,3,minSize,1000000));

CC = bwconncomp(data_b_thr);

%delete elements smaller than minSize^3
%minSize for labeling (experimental value)
% minSize=sz(3);
% ct = 1;
% for i=1:CC.NumObjects
%     if numel(CC.PixelIdxList{i}) < minSize^3
%         deleteList(ct) = i;
%         ct = ct + 1;
%     end
% end
% CC.PixelIdxList(deleteList) = [];
% CC.NumObjects = numel(CC.PixelIdxList);

%produce labeled data set
data_l = zeros(size(data,1),size(data,2),size(data,3), 'uint16');
for i=1:CC.NumObjects
    data_l(CC.PixelIdxList{i}) = i;
end
for i = 1:size(data,3)
    data_l(:,:,i) = data_l(:,:,1);
end
%% Extract beads and produce a mean bead

%iterate through all labelnumbers
%Preallocate some arrays
numBeads = max(data_l(:));
beadCenter = zeros([3 numBeads],'uint16');
sz = cast(sz, 'uint16');
%get beadCenter
for i=1:numBeads
    waitbar (i/numBeads,wbar,['Getting Bead Centers... ',num2str(i), ' of ',...
        num2str(numBeads)]);
    % Find center as point of max intensity (in non binned data)
    [x,y,z]=ind2sub(size(data),find(data_l==i));
    v=zeros(1,length(x));
    for p=1:length(x)
        v(p)=data(x(p),y(p),z(p));
    end
    [~,maxpixel]=nanmax(v);
    maxpixel=[x(maxpixel) y(maxpixel) z(maxpixel)];
    
    binningfactorTemp = binningfactor;
    binningfactor = 1;
    beadCenter(:,i)=binningfactor*maxpixel-round(binningfactor/2)+1;
    binningfactor = binningfactorTemp;
    
end

%sometimes sz(3)is way too big, because the center of the beads isnt even
%near the center of the measurementfile, in this case we should truncate
%sz(3) even more
%if the mean of all bead centers is too near to the borders sz(3) has to be
%recalculated
if mean(beadCenter(3,:))+sz(3)>size(data,3) || mean(beadCenter(3,:))-sz(3)<0
    %check where mean of beadCenters is in the 3rd dimension
    sz(3)=cast(nanmin(abs(mean(beadCenter(3,:))+[-1 -size(data,3)])),datatype);
end

%preallocate
beads = zeros([(2*sz+1) max(data_l(:))]);
beads_l=zeros([(2*sz+1) max(data_l(:))]);
beadCenterCt = 1;
k=1;

thrownOutCounter = zeros(1,4); %variable to count how many detected beads got thrown out and why
%(beads too near to border, beads too near to each other, saturated intensity, uncommon intensity)
for i=1:max(data_l(:))
    waitbar (i/max(data_l(:)),wbar,['Extracting PSFs... ',num2str(i), ' of ',...
        num2str(max(data_l(:)))]);
    %extract bead from data into beads, but only if it is not touching or
    %overlapping any border, and there is no other labled bead in this
    %array
    %take maximum in labled data and not cog for huge performance
    %optimization
    
    %if bead is not too near to a border
    if all(beadCenter(:,i)'>sz) && all(beadCenter(:,i)'+sz<=size(data))
        %extract labeled beads in extraction volume
        beads_l(:,:,:,i)=data_l(beadCenter(1,i)-sz(1):beadCenter(1,i)+sz(1),...
            beadCenter(2,i)-sz(2):beadCenter(2,i)+sz(2),...
            beadCenter(3,i)-sz(3):beadCenter(3,i)+sz(3));
        %if the wanted labelnumber is the only one in this volume copy the
        %real data in this volume
        if all(beads_l(:,:,:,i)==i | beads_l(:,:,:,i)==0)
            beads(:,:,:,k) = data(beadCenter(1,i)-sz(1):beadCenter(1,i)+sz(1),...
                beadCenter(2,i)-sz(2):beadCenter(2,i)+sz(2),...
                beadCenter(3,i)-sz(3):beadCenter(3,i)+sz(3));
            k=k+1;
            beadCentersUsedBeads{beadCenterCt} = beadCenter(:,i);
            beadCenterCt = beadCenterCt + 1;
        else
            %increase counter for beads too near to each other
            thrownOutCounter(2) = thrownOutCounter(2) + 1;
        end
    else
        %increase counter for beads too near to edges of data
        thrownOutCounter(1) = thrownOutCounter(1) + 1;
    end
end

%delete empty data
beads=beads(:,:,:,1:k-1);

%check if there are beads left after filtering
if isempty(beads)
    error(['There are no beads left after filtering. Too close to border: ' num2str(thrownOutCounter(1)) ', too close to each other ' num2str(thrownOutCounter(2)) ', saturated: ' num2str(thrownOutCounter(3)) ', unusual intensities: ' num2str(thrownOutCounter(4))]);
end

%delete beads with saturated pixel
saturationValue = intmax(data_info.class);
for i=1:size(beads,4)
    if any(any(any(beads(:,:,:,i)==saturationValue)))
        saturationList(i) = 1;
        %increase counter for saturated intensity
        thrownOutCounter(3) = thrownOutCounter(3) + 1;
    else
        saturationList(i) = 0;
    end
end
beads(:,:,:,logical(saturationList))=[];
%filter positions too
for i = 1:size(saturationList)
    if saturationList(i)
        beadCentersUsedBeads{i} = [];
    end
end

%check for brightness and throw out the ones with much higher or lower
waitbar (0,wbar,'Delete bad PSFs');
sums=squeeze(sum(sum(sum(beads))));
referenceInt = median(sums);
sumsLogic = (sums<(1+brightnessTolerance)*referenceInt & sums>(1-brightnessTolerance)*referenceInt);
beads(:,:,:,~sumsLogic)=[];
%filter out positions that are thrown out due intensity limits
for i = 1:size(sumsLogic)
    if ~sumsLogic(i)
        thrownOutCounter(4) = thrownOutCounter(4) + 1;
        beadCentersUsedBeads{i} = [];
    end
end

%tell user how many beads got filtered out
disp(['Averaged beads: ' num2str(size(beads,4))  '; filtered beads: Too close to border: ' num2str(thrownOutCounter(1)) ', too close to each other ' num2str(thrownOutCounter(2)) ', saturated: ' num2str(thrownOutCounter(3)) ', unusual intensities: ' num2str(thrownOutCounter(4))]);

%align beads with subpixel precission
alignedBeads=align_psf(beads,[res_lat/metadata.scaleSize(1)/1000, res_lat/metadata.scaleSize(2)/1000, res_ax/metadata.scaleSize(3)/1000]);
alignedBeads(alignedBeads(:,:,:,:)<0)=0;

%produce mean bead
waitbar (0,wbar,'Produce mean PSF');
meanBead=nanmean(alignedBeads,4);

%meanBead probably contains NaN planes, we remove them:

%xy planes
delI=zeros(1,1);
k=1;
for i=1:size(meanBead,3)
    if all(isnan(meanBead(:,:,i)))
        delI(k)=i;
        k=k+1;
    end
end
if delI~=0
    meanBead(:,:,delI)=[];
end
%yz-planes
delI=zeros(1,1);
k=1;
for i=1:size(meanBead,1)
    if all(isnan(squeeze(meanBead(i,:,:))))
        delI(k)=i;
        k=k+1;
    end
end
if delI~=0
    meanBead(delI,:,:)=[];
end
%xz-planes
delI=zeros(1,1);
k=1;
for i=1:size(meanBead,2)
    if all(isnan(squeeze(meanBead(:,i,:))))
        delI(k)=i;
        k=k+1;
    end
end
if delI~=0
    meanBead(:,i,:)=[];
end
%% perform pca on meanBead to get the rotation
%get threshold of meanBead
%[~, meanThr]=threshold(meanBead,'triangle');
meanThr = 2*mean(meanBead(:));
%center of gravity/point of highest intensity
%cog will be the point of rotation
[~, cog]=nanmax(meanBead(:));
[cog(1), cog(2), cog(3)]=ind2sub(size(meanBead),cog);

%preallocate memory for PCAData
PCAData=zeros(size(meanBead,1)*size(meanBead,2)*size(meanBead,3),4);

[xind, yind, zind]= ind2sub(size(meanBead),1:numel(meanBead));
PCAData(:,1:3) = cat(1,xind,yind,zind)';
PCAData(:,4) = meanBead(:);

%Delete data pointss with NaN in it
PCAData(isnan(PCAData(:,4)), :) = [];

%Threshold PCAData for less noise influence
PCAData=PCAData(PCAData(:,4)>1*meanThr,:);

%weigth the xyz values by intensity and scaling
PCADataWeights=PCAData(:,4);

%save the scales in scales for easier use
scales=[metadata.scaleSize(1) metadata.scaleSize(2) metadata.scaleSize(3)];

%calculate the PCA data points
PCAData=[PCAData(:,1)*scales(1), PCAData(:,2)*scales(2), PCAData(:,3)*scales(3)];

%substract cog of coordinates
for i=1:3
    PCAData(:,i)=PCAData(:,i)-cog(i)*scales(i);
end

%calculate covariance matrix
[U,~,PCAVar]=pca(PCAData,'Weights',PCADataWeights);

%plot pca axes
fig=figure('Visible','off');
ax(1)=subplot(2,3,4);
hold all
isosurface(((1:size(meanBead,2))-cog(2))*scales(2),((1:size(meanBead,1))-cog(1))*scales(1),((1:size(meanBead,3))-cog(3))*scales(3),meanBead,1*meanThr); %1*meanThr
lineLength=3*([res_ax, res_lat, res_lat])/1000;

lLfac=1;
for i=1:3
    line([-lineLength(i)*U(2,i)*lLfac lineLength(i)*lLfac*U(2,i)],[-lineLength(i)*lLfac*U(1,i) lineLength(i)*lLfac*U(1,i)],[-lineLength(i)*lLfac*U(3,i) lineLength(i)*lLfac*U(3,i)],'Color',[i==1 i==2 i==3]);
end
axis equal
hold off
%% linetrace through the pca lines / right scaling
%number of interpolated points in the trace (whole image)
steps=max(size(data));
intVal=zeros(steps,6);
for j=1:3
    %strongest component in U(:,j)
    [~,maxIdx]=max(abs(U(:,j)));
    tmin=-(cog(maxIdx)*scales(maxIdx)/U(maxIdx,j));
    tmax=size(meanBead,maxIdx)*scales(maxIdx)/U(maxIdx,j)+tmin;
    %interpolate data points on line
    intVal(:,2*j)=interp3(((1:size(meanBead,2))-cog(2))*scales(2),...
        ((1:size(meanBead,1))-cog(1))*scales(1),...
        ((1:size(meanBead,3))-cog(3))*scales(3),meanBead,...
        linspace(tmin*U(2,j),tmax*U(2,j),steps),...
        linspace(tmin*U(1,j),tmax*U(1,j),steps),...
        linspace(tmin*U(3,j),tmax*U(3,j),steps),...
        'cubic');
    %add scaling to intVal
    intVal(:,1+2*(j-1))=linspace(tmin,tmax,steps);
end
%% fit 1D Gaussian PCA/ICA
sigmasTH=[res_ax,res_lat,res_lat]*1.2/1000;
sigmas=zeros(1,3);
ax(2)=subplot(2,3,5);
hold all
axisscaling=zeros(3,4);

for i=1:3
    % Set up fittype and options.
    [xData, yData] = prepareCurveData( intVal(:,1+2*(i-1)), intVal(:,2*i));
    %ft = fittype( 'a*exp(-(x-mu)^2/2/sigma^2)+c', 'independent', 'x', 'dependent', 'y' );
    ft = fittype( 'a*exp(-(x-mu)^2/2/sigma^2)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    
    opts.StartPoint = [nanmax(yData) 0 sigmasTH(i)];
    %opts.StartPoint = [nanmax(yData) nanmean(meanBead(:)) 0 sigmasTH(i)];
    if beadSize == 0
        [fitresult, ~] = fit( xData, yData, ft, opts );
    else
        [fitresult, ~] = createConvolutedPsfFit(intVal(:,1+2*(i-1)), intVal(:,2*i), beadSize, sigmasTH(i), nanmax(yData),nanmin(yData));
    end
    % Plot fit with data.
    h_pca = plot( fitresult, xData, yData );
    set(h_pca(1),'Color',[i==1 i==2 i==3]);
    set(h_pca(2),'Color',[i==1 i==2 i==3]);
    legend('off');
    % Label axes
    xlabel x/µm
    ylabel intensity/a.u.
    grid on
    %save sigma
    sigmas(i)=abs(fitresult.sigma);
    %save maximum and minimum of xData and yData for goo axis scaling
    axisscaling(i,:)=[min(xData) max(xData) min(yData) max(yData)];
end
axis([min(axisscaling(:,1))/3 max(axisscaling(:,2))/3  min(axisscaling(:,3)) max(axisscaling(:,4))]);
FWHMs=sigmas*2*sqrt(2*log(2));
FWHM_pca = FWHMs;
title({['FWHM-PCA: ' '(' num2str(FWHMs(1)) ',']; [num2str(FWHMs(2)) ',' num2str(FWHMs(3)) ') µm'];['res_{lat,th}: ' num2str(res_lat/1000) ' µm'];['res_{ax,th}: ' num2str(res_ax/1000) ' µm']});

hold off;
%% fit slice through cog
A=[[1 0 0];[0 1 0];[0 0 1]];
%number of interpolated points in the trace (whole image)
steps=max(size(data));
intVal=zeros(steps,6);
for j=1:3
    %strongest component in U(:,j)
    [~,maxIdx]=max(abs(A(:,j)));
    tmin=-(cog(maxIdx)*scales(maxIdx)/A(maxIdx,j));
    tmax=size(meanBead,maxIdx)*scales(maxIdx)/A(maxIdx,j)+tmin;
    %interpolate data points on line
    intVal(:,2*j)=interp3(((1:size(meanBead,2))-cog(2))*scales(2),...
        ((1:size(meanBead,1))-cog(1))*scales(1),...
        ((1:size(meanBead,3))-cog(3))*scales(3),meanBead,...
        linspace(tmin*A(2,j),tmax*A(2,j),steps),...
        linspace(tmin*A(1,j),tmax*A(1,j),steps),...
        linspace(tmin*A(3,j),tmax*A(3,j),steps),...
        'cubic');
    %add scaling to intVal
    intVal(:,1+2*(j-1))=linspace(tmin,tmax,steps);
end
%fit
sigmasTH=[res_ax,res_lat,res_lat]*1.2/1000;
sigmas=zeros(1,3);
axisscaling=zeros(3,4);
ax(3)=subplot(2,3,6);
hold all
if beadSize == 0
    for i=1:3
        % Set up fittype and options.
        [xData, yData] = prepareCurveData( intVal(:,1+2*(i-1)), intVal(:,2*i) );
        %ft = fittype( 'a*exp(-(x-mu)^2/2/sigma^2)+c', 'independent', 'x', 'dependent', 'y' );
        ft = fittype( 'a*exp(-(x-mu)^2/2/sigma^2)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        %opts.StartPoint = [nanmax(yData) nanmean(meanBead(:)) 0 sigmasTH(i)];
        opts.StartPoint = [nanmax(yData) 0 sigmasTH(i)];
        
        [fitresult, ~] = fit( xData, yData, ft, opts);
        if(i==1)
            opts.StartPoint = [nanmax(yData) 0 sigmasTH(2)];
            [fitresult2, ~] = fit( xData, yData, ft, opts);
            confint1 = confint(fitresult);
            confint2 = confint(fitresult2);
            diff1 = confint1(1,1) - confint1(2,1);
            diff2 = confint2(1,1) - confint2(2,1);
            if abs(diff1)>abs(diff2)
                fitresult = fitresult2;
            end
        end
        % Plot fit with data.
        
        %     subplot(2,2,i);
        h_pca = plot( fitresult, xData, yData );
        set(h_pca(1),'Color',[i==3 i==1 i==2]);
        set(h_pca(2),'Color',[i==3 i==1 i==2]);
        %     legend( h_pca, 'Measurement Data', 'Fit', 'Location', 'NorthEast','Visible','False' );
        legend('off');
        % Label axes
        xlabel x/µm
        ylabel intensity/a.u.
        grid on
        %save sigma
        sigmas(i)=abs(fitresult.sigma);
        %save maximum and minimum of xData and yData for goo axis scaling
        axisscaling(i,:)=[min(xData) max(xData) min(yData) max(yData)];9
    end
else
    %beadSize is speified, use other fit method
    for i = 1:3
        [xData, yData] = prepareCurveData( intVal(:,1+2*(i-1)), intVal(:,2*i) );
        
        [fitresult, ~] = createConvolutedPsfFit(intVal(:,1+2*(i-1)), intVal(:,2*i), beadSize, sigmasTH(i), nanmax(yData), nanmin(yData));
        % Plot fit with data.
        
        %     subplot(2,2,i);
        h_pca = plot( fitresult, xData, yData );
        set(h_pca(1),'Color',[i==3 i==1 i==2]);
        set(h_pca(2),'Color',[i==3 i==1 i==2]);
        %     legend( h_pca, 'Measurement Data', 'Fit', 'Location', 'NorthEast','Visible','False' );
        legend('off');
        % Label axes
        xlabel x/µm
        ylabel intensity/a.u.
        grid on
        %save sigma
        sigmas(i)=abs(fitresult.sigma);
        %save maximum and minimum of xData and yData for goo axis scaling
        axisscaling(i,:)=[min(xData) max(xData) min(yData) max(yData)];
    end
end
axis([min(axisscaling(:,1))/3 max(axisscaling(:,2))/3  min(axisscaling(:,3)) max(axisscaling(:,4))]);
FWHMs=sigmas*2*sqrt(2*log(2));
FWHM_ortho = FWHMs;
title({['FWHM-Orth: ' '(' num2str(FWHMs(3)) ',']; [num2str(FWHMs(2)) ', ' num2str(FWHMs(1)) ') µm'];['res_{lat,th}: ' num2str(res_lat/1000) ' µm'];['res_{ax,th}: ' num2str(res_ax/1000) ' µm']});
hold off

% %% produce planar slices though pca axises
% %define planes: 1=z,2=y,3=x
%
% % xy-plane
% i=2;
% j=3;
% %strongest component in U(:,j)
% [~,maxIdxI]=max(abs(U(:,i)));
% tminI=-(cog(maxIdxI)*scales(maxIdxI)/U(maxIdxI,i));
% tmaxI=size(meanBead,maxIdxI)*scales(maxIdxI)/U(maxIdxI,i)+tminI;
%
% [~,maxIdxJ]=max(abs(U(:,j)));
% tminJ=-(cog(maxIdxJ)*scales(maxIdxJ)/U(maxIdxJ,j));
% tmaxJ=size(meanBead,maxIdxJ)*scales(maxIdxJ)/U(maxIdxJ,j)+tminJ;
%
%
% intVal=zeros(5*round((tmaxI-tminI)/metadata.scaleSize(1)),5*round((tmaxJ-tminJ)/metadata.scaleSize(2)));
%
% for k=1:size(intVal,1);
%
%     intVal(k,:)=interp3(((1:size(meanBead,2))-cog(2))*scales(2),...
%         ((1:size(meanBead,1))-cog(1))*scales(1),...
%         ((1:size(meanBead,3))-cog(3))*scales(3),meanBead,...
%         linspace((k*(tmaxI-tminI)/(size(intVal,1)-1)+tminI-(tmaxI-tminI)/(size(intVal,1)-1))*U(2,i)+(1*(tmaxJ-tminJ)/(size(intVal,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal,2)-1))*U(2,j),(k*(tmaxI-tminI)/(size(intVal,1)-1)+tminI-(tmaxI-tminI)/(size(intVal,1)-1))*U(2,i)+(size(intVal,2)*(tmaxJ-tminJ)/(size(intVal,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal,2)-1))*U(2,j),size(intVal,2)),...
%         linspace((k*(tmaxI-tminI)/(size(intVal,1)-1)+tminI-(tmaxI-tminI)/(size(intVal,1)-1))*U(1,i)+(1*(tmaxJ-tminJ)/(size(intVal,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal,2)-1))*U(1,j),(k*(tmaxI-tminI)/(size(intVal,1)-1)+tminI-(tmaxI-tminI)/(size(intVal,1)-1))*U(1,i)+(size(intVal,2)*(tmaxJ-tminJ)/(size(intVal,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal,2)-1))*U(1,j),size(intVal,2)),...
%         linspace((k*(tmaxI-tminI)/(size(intVal,1)-1)+tminI-(tmaxI-tminI)/(size(intVal,1)-1))*U(3,i)+(1*(tmaxJ-tminJ)/(size(intVal,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal,2)-1))*U(3,j),(k*(tmaxI-tminI)/(size(intVal,1)-1)+tminI-(tmaxI-tminI)/(size(intVal,1)-1))*U(3,i)+(size(intVal,2)*(tmaxJ-tminJ)/(size(intVal,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal,2)-1))*U(3,j),size(intVal,2)),...
%         'cubic');
% end
%
% % zx-plane
% i=1;
% j=2;
% %strongest component in U(:,j)
% [~,maxIdxI]=max(abs(U(:,i)));
% tminI=-(cog(maxIdxI)*scales(maxIdxI)/U(maxIdxI,i));
% tmaxI=size(meanBead,maxIdxI)*scales(maxIdxI)/U(maxIdxI,i)+tminI;
%
% [~,maxIdxJ]=max(abs(U(:,j)));
% tminJ=-(cog(maxIdxJ)*scales(maxIdxJ)/U(maxIdxJ,j));
% tmaxJ=size(meanBead,maxIdxJ)*scales(maxIdxJ)/U(maxIdxJ,j)+tminJ;
%
%
% intVal2=zeros(5*round((tmaxI-tminI)/metadata.scaleSize(3)),5*round((tmaxJ-tminJ)/metadata.scaleSize(1)));
%
% for k=1:size(intVal2,1);
%
%     intVal2(k,:)=interp3(((1:size(meanBead,2))-cog(2))*scales(2),...
%         ((1:size(meanBead,1))-cog(1))*scales(1),...
%         ((1:size(meanBead,3))-cog(3))*scales(3),meanBead,...
%         linspace((k*(tmaxI-tminI)/(size(intVal2,1)-1)+tminI-(tmaxI-tminI)/(size(intVal2,1)-1))*U(2,i)+(1*(tmaxJ-tminJ)/(size(intVal2,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal2,2)-1))*U(2,j),(k*(tmaxI-tminI)/(size(intVal2,1)-1)+tminI-(tmaxI-tminI)/(size(intVal2,1)-1))*U(2,i)+(size(intVal2,2)*(tmaxJ-tminJ)/(size(intVal2,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal2,2)-1))*U(2,j),size(intVal2,2)),...
%         linspace((k*(tmaxI-tminI)/(size(intVal2,1)-1)+tminI-(tmaxI-tminI)/(size(intVal2,1)-1))*U(1,i)+(1*(tmaxJ-tminJ)/(size(intVal2,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal2,2)-1))*U(1,j),(k*(tmaxI-tminI)/(size(intVal2,1)-1)+tminI-(tmaxI-tminI)/(size(intVal2,1)-1))*U(1,i)+(size(intVal2,2)*(tmaxJ-tminJ)/(size(intVal2,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal2,2)-1))*U(1,j),size(intVal2,2)),...
%         linspace((k*(tmaxI-tminI)/(size(intVal2,1)-1)+tminI-(tmaxI-tminI)/(size(intVal2,1)-1))*U(3,i)+(1*(tmaxJ-tminJ)/(size(intVal2,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal2,2)-1))*U(3,j),(k*(tmaxI-tminI)/(size(intVal2,1)-1)+tminI-(tmaxI-tminI)/(size(intVal2,1)-1))*U(3,i)+(size(intVal2,2)*(tmaxJ-tminJ)/(size(intVal2,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal2,2)-1))*U(3,j),size(intVal2,2)),...
%         'cubic');
% end
% % zy-plane
% i=1;
% j=3;
% %strongest component in U(:,j)
% [~,maxIdxI]=max(abs(U(:,i)));
% tminI=-(cog(maxIdxI)*scales(maxIdxI)/U(maxIdxI,i));
% tmaxI=size(meanBead,maxIdxI)*scales(maxIdxI)/U(maxIdxI,i)+tminI;
%
% [~,maxIdxJ]=max(abs(U(:,j)));
% tminJ=-(cog(maxIdxJ)*scales(maxIdxJ)/U(maxIdxJ,j));
% tmaxJ=size(meanBead,maxIdxJ)*scales(maxIdxJ)/U(maxIdxJ,j)+tminJ;
%
%
% intVal3=zeros(5*round((tmaxI-tminI)/metadata.scaleSize(3)),5*round((tmaxJ-tminJ)/metadata.scaleSize(2)));
%
% for k=1:size(intVal3,1);
%
%     intVal3(k,:)=interp3(((1:size(meanBead,2))-cog(2))*scales(2),...
%         ((1:size(meanBead,1))-cog(1))*scales(1),...
%         ((1:size(meanBead,3))-cog(3))*scales(3),meanBead,...
%         linspace((k*(tmaxI-tminI)/(size(intVal3,1)-1)+tminI-(tmaxI-tminI)/(size(intVal3,1)-1))*U(2,i)+(1*(tmaxJ-tminJ)/(size(intVal3,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal3,2)-1))*U(2,j),(k*(tmaxI-tminI)/(size(intVal3,1)-1)+tminI-(tmaxI-tminI)/(size(intVal3,1)-1))*U(2,i)+(size(intVal3,2)*(tmaxJ-tminJ)/(size(intVal3,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal3,2)-1))*U(2,j),size(intVal3,2)),...
%         linspace((k*(tmaxI-tminI)/(size(intVal3,1)-1)+tminI-(tmaxI-tminI)/(size(intVal3,1)-1))*U(1,i)+(1*(tmaxJ-tminJ)/(size(intVal3,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal3,2)-1))*U(1,j),(k*(tmaxI-tminI)/(size(intVal3,1)-1)+tminI-(tmaxI-tminI)/(size(intVal3,1)-1))*U(1,i)+(size(intVal3,2)*(tmaxJ-tminJ)/(size(intVal3,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal3,2)-1))*U(1,j),size(intVal3,2)),...
%         linspace((k*(tmaxI-tminI)/(size(intVal3,1)-1)+tminI-(tmaxI-tminI)/(size(intVal3,1)-1))*U(3,i)+(1*(tmaxJ-tminJ)/(size(intVal3,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal3,2)-1))*U(3,j),(k*(tmaxI-tminI)/(size(intVal3,1)-1)+tminI-(tmaxI-tminI)/(size(intVal3,1)-1))*U(3,i)+(size(intVal3,2)*(tmaxJ-tminJ)/(size(intVal3,2)-1)+tminJ-(tmaxJ-tminJ)/(size(intVal3,2)-1))*U(3,j),size(intVal3,2)),...
%         'cubic');
% end
% %% plot them
% ax(4)=subplot(2,3,3);
% hold on
% xy=contourf(...
%     linspace(((1)-cog(2))*scales(2),((size(meanBead,2))-cog(2))*scales(2),size(intVal,2)),...
%     linspace(((1)-cog(1))*scales(1),((size(meanBead,1))-cog(1))*scales(1),size(intVal,1)),...
%     intVal,linspace(nanmin(meanBead(:)),nanmax(meanBead(:)),10));
% xlabel x/µm
% ylabel y/µm
% legend('off');
% title XY-Slice
% axis equal;
% axis(5*[-sigmas(2) sigmas(2) -sigmas(1) sigmas(1)]);
% hold off
%
% ax(5)=subplot(2,3,1);
% hold on
% xz=contourf(...
%     linspace(((1)-cog(1))*scales(1),((size(meanBead,1))-cog(1))*scales(1),size(intVal2,2)),...
%     linspace(((1)-cog(3))*scales(3),((size(meanBead,3))-cog(3))*scales(3),size(intVal2,1)),...
%     intVal2,linspace(nanmin(meanBead(:)),nanmax(meanBead(:)),10));
% xlabel x/µm
% ylabel z/µm
% legend('off');
% title XZ-Slice
% axis equal;
% axis(5*[-sigmas(1) sigmas(1) -sigmas(3) sigmas(3)]);
% hold off
%
% ax(6)=subplot(2,3,2);
% hold on
% yz=contourf(...
%     linspace(((1)-cog(2))*scales(2),((size(meanBead,2))-cog(2))*scales(2),size(intVal2,2)),...
%     linspace(((1)-cog(3))*scales(3),((size(meanBead,3))-cog(3))*scales(3),size(intVal2,1)),...
%     intVal2,linspace(nanmin(meanBead(:)),nanmax(meanBead(:)),10));
% xlabel y/µm
% ylabel z/µm
% legend('off');
% title YZ-Slice
% axis equal;
% axis(5*[-sigmas(2) sigmas(2) -sigmas(3) sigmas(3)]);
% hold off
%% generate slices through main axis and center of gravity
zoomFactor = 2.5;

ax(4)=subplot(2,3,3);
hold on
pcolor(((1:size(meanBead,1))-cog(1))*scales(1),((1:size(meanBead,2))-cog(2))*scales(2),transpose(meanBead(:,:,round(cog(3)))));
shading interp;
contour(((1:size(meanBead,1))-cog(1))*scales(1),((1:size(meanBead,2))-cog(2))*scales(2),meanBead(:,:,round(cog(3)))',linspace(nanmin(meanBead(:)),nanmax(meanBead(:)),10),'Color','black');
xlabel x/µm
ylabel y/µm
legend('off');
title XY-Slice
axis equal;
if zoomIn == 1
    axis(zoomFactor*[-res_lat res_lat -res_lat res_lat]/1000);
end
hold off

ax(5)=subplot(2,3,1);
hold on
pcolor(((1:size(meanBead,1))-cog(1))*scales(1),((1:size(meanBead,3))-cog(3))*scales(3),transpose(squeeze(meanBead(:,round(cog(2)),:))));
shading interp;
contour(((1:size(meanBead,1))-cog(1))*scales(1),((1:size(meanBead,3))-cog(3))*scales(3),(squeeze(meanBead(:,round(cog(2)),:)))',linspace(nanmin(meanBead(:)),nanmax(meanBead(:)),10),'Color','black');
xlabel x/µm
ylabel z/µm
legend('off');
title XZ-Slice
axis equal;
if zoomIn == 1
    axis(zoomFactor*[-res_lat res_lat -res_ax res_ax]/1000);
end
hold off

ax(6)=subplot(2,3,2);
hold on
pcolor(((1:size(meanBead,2))-cog(2))*scales(2),((1:size(meanBead,3))-cog(3))*scales(3),transpose(squeeze(meanBead(round(cog(1)),:,:))));
shading interp;
contour(((1:size(meanBead,2))-cog(2))*scales(2),((1:size(meanBead,3))-cog(3))*scales(3),(squeeze(meanBead(round(cog(1)),:,:)))',linspace(nanmin(meanBead(:)),nanmax(meanBead(:)),10),'Color','black');
xlabel y/µm
ylabel z/µm
legend('off');
title YZ-Slice
axis equal;
if zoomIn == 1
    axis(zoomFactor*[-res_lat res_lat -res_ax res_ax]/1000);
end
hold off
%% maximum projection
maxProjection = max(data, [], 3);
minimum = double(min(maxProjection(:)));
maximum = double(max(maxProjection(:)));
maxProjection = uint8((double(maxProjection) - minimum)/(maximum-minimum)*255*2);
%make 3ch our of it
maxProjection = cat(3, maxProjection, maxProjection, maxProjection);
%draw red circles around all beaCenter
for i=1:size(beadCenter,2)
    maxProjection = insertShape(maxProjection,'circle',[beadCenter(2,i) beadCenter(1,i) sz(3)],'LineWidth',1, 'Color', 'red');
end
%draw red circles around all beaCenter
for i=1:size(beadCentersUsedBeads,2)
    pos = beadCentersUsedBeads{i};
    
    if ~isempty(pos)
        maxProjection = insertShape(maxProjection,'circle',[beadCentersUsedBeads{i}(2) beadCentersUsedBeads{i}(1) sz(3)],'LineWidth',2, 'Color', 'white');
    end
end
if showZProjection
    showZProjection = 'on';
else
    showZProjection = 'off';
end
figMaximumProjection=figure('Visible',showZProjection);
ax(8) = subplot(1,1,1);
imshow(maxProjection);
%overdraw them with white ones at beadCenterUsedBeads

%% create plot for microscopeEvaluation: Attempt to use a table (uitable) to display the values - hard to format properly
% fwhmquotient = FWHM_ortho./[res_lat, res_lat, res_ax]*1000;
% if showMicroscopeEvalPlot
%     temp = 'on';
% else
%     temp = 'off';
% end
% figOut = figure('visible', temp);
% [x,y,z] = ndgrid(linspace(1,size(meanBead,1),size(meanBead,1)),...
%                     linspace(1,size(meanBead,2),size(meanBead,2)),...
%                     linspace(1,size(meanBead,3), size(meanBead,3)/scales(2)*scales(3)));
% meanBeadI = interp3(meanBead, x,y,z, 'linear');
%
% outputImg = zeros(size(meanBeadI,1)+size(meanBeadI,3),size(meanBeadI,2)+size(meanBeadI,3));
% outputImg(1:size(meanBeadI,1),1:size(meanBeadI,2)) = max(meanBeadI,[],3);
% outputImg(size(meanBeadI,1)+1:size(meanBeadI,1)+size(meanBeadI,3),1:size(meanBeadI,2)) = squeeze(max(meanBeadI,[],2))';
% outputImg(1:size(meanBeadI,1),size(meanBeadI,2)+1:size(meanBeadI,2)+size(meanBeadI,3)) = squeeze(max(meanBeadI,[],1));
%
% panelH = uipanel('units','centimeters','Position',[0.5 0.5 10 10]);
% ax(7) = axes('Parent',panelH);
% imagesc(outputImg);
% axis image off
% set(gca,'unit','normalized','Position',[0 0 1 1])
% t = uitable(gcf,'Parent',panelH);
% t.Data = {num2str(res_lat,'%6.0f'),['<HTML><FONT color="red">' num2str(fwhmquotient(2),'%.3g') '</FONT></HTML>'],3;num2str(res_lat,'%6.0f'),num2str(fwhmquotient(2),'%.3g'),4;num2str(res_ax,'%6.0f'),num2str(fwhmquotient(3),'%.3g'),30;'0.0 °' 45 90};
% % t.ColumnName = {'<html><h6>theo</h6></html>','<html><h6>ortho</h6></html>','<html><h6>PCA</h6></html>'};
% t.ColumnName = {'theo','ortho','PCA'};
% t.RowName ={'x','y','z','ang'};
% axPos = get(gca,'Position');
% axRatio = [size(meanBeadI,1)/(size(meanBeadI,1)+size(meanBeadI,3)) size(meanBeadI,3)/(size(meanBeadI,1)+size(meanBeadI,3))];
% t.Units = 'normalized';
% t.Position = [axPos(1)+axPos(3)*axRatio(1) axPos(2) axPos(3)*axRatio(2) axPos(3)*axRatio(2)];
% t.FontSize = 7;
%% create plot for microscopeEvaluation
%interpolate meanBead in zDim
if showMicroscopeEvalPlot
    temp = 'on';
else
    temp = 'off';
end
figure('visible', temp);
ax(7) = subplot(1,1,1);

interpolationFactor = 5;
[x,y,z] = ndgrid(linspace(1,size(meanBead,1),interpolationFactor*size(meanBead,1)),...
    linspace(1,size(meanBead,2),interpolationFactor*size(meanBead,2)),...
    linspace(1,size(meanBead,3), interpolationFactor*size(meanBead,3)/scales(2)*scales(3)));
meanBeadI = interp3(meanBead, x,y,z, 'linear');
clear x y z

%truncate meanBeadI to interesting area
% pixelSize = scales(1)/interpolationFactor;
% FWHMs2Plot = 0.8;
% x = size(meanBeadI,1)/2-FWHMs2Plot*FWHM_ortho(1)/pixelSize:size(meanBeadI,1)/2+FWHMs2Plot*FWHM_ortho(1)/pixelSize;
% y = size(meanBeadI,2)/2-FWHMs2Plot*FWHM_ortho(2)/pixelSize:size(meanBeadI,2)/2+FWHMs2Plot*FWHM_ortho(2)/pixelSize;
% z = 1:size(meanBeadI,3); %size(meanBeadI,3)/2-FWHMs2Plot*FWHM_ortho(3)/pixelSize:size(meanBeadI,3)/2+FWHMs2Plot*FWHM_ortho(3)/pixelSize;
% x = round(x);
% y = round(y);
% z = round(z);
% meanBeadI = meanBeadI(x, x, z);

xWidth = size(meanBeadI,1);%2*round(interpolationFactor*8)+1;
%z-scaling factor
zWidth = size(meanBeadI,3);%2*round(interpolationFactor*scales(3)/scales(1)*8)+1; %means that we want the z-plots in the same resolution as xy plot -> we need to interpolate


%scale intensity values of meanBeadI
minimum = min(meanBeadI(:));
maximum = max(meanBeadI(:));
meanBeadI = (meanBeadI - minimum)/(maximum-minimum);

plotImage = ones(xWidth+zWidth, xWidth+zWidth, 3)*max(meanBeadI(:));
cogI = centerofgravity(meanBeadI);
cm = colormap('parula'); %create rgb colormap with 64 entries
%copy xyData
plotImage(1:xWidth, 1:xWidth, :) = grs2rgb(meanBeadI(1:end, 1:end, floor(cogI(3))), cm); %grs2rgb(meanBeadI(interpolationFactor*cog(1)-(xWidth-1)/2:interpolationFactor*cog(1)+(xWidth-1)/2,...
%interpolationFactor*cog(2)-(xWidth-1)/2:interpolationFactor*cog(2)+(xWidth-1)/2,...
%round(interpolationFactor*cog(3)*scales(3)/scales(1))),cm);
%copy xzData
plotImage(1:xWidth, xWidth+1:end, :) = grs2rgb(meanBeadI(1:end, floor(cogI(2)), 1:end), cm);%grs2rgb(meanBeadI(interpolationFactor*cog(1)-(xWidth-1)/2:interpolationFactor*cog(1)+(xWidth-1)/2,...
%interpolationFactor*cog(2),...
%round(interpolationFactor*cog(3)*scales(3)/scales(1))-(zWidth-1)/2:round(interpolationFactor*cog(3)*scales(3)/scales(1))+(zWidth-1)/2),cm);
%copy yzData
plotImage(xWidth+1:end, 1:xWidth, :) = grs2rgb(squeeze(meanBeadI(floor(cogI(2)), 1:end, 1:end))', cm);%grs2rgb(permute(squeeze(meanBeadI(interpolationFactor*cog(1),...
%interpolationFactor*cog(2)-(xWidth-1)/2:interpolationFactor*cog(2)+(xWidth-1)/2,...
%round(interpolationFactor*cog(3)*scales(3)/scales(1))-(zWidth-1)/2:round(interpolationFactor*cog(3)*scales(3)/scales(1))+(zWidth-1)/2)),[2,1]),cm);

availableWidth = zWidth;
spaceToBorders = 0.01*availableWidth;
spaceToBorders = round(spaceToBorders);
%we have 5 rows and 5 columns
colPos = round(xWidth + spaceToBorders + (availableWidth-2*spaceToBorders)/4*(0:3));
rowPos = round(xWidth + spaceToBorders + (availableWidth-2*spaceToBorders)/5*(0:4));

fontSize = (colPos(2)-colPos(1) + 6.8445)/3.2;%3 should be normally be 2.,6247, cjhose 3 for smaller text size and better fitting%(colPos(2)-colPos(1)-2*spaceToBorders+3.9692)/16.507; %neededWidth = 1.985 * fontSize + 0.8418 -> fontsize = (neededWidth - 0.8418)/1.985
fontSize=round(fontSize);

%headlines:
%theo:
plotImage = insertText(plotImage, [colPos(2) rowPos(1)], 'theo', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%ortho
plotImage = insertText(plotImage, [colPos(3) rowPos(1)], 'ortho', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%pca
plotImage = insertText(plotImage, [colPos(4) rowPos(1)], 'pca', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%xy:
plotImage = insertText(plotImage, [colPos(1) rowPos(2)], 'x', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%xz
plotImage = insertText(plotImage, [colPos(1) rowPos(3)], 'y', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%yz
plotImage = insertText(plotImage, [colPos(1) rowPos(4)], 'z', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%angle
plotImage = insertText(plotImage, [colPos(1) rowPos(5)], 'ang', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);

%fill with values
%x_theo:
plotImage = insertText(plotImage, [colPos(2) rowPos(2)], num2str(res_lat,'%.4g'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%y_theo:
plotImage = insertText(plotImage, [colPos(2) rowPos(3)], num2str(res_lat,'%.4g'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
%z_theo:
plotImage = insertText(plotImage, [colPos(2) rowPos(4)], num2str(res_ax,'%.4g'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);

%angle_theo:
plotImage = insertText(plotImage, [colPos(2) rowPos(5)], '0.0 °', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);
fwhmquotient = FWHM_ortho./[res_lat, res_lat, res_ax]*1000;
%x_ortho:
if fwhmquotient(1) > redThreshold
    color = 'red';
else
    color = 'black';
end
plotImage = insertText(plotImage, [colPos(3) rowPos(2)], num2str(fwhmquotient(1),'%.2f'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', color, 'BoxOpacity', 0);
%y_ortho:
if fwhmquotient(2) > redThreshold
    color = 'red';
else
    color = 'black';
end
plotImage = insertText(plotImage, [colPos(3) rowPos(3)], num2str(fwhmquotient(2),'%.2f'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', color, 'BoxOpacity', 0);
%z_ortho:
if fwhmquotient(3) > redThreshold
    color = 'red';
else
    color = 'black';
end
plotImage = insertText(plotImage, [colPos(3) rowPos(4)], num2str(fwhmquotient(3),'%.2f'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', color, 'BoxOpacity', 0);
%angle_ortho:
plotImage = insertText(plotImage, [colPos(3) rowPos(5)], '0.0 °', 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);

fwhmquotient = [FWHM_pca(2), FWHM_pca(3), FWHM_pca(1)];
fwhmquotient = fwhmquotient./[res_lat, res_lat, res_ax]*1000;
%x_pca:
if fwhmquotient(1) > redThreshold
    color = 'red';
else
    color = 'black';
end
plotImage = insertText(plotImage, [colPos(4) rowPos(2)], num2str(fwhmquotient(1),'%.2f'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', color, 'BoxOpacity', 0);
%y_pca:
if fwhmquotient(2) > redThreshold
    color = 'red';
else
    color = 'black';
end
plotImage = insertText(plotImage, [colPos(4) rowPos(3)], num2str(fwhmquotient(2),'%.2f'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', color, 'BoxOpacity', 0);
%z_pca:
if fwhmquotient(3) > redThreshold
    color = 'red';
else
    color = 'black';
end
plotImage = insertText(plotImage, [colPos(4) rowPos(4)], num2str(fwhmquotient(3),'%.2f'), 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', color, 'BoxOpacity', 0);
%angle_pca:
plotImage = insertText(plotImage, [colPos(4) rowPos(5)], [num2str(acos((sum(U(:,1).*[0;0;1]))/(sqrt(sum(U(:,1).^2))))*180/pi, '%.2f') ' °'], 'FontSize', fontSize, 'Font', 'Arial Bold', 'TextColor', 'black', 'BoxOpacity', 0);

%calculate fontSize
% availableWidth = zWidth;
% spaceToBorders = 0.05 * availableWidth;
%
%
% %round all values
% availableWidth = round(availableWidth);
% fontSize = round(fontSize);
% spaceToBorders = round(spaceToBorders);
%
% %insert text
% plotImage = insertText(plotImage, [xWidth + spaceToBorders, xWidth + zWidth/5], ['FWHM_theo = (' num2str(round(res_lat)) ', ' num2str(round(res_lat)) ', ' num2str(round(res_ax)) ') nm'], 'FontSize', fontSize, 'BoxOpacity', 0, 'TextColor', 'black', 'Font', 'Arial Bold');
% fwhmquotient = FWHM_ortho./[res_lat, res_lat, res_ax]*1000;
% %insert quotients, write them red if redThreshold is exceeded
% if any(round(fwhmquotient,2)> redThreshold)
%     textColor = 'red';
% else
%     textColor ='black';
% end
% plotImage = insertText(plotImage, [xWidth + spaceToBorders, xWidth + 2 * zWidth/5], ['ortho/theo = (' num2str(round(fwhmquotient(1),2)) ', ' num2str(round(fwhmquotient(2),2)) ', ' num2str(round(fwhmquotient(3),2)) ')'], 'FontSize', fontSize,...
%     'BoxOpacity', 0, 'TextColor', textColor, 'Font', 'Arial Bold');
% fwhmquotient = [FWHM_pca(2), FWHM_pca(3), FWHM_pca(1)];
% fwhmquotient = fwhmquotient./[res_lat, res_lat, res_ax]*1000;
% if any(round(fwhmquotient,2)> redThreshold)
%     textColor = 'red';
% else
%     textColor ='black';
% end
% plotImage = insertText(plotImage, [xWidth + spaceToBorders, xWidth + 3*zWidth/5], ['pca/theo = (' num2str(round(fwhmquotient(1),2)) ', ' num2str(round(fwhmquotient(2),2)) ', ' num2str(round(fwhmquotient(3),2)) ')'], 'FontSize', fontSize,...
%    'BoxOpacity', 0, 'TextColor', textColor, 'Font', 'Arial Bold');
% plotImage = insertText(plotImage, [xWidth + spaceToBorders, xWidth + 4*zWidth/5], ['ang = ' num2str(acos((sum(U(:,1).*[0;0;1]))/(sqrt(sum(U(:,1).^2))))*180/pi), ' °'], 'FontSize', fontSize, 'BoxOpacity', 0, 'TextColor', 'black', 'Font', 'Arial Bold');


%plot it
image(plotImage(end:-1:1,:,:));
% pcolor(linspace(0,size(plotImage,1)*scales(1)/interpolationFactor, size(plotImage,1)), linspace(0,size(plotImage,1)*scales(1)/interpolationFactor, size(plotImage,1)), plotImage);
hold on
% set(gca, 'Ydir', 'reverse');
% shading interp;
legend('off');
axis equal
% xlabel('µm');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'YDir', 'normal');
axis([0 size(plotImage, 1) 0 size(plotImage,2)]);
% close waitbar
close(wbar);
% output plot handles h

h{1}=fig;
for i=1:length(ax)
    h{i+1}=findobj(ax(i));
end
%calculate quotient of FWHM measurement to theoretical value
h{length(ax)+2} = FWHM_ortho./[res_lat, res_lat, res_ax]*1000;
FWHM_pca = [FWHM_pca(2), FWHM_pca(3), FWHM_pca(1)];
h{length(ax)+3} = FWHM_pca./[res_lat, res_lat, res_ax]*1000;
h{length(ax)+4} = [res_lat, res_lat, res_ax];
%% generate output figure
figSize = [21 29.7];  % figure size in cm (A4)
if showPlot == 1
    showPlot = 'on';
else
    showPlot = 'off';
end
outputfig = figure('PaperOrientation','Portrait', ...
    'PaperType', 'A4',...
    'PaperUnits', 'centimeter',...
    'PaperSize', figSize,...
    'Units','centimeter',...
    'Position',[0.1 0.1 figSize],...
    'Color','w','Visible',showPlot);
for i=1:6
    ax2(i)=copyobj(h{i+1}(1),outputfig);
end

end