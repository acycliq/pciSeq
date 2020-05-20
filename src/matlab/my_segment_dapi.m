function [CellMap, Boundaries] = my_segment_dapi(BigDapiFile)
% [0 CellMap BoundaryImage] = o.segment_dapi(DapiIm)
%
% segments a DAPI image and assigns each pixel to a cell. Input the path to BigDapiFile
%
% only works within region outlined by roi
%
% Output CellMap is same size as input DapiIm, with integer entries for
% each pixel, saying which cell it belongs to. (Zero if unassigned)
%
% also saved to o.CellMapFile


%% find Cell Calling Region
% y0 = min(o.CellCallRegionYX(:,1));
% x0 = min(o.CellCallRegionYX(:,2));
% y1 = max(o.CellCallRegionYX(:,1));
% x1 = max(o.CellCallRegionYX(:,2));


% % roi is [x0 x1 y0 y1]
% x0 = roi(1);
% x1 = roi(2);
% y0 = roi(3);
% y1 = roi(4);

% set the parameters
DapiThresh = 90;
DapiMinSize = 5;
DapiMinSep = 7;
DapiMargin = 10;
MinCellArea = 200;

% Mask = poly2mask(o.CellCallRegionYX(:,2)-x0+1, o.CellCallRegionYX(:,1)-y0+1, y1-y0+1, x1-x0+1);
% Dapi = imread(o.BigDapiFile, 'PixelRegion', {[y0 y1], [x0 x1]}).*uint16(Mask);
% Dapi = imread(BigDapiFile);

%%
Dapi = BigDapiFile;
if max(size(size(Dapi))) == 3
    % if dapi is a 3-dim array (m-by-n-by-k) then it is a coloured image.
    % Convert to grayscale
    Dapi = rgb2gray(Dapi(:,:,1:3));
end    
Dapi = imadjust(Dapi); % contrast enhancement
ImSz = size(Dapi);
Debug = 1;
%% threshold the map
% ThreshVal = prctile(Dapi(Mask), DapiThresh);
% ThreshVal = 15948; % I added that!
ThreshVal = prctile(Dapi(Dapi>0), DapiThresh);
bwDapi = imerode(Dapi>ThreshVal, strel('disk', 2));

if Debug
    figure(300)
    subplot(2,1,1);
    imagesc(Dapi); 
    subplot(2,1,2);
    imagesc(bwDapi);
    colormap bone
    fprintf('Threshold = %f\n', ThreshVal);
    
end
%% find local maxima 
dist = bwdist(~bwDapi);
dist0 = dist;
dist0(dist<DapiMinSize)=0;
ddist = imdilate(dist0, strel('disk', DapiMinSep));
%clear dist 
impim = imimposemin(-dist0, imregionalmax(ddist));
clear dist0
if Debug
    figure(301);
    subplot(2,1,1)
    imagesc(dist);
    subplot(2,1,2)
    imagesc(impim);
end
%% segment
% remove pixels at watershed boundaries
bwDapi0 = bwDapi;
bwDapi0(watershed(impim)==0)=0;

% assign all pixels a label
labels = uint32(bwlabel(bwDapi0));
[d, idx] = bwdist(bwDapi0);

% now expand the regions by a margin
CellMap0 = zeros(ImSz, 'uint32');
Expansions = (d<DapiMargin);
CellMap0(Expansions) = labels(idx(Expansions));

% get rid of cells that are too small
rProps0 = regionprops(CellMap0); % returns XY coordinate and area
BigEnough = [rProps0.Area]>MinCellArea;
NewNumber = zeros(length(rProps0),1);
NewNumber(~BigEnough) = 0;
NewNumber(BigEnough) = 1:sum(BigEnough);
CellMap = CellMap0;
CellMap(CellMap0>0) = NewNumber(CellMap0(CellMap0>0));

if Debug
    figure(302)
    image(label2rgb(CellMap, 'jet', 'w', 'shuffle'));

end

%%
% CellYX = fliplr(vertcat(rProps(BigEnough).Centroid)); % because XY

%% make image with boundaries
Boundaries = (CellMap ~= imdilate(CellMap,strel('disk', 1)));
% DapiBoundaries = Dapi;
% 
% OddPix = mod((1:size(Dapi,2)) + (1:size(Dapi,1))', 2);
% DapiBoundaries(Boundaries & OddPix) = .3 * max(Dapi(:));
% DapiBoundaries(Boundaries & ~OddPix) = 0;
% 
% imwrite(DapiBoundaries, 'my_segmentation.tif');

% o.CellMapFile = fullfile(o.OutputDirectory, 'CellMap.mat');
% save(o.CellMapFile, 'CellMap', 'DapiBoundaries', 'y0', 'y1', 'x0', 'x1');


end

