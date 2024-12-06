%% Code written by Anne Lyons 2023 11 17 with modifications by Michelle Frei
%% Dependencies: Install the Image Processing Toolbox
%% Load images: change paths + names of files 
sted0 = imread("STED\2023_12_08_136K_JF635_D2_well1_STED.lif - Series003 - C=1.tif");
sted1 = imread("STED\2023_12_08_136K_JF635_D2_well1_STED.lif - Series004 - C=1.tif");

%Resolution from image meta data 
imres = 640/22.56;

%% Segment Images
% Functions are at end of script
sted0_obj = SegmentCLC(sted0);
sted1_obj = SegmentCLC(sted1);

% Label objects to give IDs 
sted0_Lb = bwlabel(sted0_obj);
sted1_Lb = bwlabel(sted1_obj);

sted0_mask = boundarymask(sted0_Lb);
sted1_mask = boundarymask(sted1_Lb);

newFolderName = 'Results_2';
mkdir(newFolderName);

figure('Name','STED before');
imshow(labeloverlay(imadjust(sted0),sted0_mask,'Transparency',0));
saveas(gcf, 'Results_2\STED-before_ROI.tiff', 'tiff');
figure('Name','STED after');
imshow(labeloverlay(imadjust(sted1),sted1_mask,'Transparency',0));
saveas(gcf, 'Results_2\STED-after_ROI.tiff', 'tiff');

%% Distance measures to partner pits
[stedConn,stedCent0,stedCent1] = LinkCLC(sted0,sted1,sted0_obj,sted1_obj);
title('STED Pair');
saveas(gcf, 'Results_2\STED pair', 'tiff');

%% Define background corrected images
imshow(sted1);
imcontrast;
h = imrect;
roiMask = createMask(h);
imwrite(roiMask, 'Results_2\ROI.tiff');
%% Define the same background ROI if it was drawn previously
%roiMask = imread('Results_2\ROI.tiff');
%% BG corr
bg_s0= mean(sted0(roiMask));
bg_s1= mean(sted1(roiMask));

sted0_bg = double(sted0)-bg_s0;
sted1_bg = double(sted1)-bg_s1;


%% F1/F0 measures of partner pits

for conn = 1:length(stedConn)
    rawim0 = sted0_bg;
    % set everything where object label is not true == 0
    % 0 before image is LocB in Column2
    rawim0(sted0_Lb ~= stedConn(conn,2)) = 0;
    % get all nonzero intensity values
    F0sted{conn} = double(nonzeros(rawim0(:)));
    F0_mean_sted(conn,:) = mean(F0sted{conn});

    rawim1 = sted1_bg;
    % set everything where connected label is not true == 0
    % 1 afterimage is LocObj in Column1
    rawim1(sted1_Lb ~= stedConn(conn,1)) = 0;
    % get all nonzero intensity values
    F1sted{conn} = double(nonzeros(rawim1(:)));
    F1_mean_sted(conn,:) = mean(F1sted{conn});
end
Ratio_sted = F1_mean_sted./F0_mean_sted;
deltaF_sted = F1_mean_sted./F0_mean_sted-1;


stedConn = [stedConn,F0_mean_sted,F1_mean_sted,Ratio_sted];
%Filter the files for Distances below sqrt(200) pixels = 14 pixels = 500 nm  
stedConn_filt = stedConn(stedConn(:,3)<200,:);
% so stedConn have the follwing columns: ID after, ID before,
% Distance, Mean before, mean after, ratio

%% Colour coded ratio image for all pits not filtered 
ratioim = zeros(length(sted1),length(sted1)); 
for conn2 = 1:length(stedConn)
    rawim_sted1 = sted1_obj;
    % set everything where object label is not true == 0
    % 0 before image is LocB in Column2
    rawim_sted1(sted1_Lb ~= stedConn(conn2,1)) = 0;
    % get all nonzero intensity values
    rawim_sted1_nw = rawim_sted1*stedConn(conn2,6);
    ratioim = ratioim + rawim_sted1_nw; 
end 
ratioim = ratioim*1000;

figure;
imagesc(ratioim);
colormap ("jet");
colorbar;
title('Ratioimage');
saveas(gcf, 'Results_2\Ratio_image.tiff', 'tiff');
ratioim_16 = uint16(ratioim);
imwrite(ratioim_16, 'Results_2\Ratio_image_16_scaled1000.tiff');

%% Output Files and Plotting
%make a figure correlatinge ratio and distance 
figure;
scatter(stedConn_filt(:,3),stedConn_filt(:,6),'blue')
ylabel('Mean F1/F0'); xlabel('Pit Partner Distance (Pixels)');
legend('STED');
saveas(gcf, 'Results_2\Distance-Ratio.tiff', 'tiff');

figure; 
histogram(stedConn_filt(:,6), 'FaceColor','blue','BinWidth', 0.1)
xlabel('Mean F1/F0'); ylabel('Frequency');
legend('STED');
saveas(gcf, 'Results_2\RatioHisto', 'tiff');

stedT = array2table(stedConn_filt,...
    'VariableNames',{'Im1_ObjID','Im0_ObjID','ObjDist','MeanF0','MeanF1','Mean Ratio'});
writetable(stedT, 'Results_2\STED-Anlaysis.xlsx');
close all;
clear all;

%% Subroutine functions

function CLCobj = SegmentCLC(input_image)
% A funtion with difference of gaussians filtering to enhance dot
% recognition in segmentation

    %Filter settings determined in DotSpotter/ROInavigator
    %Segmentation Parameters 
    hsize = 5;
    sigma = 5;
    diffOG = 5;
    blurSig = 2;

    % Filtering parameters 
    CircMin = 0.2;
    CircMax = 1;
    AreaMin = 25;
    % Futher parameter
    %AreaMax = 10181;
    
    
    % Difference of Gaussians to search enhance gaussian-like puncta
    gaussian1 = fspecial('Gaussian',hsize,sigma);
    gaussian2 = fspecial('Gaussian',hsize,sigma+diffOG);
    
    %Defintion of the Kernel for DOG 
    DoG = gaussian1-gaussian2;
    % search for the same shape as DoG filter 
    DoG_filt = conv2(double(input_image),DoG,'same');

    % Noise reduction with gaussian blur sigma
    DoG_blur = imgaussfilt(DoG_filt,blurSig);
    
    % Segment filtered image
    J = imbinarize(DoG_blur);
    JL = bwlabel(J);
    stats = regionprops(J,"Circularity","Area");
   
    keepMask1 = [stats.Circularity]>=CircMin;
    keepMask2= [stats.Circularity]<=CircMax;
    keepMask3 = [stats.Area]>=AreaMin;
    keep = keepMask1 + keepMask2 + keepMask3;
    CLCobj = ismember(JL,find(keep == 3)) > 0;
    
    % Include if AreaMax is used to filter 
    % keepMask4 = [stats.Area]<=AreaMax;
    % keep = keepMask1 + keepMask2 + keepMask3 + keepMask4;
    % CLCobj = ismember(JL,find(keep == 4)) > 0;

end

function [fcheck2,puncCent0,puncCent1] = LinkCLC(im0,im1,input_binary0,input_binary1)
% A funtion to link segmented objects between two images that are the
% closest together possible, perhaps two timepoint images
    figure;
    imshowpair(im0,im1); 
    text (100,670, 'Green - before, Pink - after, yellow - all pairings, red- filtered pairing');
    hold on;
    
    stats0 = regionprops('table',input_binary0,'Centroid');
    puncCent0 = round(stats0.Centroid);
    stats1 = regionprops('table',input_binary1,'Centroid');
    puncCent1 = round(stats1.Centroid);
    % For object ID overlay onto image, uncomment for loops below:
    % for i = 1:length(puncCent0)
    %     text(puncCent0(i,1),puncCent0(i,2),num2str(i),'Color','g','FontSize',8);
    % end
    % for i = 1:length(puncCent1)
    %     text(puncCent1(i,1),puncCent1(i,2),num2str(i),'Color','m','FontSize',8);
    % end
 
    
    % Measure distances from after image (1) to before image (0)... the number
    % of measures is equal to the number of objects in the after image (1)
    [punc2punc,LocObj,LocB] = EucDistObj(puncCent1,puncCent0);

    % LocObj is the list of object IDs in afterimage(1) (ordered 1:object#)
    % LocB is the list of objectIDs in beforeimage(0) index matched to LocObj
    
    % While LocObj has no duplicates and is ordered numerically, it is 
    % possible for LocB to have duplicates... more on dealing with this
    % below
   
    % plots the distances over the overlay image in 
    for pln = 1:length(punc2punc(:,1)) 
        p3 = plot([puncCent1(LocObj(pln),1) puncCent0(LocB(pln),1)], ...
        [puncCent1(LocObj(pln),2) puncCent0(LocB(pln),2)],'y',...
        'MarkerSize',3,'DisplayName','punc2punc'); hold on;
    end
    
    % in the initial approach, the partnering errors where
    % there are pits with several partners (duplicates in LocB)
    % let's fix this by only keeping pairing with the closest neighbor
    
    % Part A
    %new matrix with the after, before, distance
    check = [LocObj,LocB,punc2punc];
    %Sort by column2 (LocB) to filter multipartners
    check2 = sortrows(check,2); 
    fcheck2 = [];
    % this routine goes through all the indices in the ordered matrix
    for pln = 1:max(LocB(:,1))
        ind = find(check2(:,2)==pln);
        if length(ind) == 1 %if only one partner exists, it adds this to the fcheck2 version
            check_app = check2(ind,:);
            if isempty(fcheck2) == 1
                fcheck2 = check_app;
            else
                fcheck2 = [fcheck2; check_app];
            end
        elseif length(ind) > 1 % if multiple indices in LocB contain same value
            %this is the subsection of check2 that has the same indicies
            %(also contains the Lobj, Bobj, distance
            ddist = check2((min(ind):max(ind)),:);
            [~,I] = min(ddist(:,3)); %get indexed minimum, the row that is the shortest
            check_app = ddist(I,:); %and use for appending matrix
            if isempty(fcheck2) == 1
                fcheck2 = check_app;
            else
                fcheck2 = [fcheck2; check_app];
            end
        end
    end
    
    % define the new selected pairings
    LocObjf = fcheck2(:,1); 
    LocBf= fcheck2(:,2); 
    punc2puncf= fcheck2(:,3);

    for pln = 1:length(LocBf(:,1))
        p3 = plot([puncCent1(LocObjf(pln),1) puncCent0(LocBf(pln),1)], ...
        [puncCent1(LocObjf(pln),2) puncCent0(LocBf(pln),2)],'c',...
        'MarkerSize',3,'DisplayName','punc2punc', 'Color', 'red'); hold on;
    end
    
    hold off

end

function [minEdist,LocObj,LocBoundary] = EucDistObj(Object,Boundary,varargin)
% A function to take the squared Euclidean distance of an "object" relative 
% to a "boundary", which are These can just be two coordinate systems
% Obj2Obj varargin option allows measurements within the same image 
% without having objects measure to themself.
% the variables returned are the minimum distances minEdist (column vector
% with the minima, LocObj the indicies of the objects in the after image
% (1,2,3,4) Loc Boundary indicies in the before image

% matrix with all the distances between the objects 
Edist = dist2(Object,Boundary);

if nargin > 2
    if strcmp(varargin{1},'obj2obj') == 1 % Distance between objects
        %if we use the same method btwn objects, they measure to themselves
        %making a zero value minimum... instead, sort the columns of Edist to
        %numerically get the minimums in the second (nonzero) row
        sortEdist = sort(Edist);
        minEdist = sortEdist(2,:)';
    end
else %distance from object to boundary
    minEdist = min(Edist,[],2);
end

LocBoundary = zeros(length(Edist(:,1)),1); 
LocObj = LocBoundary;

% LocB will have the indices from Edist corresponding to the minimal
% function  hence that is the object in the before image
%Loc Boundary will only store the first value (at the row)
%Loc object will save the indices of the current object
for i = 1:length(Edist(:,1))
    LocB = find(minEdist(i,:)==Edist(i,:));
    LocBoundary(i,:) = LocB(1); %in case multiple elements are found
    LocObj(i,:) = i;
end

end

function n2 = dist2(x, c)
% DIST2	Calculates squared distance between two sets of points.
%
%	Description
%	D = DIST2(X, C) takes two matrices of vectors and calculates the
%	squared Euclidean distance between them.  Both matrices must be of
%	the same column dimension.  If X has M rows and N columns, and C has
%	L rows and N columns, then the result has M rows and L columns.  The
%	I, Jth entry is the  squared distance from the Ith row of X to the
%	Jth row of C.
%
%	See also
%	GMMACTIV, KMEANS, RBFFWD
%
%
%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

% checks how big my two lists of centroids are 
[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
	error('Data dimension does not match dimension of centres')
end

n2 = (ones(ncentres, 1) * sum((x.^2)', 1))' + ...
  		ones(ndata, 1) * sum((c.^2)',1) - ...
  		2.*(x*(c'));
end
