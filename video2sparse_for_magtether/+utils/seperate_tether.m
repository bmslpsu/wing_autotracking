function img_tether_dilated = seperate_body_tether(Image_bw, pts_teth, bottomFlag)
%% Seperates body & tether masks
%
% INPUT:
%   Image_bw: binarized image or body and tether mask
%   bottomFlag: is this the bottom view
%   pts_teth: location of tether
%
% OUTPUT
%   img_body: body image or body mask
%   img_tether: tether image or mask
%

% Extract body and tether mask
y = round(pts_teth(:,2)); % y and x location of tether points
if ~bottomFlag % side view
    img_tether = Image_bw;
    img_tether(y(1):end,:) = 0;

elseif bottomFlag % bottom view which needs special methods
    Image_bw(y(2):end,:) = 0; % removes the brass rod(not particularly needed now)
    Image_bw = bwareafilt(Image_bw,3); % removes the noise and small objects
    Image_bw_noBody = Image_bw;
    Image_bw_noBody(1:y(1),:) = 0;
%     Vector_sum = sum(Image_bw_noBody,2);
%     index_sepPoint = find(Vector_sum > 10, 1, 'first'); % finds the first point at which the body enters the tether mask
%     if ~isempty(index_sepPoint)
%         Image_bw_noBody(1:index_sepPoint,:) = 0; % removes the body
%     end
    img_tether = Image_bw_noBody;
end
se = strel('square',3);
img_tether_dilated = imdilate(img_tether,se);

end