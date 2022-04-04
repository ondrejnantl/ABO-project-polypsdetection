function [binaryMap] = detectPolyps(inputImage)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Smazani odlesku - alternativa 2
pmRefl = rangefilt(inputImage,true(7));
reflMask = imbinarize(imfill(pmRefl,'holes'));
noReflImage = regionfill(inputImage,reflMask);

% Detekce hran v obraze
edgedImage = edge(noReflImage,'canny',std(noReflImage(:)));
% [edgeRows,edgeCols] = find(edgedImage == 1);
% labelImage = bwlabel(edgedImage);
% for i = 1:length(edgeRows)
%     for j = 1:length(edgeCols)
%         if sum(edgedImage((edgeRows(i)-1):(edgeRows(i)+1),(edgeCols(i)-1):(edgeCols(i)+1)),'all') < 3
%             edgedImage((edgeRows(i)-1):(edgeRows(i)+1),(edgeCols(i)-1):(edgeCols(i)+1)) = 0;
%         elseif sum(edgedImage((edgeRows(i)-1):(edgeRows(i)+1),(edgeCols(i)-1):(edgeCols(i)+1)),'all') > 6
%             edgedImage((edgeRows(i)-1):(edgeRows(i)+1),(edgeCols(i)-1):(edgeCols(i)+1)) = 1;
%         end
%     end
% end

end

