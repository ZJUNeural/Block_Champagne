addpath(genpath('./'));
%% Load DBS data
subj = 'sub-01'; run = 'run-01';
load([subj,'_',run,'_StimulatedData.mat']);
load([subj,'_',run,'_Forward.mat']);
eT = [-0.01,0.01];
[~,snaps] = min(abs(epoTime-eT(1)));
[~,snape] = min(abs(epoTime-eT(2)));

%% Select Good Channel
badlist = [];
for it = 1:size(epochInfo,1)
   tmp = epochInfo(it,:);
   tmp = str2double(tmp(2:end));
   badlist = [badlist,tmp]; 
end
goodlist = setdiff(1:size(avgData,1),badlist);
L = double(Gain_fix(goodlist,:));
B = avgData(goodlist,snaps:snape)*1e6;
%% Get Downsampled Cortex
% the index start from 0 in python and 1 in matlab, so the Index should +1
Leftcv = denseCortexVertices{1}(CortexVerticesIdx(1,:)+1,:);
Rightcv = denseCortexVertices{2}(CortexVerticesIdx(2,:)+1,:);

% transform the face index from dense cortex to downsampled cortex 
LeftTri = double(squeeze(CortexFaces(1,:,:)));
for it = 1:numel(CortexVerticesIdx(1,:))
    vtIdx = CortexVerticesIdx(1,it);
    LeftTri(LeftTri==vtIdx) = it;    
end
RightTri = double(squeeze(CortexFaces(2,:,:)));
for it = 1:numel(CortexVerticesIdx(2,:))
    vtIdx = CortexVerticesIdx(2,it);
    RightTri(RightTri==vtIdx) = it;    
end
%% Merge Cortex
cv = [Leftcv;Rightcv];
LeftIdx = 1:size(Leftcv,1);
RightIdx = size(Leftcv,1)+(1:size(Rightcv,1));
tri = [LeftTri;RightTri+size(Leftcv,1)];
Cortex_low.Vertices = cv+trans(1:3,4)';
Cortex_low.Faces = tri;
Cortex_low.VertConn = tess_vertconn(cv,tri);
Scouts(1).Vertices = LeftIdx; Scouts(2).Vertices = RightIdx;
Cortex_low.Atlas(1).Name='Structures';
Cortex_low.Atlas(1).Scouts = Scouts;

[nSensor,nSource] = size(L);

%% ESI
J = BlockChampagne(L,B,Cortex_low.VertConn,[],'knb',1);
[~,maxT] = max(std(B,[],1));
plotStim_source(abs(J(1:nSource,maxT)),Cortex_low,stim_coords,0.5);view([90,10]);



function plotStim_source(source,Cortex,coords,thresh)
if nargin<4
    thresh = [];
end
visual_source(source,Cortex,thresh);hold on; 
scatter3(mean(coords(:,1)),mean(coords(:,2)),mean(coords(:,3)),'filled');

scatter3((coords(1,1)),(coords(1,2)),(coords(1,3)),'filled');
scatter3((coords(2,1)),(coords(2,2)),(coords(2,3)),'filled');
hold off

end