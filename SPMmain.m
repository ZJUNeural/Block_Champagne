addpath(genpath('./'));
%% load SPM data and fMRI maps
load('SPMdata');
nSource = size(L,2);
%% visualize fMRI maps

temp = zeros(size(pQ{1}));
for iter = 1:length(pQ)
    temp = temp+pQ{iter};
end
temp = temp(:);
% visualize fMRI activation
visual_source(temp,spm_cortex);
%% ESI with fMRI prior
sgn = -1;

% ESI of T1
load('PriorStage1');
segData = fam_data(:,1:47);
J_1 = BlockChampagne(L,segData,spm_cortex.VertConn,[],'atlas',atlas,'knb',2);
visual_source(sgn*J_1(1:nSource,40),spm_cortex,0.1);
visual_source(sum(J_1(1:nSource,:).^2,2),spm_cortex,0.1);
%%
% ESI of T2
load('PriorStage2');
segData = fam_data(:,48:end);
J_2 = BlockChampagne(L,segData,spm_cortex.VertConn,[],'atlas',atlas,'knb',2);
visual_source(sgn*J_2(1:nSource,6),spm_cortex,0.1);
visual_source(sum(J_2(1:nSource,:).^2,2),spm_cortex,0.1);
