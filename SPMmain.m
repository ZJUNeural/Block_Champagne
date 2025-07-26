addpath(genpath('./'));
%% load SPM data and fMRI maps
load('SPMdata');
load('base_data_famous');
nSource = size(L,2);
%% visualize fMRI maps

temp = zeros(size(pQ{1}));
for iter = 1:length(pQ)
    temp = temp+pQ{iter};
end
temp = temp(:);
% visualize fMRI activation
visual_source(temp,spm_cortex);
%% Whiten
CovN = cov(base_data_famous');
[U,S,~] = svd(CovN);
S = diag(S);
rankn = rank(CovN);
S = S(1:rankn); U = U(:,1:rankn);
whiten = diag(sqrt(1./S))*U';
L = whiten*L;
fam_data = whiten*fam_data;
%% ESI with fMRI prior
sgn = -1;

% ESI of T1
load('PriorStage1');
temp = struct2table(atlas);
Atlas = temp.Vertices;
segData = fam_data(:,1:47);
[J_1,par1] = BlockChampagne(L,segData,spm_cortex.VertConn,[],'atlas',Atlas,'knb',1,'maxiter',100);

visual_source(sum(J_1(1:nSource,:).^2,2),spm_cortex,0.01);
%%
[J_1n,par1n] = BlockChampagne(L,segData,spm_cortex.VertConn,[],'knb',2,'maxiter',100);
visual_source(sum(J_1n(1:nSource,:).^2,2),spm_cortex,0.1)

%%
% ESI of T2
load('PriorStage2');
 temp = struct2table(atlas);
Atlas = temp.Vertices;
segData = fam_data(:,48:end);
[J_2,par2] = BlockChampagne(L,segData,spm_cortex.VertConn,[],'atlas',Atlas,'knb',1,'maxiter',100);

visual_source(-sum(J_2(1:nSource,:).^2,2),spm_cortex,0.01,'pnactive');
%%
[J_2n,par2n] = BlockChampagne(L,segData,spm_cortex.VertConn,[],'knb',2,'maxiter',100);
visual_source(-sum(J_2n(1:nSource,:).^2,2),spm_cortex,0.01,'pnactive');