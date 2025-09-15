function [S,par] = BlockChampagne(L,B,MM,Cnoise,varargin)
% solve the sparse bayesian learning source imaging problem with block-champagne algorithm 
% function [S,par] = BlockChampagne(L,B,MM,Cnoise,varargin)
%   Input:
%       L(nSensor*(nSource*nd)):leadfield matrix
%       B(nSensor*nSnap): E/MEG measurement
%       MM:adjacent matrix for each dipole
%       Cnoise(nSensor*nSensor or []): noise covariance, can be learned from data, default:[] 
%   Optional:
%       maxiter:maximum iterations, default:50
%       epsilon:tolerance for stopping algorithm, default:1e-8
%       isextended(0 or 1):learn noise series with augmented leadfield matrix, default:1
%       updatecn(0 or 1):learn noise covariance from data, default:1
%       knb(int):construct blocks with preceding k order neighbors for each
%           dipole, default:1 (used if blockOption.type == 'Conn')
%       atlas:prior constraints
%% Initialize
Knb = 1;
[nSensor,nSource] = size(L); % number of channels & source dipoles
updateCn = 0;
updateR = 0; % use fixed r here
maxIter = 50;
epsilon = 1e-8;
nSamp = size(B,2); % snapshots
if norm(B,'fro')<sqrt(numel(B))
    B = B./norm(B,'fro')*sqrt(numel(B));
end
if nargin<4||isempty(Cnoise)
    Cnoise = 1e-4*eye(nSensor);
    updateCn = 1;
end

Covy = B*B'/nSamp;
isextended = 1;
isAtlas = 0;
Rcorr = 0.99;newRcorr = Rcorr;
if nargin>4
    for inar = 1:2:length(varargin)
        Param = lower(varargin{inar});
        Value = varargin{inar+1};
        switch Param
            case 'maxiter'
                maxIter = Value;
            case 'epsilon'
                epsilon = Value;
            case 'updater'
                updateR = Value;
            case 'knb'
                Knb = Value;
            case 'isextended'
                isextended = Value;
            case 'atlas'
                Atlas = Value;
                isAtlas = 1;                
        end
    end
end

if isextended
    LF = [L eye(nSensor)]; % leadfield: source & noise gain
    naSource = nSource+nSensor;
    M = blkdiag(MM^Knb,zeros(nSensor));
    M = M+eye(size(M));
else
    LF = L; % leadfield: source & noise gain
    naSource = nSource;
    M = MM^Knb;
    M = M+eye(size(M));
end
cls = sum(M>0,2);
G = sparse(naSource,sum(cls));
clusters = cell(naSource,1);
clusterBlk = clusters;
for it = 1:numel(cls)
    if it == 1
        clusters{it} = 1:sum(cls(1:it));
    else
        clusters{it} = sum(cls(1:it-1))+1 : sum(cls(1:it));
    end
    tidx = sort(unique([it,find(M(it,:)>0)]));
%     G(tidx,clusters{it}) = eye(numel(tidx));
    if cls(it) == 1
        clusterBlk{it} = 1;
        G(tidx,clusters{it}) = eye(numel(tidx));
    else
        tempBlk = Rcorr*ones(cls(it)) + eye(cls(it))*(1-Rcorr);
        [U,S] = eig(tempBlk);
        Blk.U = U;
        Blk.S = diag(S);
        clusterBlk{it} = Blk;
        G(tidx,clusters{it}) = U*diag(sqrt(Blk.S));
    end
end
% LF = LF*G;
LFG = LF*G;
Cgamma = (ones(naSource,1)+1e-3*randn(naSource,1))*trace(Covy)/trace(LF*LF.');
gamma = [];Ncls = zeros(naSource,1);
for it = 1:naSource
    Ncls(it) = numel(clusters{it});
    gamma = [gamma;Cgamma(it)*ones(Ncls(it),1)];
end

if isAtlas
    nAtlas = numel(Atlas);
    AtlasSet = cell(nAtlas,1);
    for iter = 1:nAtlas
        als_cls = []; tgamma = [];
        curVert = Atlas{iter};
        AtlasCluster = cell(numel(curVert),1);
        for idp = 1:numel(curVert)
            als_cls = union(als_cls,clusters{curVert(idp)});
            tgamma = [tgamma;Cgamma(curVert(idp))*ones(Ncls(curVert(idp)),1)];
            if idp == 1
                AtlasCluster{idp} = 1:numel(als_cls);
            else
                AtlasCluster{idp} = (AtlasCluster{idp-1}(end)+1):numel(als_cls);
            end
        end
        als.alscls = als_cls;
        als.resGam = zeros(numel(als_cls));
        als.gamma = diag(tgamma);
        als.AtlasCluster = AtlasCluster;
        AtlasSet{iter} = als;        
    end
end

dipIdx = 1:naSource;

cost = -1;
costlist = []; 
Gram = LFG.*repmat(gamma',nSensor,1)*LFG'+ Cnoise;
for iter = 1:maxIter    
    % gamma update
    
    % inverse of Gram
    [Ug,Sg] = svd(Gram);
    Sg = max(real(diag(Sg)),0);
    invSg = zeros(nSensor,1);
    ff = find(Sg>1e-6*max(Sg));
    invSg(ff) = 1./Sg(ff);
    invGram = Ug*diag(invSg)*Ug';

    W = repmat(gamma,1,nSensor).*LFG'*invGram;
    for it = 1:numel(clusters)
        if it == 1
            col = 1:sum(Ncls(1:it));
        else
            col = sum(Ncls(1:it-1))+1 : sum(Ncls(1:it));
        end
        ss = trace(W(col,:)*Covy*W(col,:)');
        z = trace(LFG(:,col)'*invGram*LFG(:,col));
        Cgamma(it) = (sqrt(z)./max(z,1e-8))*sqrt(ss);
    end  

    Cgamma(Cgamma>1e16) = 1e16;
    Cgamma(Cgamma<1e-12) = 1e-12;
    gamma = [];
    for it = 1:naSource
        gamma = [gamma;Cgamma(it)*ones(Ncls(it),1)];
    end     
    Gram = LFG.*repmat(gamma',nSensor,1)*LFG';
    if isAtlas % update the symmetric part
        for ials = 1:nAtlas
            als = AtlasSet{ials};
            als_cls = als.alscls;
            atlas_gamma = als.gamma;
        
            W_als = atlas_gamma*LFG(:,als_cls)'*invGram;
            ss = W_als*Covy*W_als';
            z = LFG(:,als_cls)'*invGram*LFG(:,als_cls);
%             [uz,sz] = svd(z);
%             szd = max(real(diag(sz)),0);
%             zz = uz*diag(sqrt(szd))*uz';
%             invszd = 1./sqrt(max(szd,1e-8));
%             invzz = uz*diag(invszd)*uz';            
%             [ps,ds] = svd(zz*ss*zz);
%             dsd = max(real(diag(ds)),0);
%             zsz = ps*diag(sqrt(dsd))*ps';   
%             cov_symmetric = invzz*zsz*invzz;
            temp = atlas_gamma - atlas_gamma*z*atlas_gamma;
            cov_symmetric = (ss+temp);
            Gram = Gram + LFG(:,als_cls)*cov_symmetric*LFG(:,als_cls)'-LFG(:,als_cls)*diag(gamma(als_cls))*LFG(:,als_cls)';

            als.gamma = cov_symmetric;
            W(als_cls,:) = cov_symmetric*LFG(:,als_cls)'*invGram;
            AtlasSet{ials} = als;   
        end
    end
    if updateR == 1
         wdiag = 0; wsubdiag = 0; Numdiag = 0; Numsubdiag = 0;
         for it = 1:nSource
             tidx = sort(unique([it,find(M(it,:)>0)]));
             col = clusters{it};
%              Blk = clusterBlk{it}; 
             Wit = G(tidx,clusters{it})*W(col,:);
             ss = (Wit*Covy*Wit');
             z = LFG(:,col)'*invGram*LFG(:,col);
%              temp = diag(gamma(col)) - Cgamma(it).^2.*z;
             temp = diag(gamma(col)) - diag(gamma(col))*z*diag(gamma(col));
             temp = G(tidx,clusters{it})*temp*G(tidx,clusters{it})';
             E_zz = (ss + temp)/nSamp;
             wdiag = wdiag+sum(diag(E_zz));
             Numdiag = Numdiag + numel(diag(E_zz));
             wsubdiag = wsubdiag+sum(diag(E_zz,1)); 
             Numsubdiag = Numsubdiag + numel(diag(E_zz,1));
         end
         newRcorr = (wsubdiag/Numsubdiag)/(wdiag/Numdiag);
        newRcorr = min(0.9999,newRcorr);
        newRcorr = max(newRcorr,0.6);
        kratio = newRcorr/Rcorr;
        Rcorr = newRcorr;
        for it = 1:nSource
            Blk = clusterBlk{it};          
            G(:,clusters{it}) = G(:,clusters{it})*diag(sqrt(Blk.S*kratio+1-kratio)./sqrt(Blk.S));
            clusterBlk{it}.S = Blk.S*kratio+1-kratio;
        end
        LFG = LF*G;         
    end

    %% Cnoise update
    if updateCn == 1        
        fw = eye(nSensor) - LFG*W;
        lam = sqrt(sum(fw*Covy.*fw,2)./diag(invGram));
        Cnoise = diag(max(real(lam),1e-8));
    end
    Gram = Gram + Cnoise;
    %% Cost
    cost_old = cost;
    cost = -0.5*sum(sum(Covy.*invGram))-0.5*(sum(log(max(Sg,1e-8)))+nSensor*log(2*pi));
    costlist = [costlist,cost];
    if abs((cost-cost_old)/cost)<epsilon
        break;
    end
end

gtW = G*W;
S = gtW*B;
par.W = gtW;
par.dipIdx = dipIdx;
par.gamma = Cgamma;
par.noise = Cnoise;
par.cost = costlist;
par.blkcorr = newRcorr;
end
