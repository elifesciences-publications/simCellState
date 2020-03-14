function sysProp = sysProperties(typeIn, In1, In2, In3)
%  sysProp = sysProperties(typeIn, In1, In2, In3)
% This function compute the system properties for a given model. Input
% model can be specified as Cd-Ct matrices or T-L.
% The output structure contains:
% - strongly connected component
% - loops and loop connection
% - stability of the system
% - dgraph
% - classification of nodes: R (renewing), C (committed), D (developing), O (other), 0 (death)
% 
% Examples:
%  sysProp = sysProperties('TL', T, L, N)
%  sysProp = sysProperties('CdCt', Cd, Ct, [])

% Input type and conversions
switch typeIn
    case 'TL'
        T = In1; L = In2; N = In3;
        [A, J] = TL2AJ(T, L, N);
        [Cd, Ct] = TL2CdCt(T, L, N);
    case 'CdCt'
        Cd = In1; Ct = In2;
        [A, J] = CdCt2AJ(Cd, Ct);
        N = size(A, 1);
end

% strongly connected components
R = sparse(A'); % sparse inputs
[S, C, sSCC] = findSCC(R);
[eigA, eigJ] = eigSCC(A, J, S, C);
% graph
% G = digraph(R);
Ginv = digraph(R');
mut = diag(A-J)'; mut(mut==0) = 1;
Ap = A./repmat(mut, N, 1);
G = digraph(sparse(Ap')); % sparse((full(R')./repmat(mut, N, 1))')); % probability

% Nodes classification: renewing, committed and developing
indxR = find(C);
indx = setdiff(1:1:N, indxR);
indx0 = find(sum(A) == 0); % death
indx = setdiff(indx, indx0);
indxC_all = findOpenEnd(G, indx, indxR); % committed
indxD_all = findOpenEnd(Ginv, indx, indxR); % developing
indxC = setdiff(indxC_all, indxD_all);
indxD = setdiff(indxD_all, indxC_all);
indxO = setdiff(setdiff(indx, indxC), indxD); % all the other

% check stability
isct = sum(Ct~=0);
detFlag = all(isct(sum(Cd)~=0)==0) & all(isct <=1);
if detFlag % deterministic model
    stability = detModelStability(A./repmat(mut, N, 1), S, C, G);
else % stocastic  model
    stability = stoModelStability(eigJ);
end
stabFlags = stabilityFlags(eigJ); % at least one SSC is stable-growing-decaying

% new vertices and SCCs classification
% single nodes that are not death
indxSN = setdiff(find(C==0), indx0);
sccType = NaN(1, S+length(indxSN));
sccIndx = cell(1, S+length(indxSN));
for is = 1:length(indxSN)
    sccType(is) = sign(J(indxSN(is), indxSN(is)));
    sccIndx{is} = indxSN(is);
end
for is = 1:S % all 
    maxEig = max(real(eigJ{is})); maxEig(abs(maxEig)<1e-4) = 0;
    sccType(length(indxSN)+is) = sign(maxEig);
    sccIndx{length(indxSN)+is} = find(C == is);
end
vertexType = -ones(1, N); % -1 is other; 0 is death; 1 is developing; 2 is renewing and 3 is committed
vertexType(indx0) = 0; % death
indxCrit = find(sccType == 0);
for ii = 1:length(indxCrit)
    vertexType(sccIndx{indxCrit(ii)}) = 2; % renewing
end
indx1 = findOpenEnd(G, find(vertexType == -1), find(vertexType == 2));
indx2 = findOpenEnd(Ginv, find(vertexType == -1), find(vertexType == 2));
vertexType(setdiff(indx1, indx2)) = 3; % committed
vertexType(setdiff(indx2, indx1)) = 1; % developing

% Save system properties
sysProp = struct('N', N, 'Cd', Cd, 'Ct', Ct, 'A', A, 'Ap', Ap, 'J', J, 'digraph', G, ...
    'indxR', indxR, 'indxD', indxD, 'indxC', indxC, 'indx0', indx0, 'indxO', indxO,...
    'NSSC', S, 'SCC', C, 'sizeSCC', sSCC, ...
    'stability', stability, 'detFlag', detFlag, 'mut', mut, ...
    'sccType', sccType, 'vertexType', vertexType);
sysProp.sccIndx = sccIndx;
sysProp.eigA = eigA;
sysProp.eigJ = eigJ;
sysProp.stabFlags = stabFlags;

end

% buildin functions
% -----------------
function indxOut = findOpenEnd(G, indx, indxR)

indxOut = []; % initialization
for ii = 1:length(indx)
    T = dfsearch(G, indx(ii));
    if isempty(intersect(T, indxR))
        indxOut = [indxOut indx(ii)];
    end
end
indxOut = unique(indxOut);

end


function stability = detModelStability(M, S, C, G)

% find loops and connections
sumLoop = cell(1, S);
connectedLoop = cell(S,S);
for isi = 1:S
    idinode = find(C == isi);
    % loops
    scFi = M(idinode,idinode);
    sumLoop{isi} = sum(scFi);
    % connections
    T = dfsearch(G, idinode(1));
    for iso = setdiff(1:S, isi)
        idonode = find(C == iso);
        connectedLoop{iso, isi} = intersect(T,idonode);
    end
end

% stability
if S == 0 % no scc
    istab = false;
else
    istab = true;
end
for isi = 1:S % check of loops within SCC
    istab = istab & all(sumLoop{isi} == 1);
    for iso = 1:S % check disjoint loops
        istab = istab & isempty(connectedLoop{iso, isi});
    end
end

% output
stability = struct('flag', istab);
stability.sumLoop = sumLoop;
stability.connectedLoop = connectedLoop;

end


function stability = stoModelStability(eigJ) % TBD in case of more than one SCC

tol = 1e-14;

istab = true;
for is = 1:length(eigJ)
    if stabilitySCC(eigJ{is}) > tol
        istab = false;
    end
end

% output
stability = struct('flag', istab);

end

function stabFlags = stabilityFlags(eigJ)

tol = 1e-13;

eigJmaxReal = cellfun(@(x) max(real(x)), eigJ);
eigJmaxAbs = cellfun(@(x) min(abs(x)), eigJ);

% initialization
stabFlags = zeros(1,3); % stable-growing-decaying
if any(eigJmaxAbs < tol)
    stabFlags(1) = 1;
end
if any(eigJmaxReal > tol)
    stabFlags(2) = 1;
    if length(eigJ) == 1 % careful with this in case of multiple SCC
        % disp('Not stable (1 SCC with both 0 and positive eigenvalues)!')
        % stabFlags(1) = 0;
    end
end
if any(eigJmaxReal < -tol)
    stabFlags(3) = 1;
end

end