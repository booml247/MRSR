 
function [indx stump,error]=evalute(p,X,L,method,feature_num,lambda,lambda1,tau,DictSize)
stump = 0;
lab_val=unique(L);
mm=length(lab_val);
dat=[];
Lab=[];
for i=1:mm
      dat=[dat;X(L==lab_val(i),:)];
      Lab=[Lab;i*ones(sum(L==lab_val(i)),1)];
end
X=dat;
L=Lab;

switch method
    case 'RSR'
      [Z,~,obj] = RSR_l21(X,lambda) ;
      score= sum(Z.*Z,2);
      [~,indx] = sort(score,'descend');
    case 'MRSR'
      [W] = MRSR_l21(X,lambda,lambda1);  
      score= sum(W.*W,2);
      [~,indx] = sort(score,'descend');
    case 'FRFS'
      [indx stump] = FRFS(p,X,lambda); 
    case 'DPFS'
      [V,U]=DPFS(X',lambda,tau,DictSize);
       W=V;
       score= sum(W.*W,2);
       [~,indx] = sort(score,'descend');
    case 'DPFSP'
      [V,~,error]=DPFS_p(X',lambda,tau,DictSize,p);
       W=V;
       score= sum(W.*W,2);
       [~,indx] = sort(score,'descend');
    case 'L2FS'
      indx=L2FS(X,lambda);
      
    case 'EUFS'
      maxiter=20;
      X=normcol_equal(X);
      indx=EUFS_v1(X',DictSize,lambda,tau, maxiter);
    case 'EUFSW'
      maxiter=20;
      X=normcol_equal(X);
      [~,~,G]=EUFS_v1(X',DictSize,lambda,tau, maxiter);
      indx = FeatureSelection(X',G, size(X,2),tau);
    case 'AML21'
       W= lrra(X,X,lambda);
       score= sum(W.*W,2);
       [~,indx] = sort(score,'descend');
    case 'FRFSF'
      indx= RFS(X',X,lambda);
    case 'Variance' 
      indx = VAFS(X);
    case 'SparseFS'
      W=SparseFS(X,lambda);
      score= sum(W.*W,2);
      [~,indx] = sort(score,'descend');
    case  'MCFS' 
      options = [];
      options.k = 5; %For unsupervised feature selection, you should tune
      options.nUseEigenfunction = 5;  %You should tune this parameter.
      feature_num=150;
      indx = MCFS_p(X,feature_num,options);
    case  'Laplacian' 
     W = constructW(X);
     Y = LaplacianScore(X, W);
    [~,indx] = sort(-Y);
    case  'UDFS'
     para.k=5;
     para.lambda=lambda;
     LL = LocalDisAna(X', para);
     A = X'*LL*X;
     W=LquadR21_reg(A,mm,lambda);
     score= sum(W.*W,2);
     [~,indx] = sort(score,'descend');
    case  'SPEC'     
     W = X*X';
     w=fsSpectrum(W,X,-1); 
     [~,indx]=sort(w,'descend');  
    case  'MRSF'
     W = constructW(X);    
     [eigenvectors,Eigenvalues] = eig(full(W));
     Y = eigenvectors*sqrt(Eigenvalues);
     indx= mrsf( X, Y, 150);  
    case  'MRSFP'
%      options.Metric = 'Euclidean';
     options.NeighborMode = 'KNN';
     options.k = 5;
     options.WeightMode = 'HeatKernel';
     options.t = 1;
     options.PCARatio = 0.6;  % for LPP
     options.lambda =lambda;
     indx=MRSF_p(X,options);
    case 'fsfs'
     indx=fsfs(X,size(X,2),size(X,2)-feature_num);
    case 'MCFR'
     options = [];
     options.k = 5; %For unsupervised feature selection, you should tune
     options.nUseEigenfunction = 5;  %You should tune this parameter.
     options.lambda = lambda;
     indx = MCFR(X,feature_num,options);
    case 'MRFR'
     indx = MRFR(X,lambda);
     
    case 'RUFS'
     options.GraphDistanceFunction = 'cosine';
     options.Kernel = 'cosine';
     options.KernelParam = 0; % do not matter for consine kernel
     options.GraphWeights = 'distance';

    NN = 5;
    DISTFUNC = options.GraphDistanceFunction;
    KernelType = options.Kernel;
    KernelParam = options.KernelParam;
    llambda = 0.001;
    tic;
    L_LLR = calcLLR(X', NN, DISTFUNC, KernelType, KernelParam, llambda);
    t = toc;
    time_LLR = t;

    label = litekmeans(X,mm,'Replicates',10);
    for i=1:size(X,1)
        temp=zeros(1,mm);
        temp(1,label(i))=1;
        G0(i,:)=temp;
    end

    options.gamma = 1e-4;
    options.mu = 1e-1;
    options.beta = 1e-0;
    options.MaxIter = 5;
    options.epsilon = 1e-5;
    options.verbose = 1;

    options.nu = 1;
    options.alpha = 1;
    options.beta = lambda;
    tic;
    [ W, ~, ~ ] = RUFS(X, L_LLR, G0, options);
    t = toc;
    score = sum(W .* W, 2);
    [~, indx] = sort(score, 'descend');

   case 'DPFSM'
    options.GraphDistanceFunction = 'cosine';
    options.Kernel = 'cosine';
    options.KernelParam = 0; % do not matter for consine kernel
    options.GraphWeights = 'distance';

    NN = 5;
    DISTFUNC = options.GraphDistanceFunction;
    KernelType = options.Kernel;
    KernelParam = options.KernelParam;
    llambda = 0.001;
    tic;
    L = calcLLR(X', NN, DISTFUNC, KernelType, KernelParam, llambda);
    V=DPFS_p(X',L,lambda,tau,1,DictSize,p);
    W=V;
    score= sum(W.*W,2);
    [~,indx] = sort(score,'descend');
    case 'SDPFS'
     %Parameter setting
     DictSize = 30;
     tau    = 0.05;
     lambda = 0.003;
     gamma  = 0.0001;
     indx=SDPFS(X',L,DictSize, tau, lambda, gamma );
    case 'CUFS'
     lambda0=1;
     lambda3=1; 
     [W,obj]=CUFS(X,DictSize,lambda0,lambda,tau,lambda3);
     score= sum(W.*W,2);
     [~,indx] = sort(score,'descend');
     %Parameter setting
end

