function [W] = MRSR_l21(X,lambda,lambda1)
% This routine solves the following nuclear-norm optimization problem,
% min |E|_2,1+lambda*tr(W'*X'*L*X*W)+lambda1*|W|_2,1
% s.t., X = X*W+E
% inputs:
%        X -- N*M data matrix, M is the data dimension, and N is the number
%             of data vectors.

maxIter = 20;
[NClus,NFea] = size(X);
%A=eye(NClus);
%for i=1:NClus
%    for j=1:NClus
%        A(i,j)=exp(-norm(X(i,:)-X(j,:))/0.5);
%    end
%end
%L=diag(sum(A))-A;
%L=Laplacian( X' );
L = Laplacian_LRGA(X');
epsilon=10^(-6);
deta=0.0001;
v1=1:NClus;
v2=1:NClus;
Q=0;

%% Initializing optimization variables
% intialize
GL = eye(NClus);
GLt=eye(NClus);
GR=eye(NFea);
GRt=eye(NFea);
W = eye(NFea);
%% Start main loop
iter = 0;
temp= X'* L * X ;

while iter<maxIter
%     tic
    iter = iter + 1;
    P=Q;
    %udpate W
    Wt=W;
    W = (X' *GL* X + lambda* temp + lambda1*GR) \ X' * GL * X;
    %update GL,GR
    
    for ii=1:NClus
       v1(ii)=1/(2 * max(norm(X(ii,:)-X(ii,:) * W,2),epsilon));
    end
    for jj=1:NFea
       v2(jj)=1/(2 * max(norm(W(jj,:),2),epsilon));
    end
    GL=diag(v1);
    GR=diag(v2);
%   Q=trace((X-X*W)'*GLt*(X-X*W))+lambda0*trace(W'*X'*L*X*W)^2/(2*trace(Wt'*X'*L*X*Wt))+lambda1*trace(W'*GRt*W);
    GLt=GL;
    GRt=GR;
%     if abs(Q-P)<deta
%         break
%     else
%         continue
%     end
    
        
end
%toc





