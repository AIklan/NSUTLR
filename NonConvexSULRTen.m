function [X,obj1,res,res_p1,res_d1] = NonConvexSULRTen(Y, A, para)

%%
%--------------------------------------------------------------
% load and test required parameters
%--------------------------------------------------------------
[L, K] = size(Y);
m = size(A, 2);
if size(A, 1) ~= L;error(['The sizes of hyperspectral data matrix Y and spectral ','library matrix A are inconsistent!']);end
if isfield(para, 'maxiter');MaxIter = para.maxiter;else MaxIter = 500;end
if isfield(para, 'epsilon');epsilon = para.epsilon;else epsilon = 1e-5;end
if isfield(para, 'epsilon');XT = para.xt;else error('The GT is not available'); end
if isfield(para, 'mu'); mu = para.mu;else mu = 0.1;end
if isfield(para, 'lambda_p');lambda1 = para.lambda_p;else error('The parameter lambda_p is missing!');end
%if isfield(parameter, 'seg');seg = parameter.seg;else error('The parameter seg is missing!');end
if isfield(para, 'imgsize');imgsize = para.imgsize;else error('The parameter lambda_s is missing!');end
if isfield(para, 'lambda_s');lambda2 = para.lambda_s;else error('The parameter lambda_s is missing!');end
if isfield(para, 'verbose');verbose = para.verbose;else verbose = 1;end
if verbose==1;figure; end



%%
%---------------------------------------------
%  Initializations
%---------------------------------------------


IF = (A'*A + 3*eye(m))^-1;
U = IF*A'*Y;


[U] =  sunsal(A,Y,'lambda',0,'ADDONE','no','POSITIVITY','yes', ...
    'TOL',1e-4, 'AL_iters',200,'verbose','yes');

V1 = A*U;
V2 = U;
V3 = U;
V4 = U;

D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
D4 = V4*0;

%current iteration number
i = 1;

%primal residual 
res_p = inf;

%dual residual
res_d = inf;

%error tolerance
tol = sqrt((3*m + L)/2*K/2)*epsilon;
fun='logarithm';
hfun_sg = str2func( [ fun '_sg'] ) ;
%%
%---------------------------------------------
%  ADMM iterations
%---------------------------------------------
j=1;
while j<=MaxIter
    Au=U - D2;
    W1=repmat(1./(sum(Au.*Au,2)+eps),1,size(Au,2));
    %W1 = repmat(1./(sqrt(sum(abs(U-D3),2))+eps),1, size((U-D3),2));
    W2=1./(sqrt(abs(Au))+eps);
    while (i <= 10) && ((abs(res_p) > tol) || (abs(res_d) > tol))

%while (i <= MaxIter) && ((abs(res_p) > tol) || (abs(res_d) > tol))
    if mod(i, 10) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
        V40 = V4;
    end
    %update U and V
    U = IF*(A'*(V1 + D1) + (V2 + D2) + (V3 + D3) + (V4 + D4));
    V1 = 1/(1+mu)*(Y + mu*(A*U - D1));
    Au=reshape((U - D2)',imgsize(1),imgsize(2),size(A,2));
    V2=Log_prox_tnn( Au,lambda1/mu );
    V2 = reshape(V2,imgsize(1)*imgsize(2),size(A,2))';
    V3 = solve_Lp(U - D3, lambda2/mu.*W1,1/2,1);
    V4 = max(U - D4, 0);
    %update D
    D1 = D1 - A*U + V1;
    D2 = D2 - U + V2;
    D3 = D3 - U + V3;
    D4 = D4 - U + V4;
    
    if mod(i, 10) == 1
        %object function
        obj = 1/2*norm(A*U - Y, 'fro');
        %primal residual
        res_p = norm([V1; V2; V3; V4] - [A*U; U; U; U], 'fro');
        %dual residual
        res_d = norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');
  
        if res_p > 10*res_d
            mu = mu*2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
            D4 = D4/2;
        elseif res_d > 10*res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
            D4 = D4*2;
        end
    end
    sre(i)=20*log10(norm(XT,'fro')/norm(U-XT,'fro'));
    xreds(i)=norm(U-XT,'fro');
    obj1(i) = 1/2*norm(A*U - Y, 'fro')+lambda2*sum(sum(abs(U)));
    res(i)=norm(A*U - Y, 'fro');
    %primal residual
    res_p1(i) = norm([V1; V2; V3; V4] - [A*U; U; U; U], 'fro');
    %dual residual
     res_d1(i) = norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');


    if verbose==1
    subplot_tight(5, 3, 1,[.08 .08]); plot(obj1,'-ro'); xlim([0 MaxIter]); title('Objective Values','fontsize',8);
    subplot_tight(5, 3, 2,[.08 .08]); plot(res_p1,'-r.');hold on;plot(res_d1,'-g.'); xlim([0 MaxIter]); title('Residuals','fontsize',8);
    subplot_tight(5, 3, 3,[.08 .08]); plot(sre,'-ro'); xlim([0 MaxIter]); title('SRE Values','fontsize',8);
  
    subplot_tight(5, 3, 4,[.08 .08]); imagesc(XT); xlim([0 MaxIter]); title('XT','fontsize',8);axis off
    subplot_tight(5, 3, 5,[.08 .08]); imagesc(V3); xlim([0 MaxIter]); title('Global Pers.','fontsize',8);axis off
    subplot_tight(5, 3, 6,[.08 .08]); imagesc(V2); xlim([0 MaxIter]); title('Local Pers.','fontsize',8);axis off

    subplot_tight(5, 3, 7,[.08 .08]); imagesc(reshape(U(2,:),imgsize)); title('Abundance #2','fontsize',8);axis off
    subplot_tight(5, 3, 8,[.08 .08]); imagesc(reshape(XT(2,:),imgsize)); title('GT #2','fontsize',8);axis off
    subplot_tight(5, 3, 9,[.08 .08]); imagesc(abs(reshape(XT(2,:),imgsize)-reshape(U(2,:),imgsize))); title('Diff #2','fontsize',8);axis off


    subplot_tight(5, 3, 10,[.08 .08]); imagesc(reshape(U(7,:),imgsize)); title('Abundance #7','fontsize',8);axis off
    subplot_tight(5, 3, 11,[.08 .08]); imagesc(reshape(XT(7,:),imgsize)); title('GT #7','fontsize',8);axis off
    subplot_tight(5, 3, 12,[.08 .08]); imagesc(abs(reshape(XT(7,:),imgsize)-reshape(U(7,:),imgsize))); title('Diff #7','fontsize',8);axis off
   
    subplot_tight(5, 3, 13,[.08 .08]); imagesc(reshape(U(9,:),imgsize)); title('Abundance #9','fontsize',8);axis off
    subplot_tight(5, 3, 14,[.08 .08]); imagesc(reshape(XT(9,:),imgsize)); title('GT #9','fontsize',8);axis off
    subplot_tight(5, 3, 15,[.08 .08]); imagesc(abs(reshape(XT(9,:),imgsize)-reshape(U(9,:),imgsize))); title('Diff #9','fontsize',8);axis off
    drawnow;
    end
    fprintf('i = %d, obj = %.4f,res_p = %.4f, res_d = %.4f, mu = %.1f, SRE=%.2f, ||X-XT||=%.2f, [lambda=%.1e, beta=%.1e]\n',...
        i, obj1(i), res_p1(i), res_d1(i), mu,sre(i),xreds(i),lambda2,lambda1);

    i = i + 1;
end
if i == MaxIter + 1
    display('Maximum iteration reached!');
end


j=j+1;
i=1;
end

X = U;
end

function W = SSpasity(S,imgsize)
nr=imgsize(1);nc=imgsize(2);
[p,N]=size(S);
S=reshape(S',nr,nc,p);
W=zeros(p,1);
for i=1:p; 
    Saux=imfilter(S(:,:,i),fspecial('average',3));W(i,1)=1./norm(Saux,'fro');
end
W=repmat(W,1,nr*nc);
end

function [X, tnn, trank] = prox_tnn(Y, rho)

dim = ndims(Y);
[n1, n2, n3] = size(Y);
n12 = min(n1, n2);
Yf = fft(Y, [], dim);%傅里叶变换
Uf = zeros(n1, n12, n3);
Vf = zeros(n2, n12, n3);
Sf = zeros(n12,n12, n3);

Yf(isnan(Yf)) = 0;%初始化
Yf(isinf(Yf)) = 0;

trank = 0;
for i = 1 : n3
    [Uf(:,:,i), Sf(:,:,i), Vf(:,:,i)] = svd(Yf(:,:,i), 'econ');%矩阵的svd
    s = diag(Sf(:, :, i));
    s = max(s - rho, 0);%阈值处理
    Sf(:, :, i) = diag(s);
    temp = length(find(s>0));
    trank = max(temp, trank);
end
Uf = Uf(:, 1:trank, :);
Vf = Vf(:, 1:trank, :);
Sf = Sf(1:trank, 1:trank, :);

U = ifft(Uf, [], dim);
S = ifft(Sf, [], dim);
V = ifft(Vf, [], dim);

X = tprod( tprod(U,S), tran(V) );

tnn = sum( diag( Sf(:,:,1) ) );

end
function [X, tnn, trank] = prox_tnn_w(Y, rho)

% The proximal operator of the tensor nuclear norm of a 3-order tensor
%
% min_X rho*||X||_*+0.5*||X-Y||_F^2
%
% Y     -    n1*n2*n3 tensor
%
% X     -    n1*n2*n3 tensor
% tnn   -    tensor nuclear norm of X
% trank -    tensor tubal rank of X

para = 5;    % sanDiego and others:5;    airport-1 and -2: 6;       airport-3and-4:10; 
dim = ndims(Y);
[n1, n2, n3] = size(Y);
n12 = min(n1, n2);
Yf = fft(Y, [], dim);
Uf = zeros(n1, n12, n3);
Vf = zeros(n2, n12, n3);
Sf = zeros(n12,n12, n3);

Yf(isnan(Yf)) = 0;
Yf(isinf(Yf)) = 0;

trank = 0;
endValue = int16(n3/2 + 1);
for i = 1 : endValue
    [Uf(:,:,i), Sf(:,:,i), Vf(:,:,i)] = svd(Yf(:,:,i), 'econ');
    s1 = diag(Sf(:, :, i));
    
    sikema = 10^-6;
    diagsj_w = 1./(s1+sikema);
    diagsj_w = diagsj_w/(diagsj_w(para));

    s=max(s1-rho*diagsj_w,0);
   
    
    Sf(:, :, i) = diag(s);
    temp = length(find(s>0));
    trank = max(temp, trank);
end
for j =n3:-1:endValue+1
    Uf(:,:,j) = conj(Uf(:,:,n3-j+2));
    Vf(:,:,j) = conj(Vf(:,:,n3-j+2));
    Sf(:,:,j) = Sf(:,:,n3-j+2);
end

Uf = Uf(:, 1:trank, :);
Vf = Vf(:, 1:trank, :);
Sf = Sf(1:trank, 1:trank, :);

U = ifft(Uf, [], dim);
S = ifft(Sf, [], dim);
V = ifft(Vf, [], dim);

X = tprod( tprod(U,S), tran(V) );
tnn =  diag( Sf(:,:,1) ) ;
end


function   x   =  solve_Lp( y, lambda, p, isWeighted )

% Modified by Dr. Weisheng Dong
if size(y,2)~=size(lambda,2)
    
    lambda=lambda(:,1);
end


if isWeighted==1
    J     =   2;
    tau   =  (2*lambda.*(1-p)).^(1/(2-p)) + p*lambda.*(2*(1-p)*lambda).^((p-1)/(2-p));
    x     =   zeros( size(y) );
    i0    =   find( abs(y) > tau );
    %lambda0=lambda(i0);
    if length(i0) >= 1
        if size(lambda,1)>1
            lambda  =   lambda(i0);
        end
        y0    =   y(i0);
        t     =   abs(y0);
        for  j  =  1 : J
            t    =  abs(y0) - p*lambda.*(t).^(p-1);
        end
        x(i0)   =  sign(y0).*t;

    end
    
else
    
    J     =   2;
    tau   =  (2*lambda.*(1-p))^(1/(2-p)) + p*lambda.*(2*(1-p)*lambda)^((p-1)/(2-p));
    x     =   zeros( size(y) );
    i0    =   find( abs(y) > tau );
    
    if length(i0) >= 1
        % lambda  =   lambda(i0);
        y0    =   y(i0);
        t     =   abs(y0);
        for  j  =  1 : J
            t    =  abs(y0) - p*lambda.*(t).^(p-1);
        end
        x(i0)   =  sign(y0).*t;
    end
    
end
end
