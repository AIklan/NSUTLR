clear all
close all
clc
%[Y,XT,A,nc,nr,nb]=Syn(50);

%%
for SNR=50
       if SNR==30
        %           1    2     3   4    5    6    7    8    9    10   11   12
        %lambda = [1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1,5,10];;
        load DC2_30dB.mat
        lambda_sunsal=1e-2;
        
        lambda_clsunsal=1e-2;
        
        lambda_sunsal_tv=1e-2;lambdatv_sunsal_tv=1e-3;
        
        lambda_lcsu=1e-1;
        
        lambda_drsu=1;lambdatv_drsu=1e-3;
        
        lambda_wcsutv=1;beta_wcsutv=1e-3;
        
        lambda_sslrsu=1e-2;beta_sslrsu=1e-3;
        
        sw=3;lambda_glo=1e-1;lambda_loc=1e-3;
        load DC_SULRTen30dB
        load DC2_sunsal_L21_NLpatch_v2_30dB
        wnl_lambda1=lambda1;wnl_lambda2=lambda2;wnl_lambda3=lambda3;
        load DC2_NL_TSUn30dB
        load DC2_FastHySU_30dB_para_D24_H17_M30
    elseif SNR==40
        %           1    2     3   4    5    6    7    8    9    10   11   12
        %lambda = [1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1,5,10];;
        load DC2_40dB.mat
        lambda_sunsal=1e-3;
        
        lambda_clsunsal=1e-2;
        
        lambda_sunsal_tv=1e-3;lambdatv_sunsal_tv=1e-4;
        
        lambda_lcsu=1e-2;
        
        lambda_drsu=1e-1;lambdatv_drsu=1e-4;
        
        lambda_wcsutv=1;beta_wcsutv=1e-4;
        
        lambda_sslrsu=1e-2;beta_sslrsu=1e-2;
        sw=3;lambda_glo=5e-2;lambda_loc=1e-4;
        load DC_SULRTen40dB
        load DC2_sunsal_L21_NLpatch_v2_40dB
        wnl_lambda1=lambda1;wnl_lambda2=lambda2;wnl_lambda3=lambda3;
        load DC2_NL_TSUn40dB
        load DC2_FastHySU_40dB_para_D24_H17_M35
    elseif SNR==50
        %           1    2     3   4    5    6    7    8    9    10   11   12
        %lambda = [1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1,5e-1,1,5,10];;
        load DC2_50dB.mat
        lambda_sunsal=1e-3;
        
        lambda_clsunsal=1e-4;
        
        lambda_sunsal_tv=1e-3;lambdatv_sunsal_tv=1e-4; 
        
        lambda_lcsu=1e-4;
        
        lambda_drsu=1e-3;lambdatv_drsu=1e-4;
        
        lambda_wcsutv=1;beta_wcsutv=1e-4;
        
        lambda_sslrsu=1e-2;beta_sslrsu=1e-4;
        load DC2_FastHySU_50dB_para_D24_H17_M41
        load DC_SULRTen50dB
        load DC2_sunsal_L21_NLpatch_v2_50dB
        wnl_lambda1=lambda1;wnl_lambda2=lambda2;wnl_lambda3=lambda3;
        load DC2_NL_TSUn50dB.mat
        sw=3;lambda_glo=1e-2;lambda_loc=1e-5;
    end
    % nonnegative obervation
    % X = reshape(Y',nr,nc,nb);  % 3D HSI row*col*L
    
    % Sw = 4; P = round(nr*nc/Sw^2); Ws = 0.5;
    % seg = slic_HSI(X, P, Ws);
    
    %nmftype='sunsaltv';
    %% varing para first round
    
    % Time=zeros(1,7);
    % tic
    % [X_ncls] =  sunsal(A,Y,'lambda',0,'ADDONE','no','POSITIVITY','yes', ...
    %     'TOL',1e-4, 'AL_iters',200,'verbose','yes');
    % Time(1)=toc;
    % SRE_ncls = 20*log10(norm(XT,'fro')/norm(X_ncls-XT,'fro'));
    % %
    % SPA_ncls=length(find(X_ncls>0.005))/((size(X_ncls,1)*size(X_ncls,2)));
    %
    % RMSE_ncls=sqrt(mean2((X_ncls-XT).^2));
    
    %%
    tic
    [X_sunsal] =  sunsal(A,Y,'lambda',0.9*lambda_sunsal,'ADDONE','no','POSITIVITY','yes', ...
        'TOL',1e-4, 'AL_iters',200,'verbose','yes');
    Time(1)=toc;
    SRE_sunsal = 20*log10(norm(XT,'fro')/norm(X_sunsal-XT,'fro'));
    
    SPA_sunsal=length(find(X_sunsal>0.005))/((size(X_sunsal,1)*size(X_sunsal,2)));
    
    RMSE_sunsal=sqrt(mean2((X_sunsal-XT).^2));
    
    %%
    
    tic
    [X_clsunsal] =  clsunsal(A,Y,'POSITIVITY','yes','VERBOSE','yes','ADDONE','no', ...
        'lambda', 0.9*lambda_clsunsal,'AL_ITERS',200, 'TOL', 1e-8);
    Time(2)=toc;
    SRE_clsunsal = 20*log10(norm(XT,'fro')/norm(X_clsunsal-XT,'fro'));
    
    SPA_clsunsal=length(find(X_clsunsal>0.005))/((size(X_clsunsal,1)*size(X_clsunsal,2)));
    %
    RMSE_clsunsal=sqrt(mean2((X_clsunsal-XT).^2));
    %%
    tic
    [X_sunsaltv,res,rmse_i] = sunsal_tv(A,Y,'MU',0.05,'POSITIVITY','yes','ADDONE','no', ...
        'LAMBDA_1',0.9*lambda_sunsal_tv,'LAMBDA_TV', 0.9*lambdatv_sunsal_tv, 'TV_TYPE','niso',...
        'IM_SIZE',[100,100],'AL_ITERS',200, 'TRUE_X', XT,  'VERBOSE','yes');
    Time(3)=toc;
    
    SRE_sunsaltv = 20*log10(norm(XT,'fro')/norm(X_sunsaltv-XT,'fro'));
    SPA_sunsaltv=length(find(X_sunsaltv>0.005))/((size(X_sunsaltv,1)*size(X_sunsaltv,2)));
    %
    RMSE_sunsaltv=sqrt(mean2((X_sunsaltv-XT).^2));
    
%     %%
%     tic
%     Sw = 5; P = round(nr*nc/Sw^2); Ws = 0.5;
%     seg = slic_HSI(X, P, Ws);
%     % for i=1:length(lambda)
%     
%     [X_lcsu] =  lcsuv2(A,Y,'POSITIVITY','yes','SEG',seg,'VERBOSE','yes','ADDONE','no', ...
%         'lambda', lambda_lcsu,'AL_ITERS',200, 'TOL', 1e-8);
%     Time(4)=toc;
%     SRE_lcsu = 20*log10(norm(XT,'fro')/norm(X_lcsu-XT,'fro'));
%     
%     SPA_lcsu=length(find(X_lcsu>0.005))/((size(X_lcsu,1)*size(X_lcsu,2)));
%     %
%     RMSE_lcsu=sqrt(mean2((X_lcsu-XT).^2));
    
    
    %%
    
    
    %%
    
    %
    %
    
    %%

    %%
    
    tic
    parameter.lambda_s = 0.9*lambda_sslrsu;% nuclear norm
    parameter.lambda_p = 0.9*beta_sslrsu;% reweighted sparsity
    parameter.epsilon = 1e-5;
    parameter.maxiter = 200;
    parameter.mu = 0.1;
    parameter.xt = XT;
    
    X_sslrsu = sslrsu(Y, A, parameter);
    Time(4)=toc;
    
    SRE_sslrsu = 20*log10(norm(XT,'fro')/norm(X_sslrsu-XT,'fro'));
    
    
    SPA_sslrsu=length(find(X_sslrsu>0.005))/((size(X_sslrsu,1)*size(X_sslrsu,2)));
    
    RMSE_sslrsu=sqrt(mean2((X_sslrsu-XT).^2));
    %%

    parameter.lambda1 = 0.9*lambda1;
    parameter.lambda2 = 0.9*lambda2;
    parameter.mu      = 0.1;
    parameter.verbose = 0;
    parameter.MaxIter = 100;

    par.patsize       =   25;   % Patch size
    par.patnum        =   10;   % similar patch numbers
    par.step          =   floor((par.patsize-1));
    par.delta         =   0.1;
    %% run NL-TSUn
    tic
    [Out.X_hat] = NL_TSUn(A, Y,parameter, par);
    Time(5)=toc;
    X_tsun=Out.X_hat;
    SRE_nlsun = 20*log10(norm(XT,'fro')/norm(X_tsun-XT,'fro'));
    SPA_nlsun=length(find(X_tsun>0.005))/((size(XT,1)*size(XT,2)));
    RMSE_nlsun=sqrt(mean2((X_tsun-XT).^2));

    tic
    Sw = sw; P = round(nr*nc/Sw^2); Ws = log10(sqrt(SNR/3));
    seg = slic3Dhsi(X, P, Ws);
    %seg = slic_HSI(X, P, Ws);
    %for i=1:length(lambda)
    %for j=1:length(beta)
    parameter.lambda_s = 0.5*lambda_glo;
    parameter.lambda_p = 0.5*lambda_loc;% local sparse prior
    % optional parameters
    parameter.epsilon = 1e-5;
    parameter.maxiter = 200;
    parameter.mu = 0.20;
    parameter.xt = XT;
    parameter.verbose = 0;
    parameter.seg=seg;
    parameter.imgsize=[nr,nc];
    [X_lgsu,residual] = lgsuv2(Y, A, parameter);
    %X_lgsu = lgsuV3gpu(Y, A, parameter);
    Time(6)=toc;
    SRE_lgsu = 20*log10(norm(XT,'fro')/norm(X_lgsu-XT,'fro'));
    SPA_lgsu=length(find(X_lgsu>0.005))/((size(X_lgsu,1)*size(X_lgsu,2)));
    %
    RMSE_lgsu=sqrt(mean2((X_lgsu-XT).^2));

 tic
    [X_fast,Xt,Ap,D,iterf,idx]=runFastHySU(Y,A,9,0.005*lambdaHySU,0.01*betaHySU,10,[nr,nc],XT);
    Time(7)=toc;
    SRE_fast = 20*log10(norm(XT(1:q,:),'fro')/norm(X_fast-XT(1:q,:),'fro'));
    SPA_fast=length(find(X_fast>0.005))/((size(X_fast,1)*size(X_fast,2)));
    %
    RMSE_fast=sqrt(mean2((X_fast-XT(1:q,:)).^2));

    tic
    parameter.lambda_s = 0.8*lambda_s;
    parameter.lambda_p = 0.8*lambda_t;
    parameter.epsilon = 1e-5;
    parameter.maxiter = 200;
    parameter.mu = 0.10;
    parameter.xt = XT;
    parameter.imgsize=[nr,nc];
    parameter.verbose=0;
    %           [X_sunlrten,U] = SUDualSpasity(Y, A, parameter);
    [X_sulora,obj,res,resd,resp]=SULRTen(Y, A, parameter);
    %X_lgsu = lgsuV3gpu(Y, A, parameter);
    Time(8)=toc;
    SRE_sulora = 20*log10(norm(XT,'fro')/norm(X_sulora-XT,'fro'));
    SPA_sulora=length(find(X_sulora>0.005))/((size(X_sulora,1)*size(X_sulora,2)));
    %
    RMSE_sulora=sqrt(mean2((X_sulora-XT).^2));



    loadpara=load(strcat('DC_NonConvSULRTen',num2str(SNR),'dB.mat'));
    parameter.lambda_s = loadpara.lambda_s;
    parameter.lambda_p = loadpara.lambda_t;
    %parameter.logtol = tol(k);
    parameter.epsilon = 1e-5;
    parameter.maxiter = 200;
    parameter.mu = 0.50;
    parameter.xt = XT;
    parameter.imgsize=[nr,nc];
    parameter.verbose=0;
    tic
    [X_nonconvSULR,obj1,res,res_p1,res_d1] = NonConvexSULRTen(Y, A, parameter);
    Time(9)=toc;
    %X_sunlrten=SULRTen(Y, A, parameter);
    SRE_nonconvSULR = 20*log10(norm(XT,'fro')/norm(X_nonconvSULR-XT,'fro'));
    SPA_nonconvSULR=length(find(X_nonconvSULR>0.005))/((size(X_nonconvSULR,1)*size(X_nonconvSULR,2)));
    RMSE_nonconvSULR=sqrt(mean2((X_nonconvSULR-XT).^2));
    %end
    

    result=[SRE_sunsal,SPA_sunsal,RMSE_sunsal,Time(1);...
        SRE_clsunsal,SPA_clsunsal,RMSE_clsunsal,Time(2);...
        SRE_sunsaltv,SPA_sunsaltv,RMSE_sunsaltv,Time(3);...
        SRE_sslrsu,SPA_sslrsu,RMSE_sslrsu,Time(4);...
        SRE_nlsun,SPA_nlsun,RMSE_nlsun,Time(5);...
         SRE_lgsu,SPA_lgsu,RMSE_lgsu,Time(6);...
          SRE_fast,SPA_fast,RMSE_fast,Time(7);...
        SRE_sulora,SPA_sulora,RMSE_sulora,Time(8);...
        SRE_nonconvSULR,SPA_nonconvSULR,RMSE_nonconvSULR,Time(9)]

        filename= strcat('DC_Total_R1',num2str(SNR),'dB.mat');
    path='D:\Program Files\MATLAB\R2016a\a.phil\4.SparseUnmixing\MoneyInTheCode\nonConvexSULR\data\';
    save([path,filename],'result','X_sunsal','X_clsunsal','X_sunsaltv','X_sslrsu','X_tsun','X_lgsu','X_fast','X_sulora','X_nonconvSULR')

    %%

    %     figure
    %     p=9;
    %     for j=1:p
    %         subplot_tight(2, p, j,[.003 .003]);
    %         imagesc(reshape(XT(j,:)',nr, nc),[0,1]);axis image;axis off;
    %         subplot_tight(2, p, j+p,[.003 .003]);
    %         imagesc(reshape(X_sunlrten(j,:)',nr, nc),[0,1]);axis image;axis off;
    %     end
    %     drawnow;



end


 figure
    p=9;q=10;
    for j=1:p
        %===============X===========================
        subplot_tight(q, p, j,[.003 .003]);
        imagesc(reshape(XT(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============NCLS===========================
        subplot_tight(q, p, j+p,[.003 .003]);
        imagesc(reshape(X_sunsal(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============SUNSAL===========================
        subplot_tight(q, p, j+2*p,[.003 .003]);
        imagesc(reshape(X_clsunsal(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============CLSUNSAL===========================
        subplot_tight(q, p, j+3*p,[.003 .003]);
        imagesc(reshape(X_sunsaltv(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============SUNSALTV===========================
        subplot_tight(q, p, j+4*p,[.003 .003]);
        imagesc(reshape(X_sslrsu(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============LCSU===========================
        subplot_tight(q, p, j+5*p,[.003 .003]);
        imagesc(reshape(X_tsun(j,:)',nr, nc),[0,1]);axis image;axis off;
        subplot_tight(q, p, j+6*p,[.003 .003]);
        imagesc(reshape(X_lgsu(j,:)',nr, nc),[0,1]);axis image;axis off;
        subplot_tight(q, p, j+7*p,[.003 .003]);
        imagesc(reshape(X_fast(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============DRSU===========================
        subplot_tight(q, p, j+8*p,[.003 .003]);
        imagesc(reshape(X_sulora(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============DRSU===========================
        subplot_tight(q, p, j+9*p,[.003 .003]);
        imagesc(reshape(X_nonconvSULR(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============DRSU===========================
%         subplot_tight(q, p, j+8*p,[.003 .003]);
%         imagesc(reshape(X_wnltdsu(j,:)',nr, nc),[0,1]);axis image;axis off;
%         %===============DRSU===========================
%         subplot_tight(q, p, j+9*p,[.003 .003]);
%         imagesc(reshape(X_sulora(j,:)',nr, nc),[0,1]);axis image;axis off;
        %===============LGSU===========================
    end
    drawnow;
    colormap(jet)
    filenameRoc= strcat('DC2AbuMapsR1.png');%filenameRocFig= strcat(dctype,'Roc.fig');
    path='D:\Program Files\MATLAB\R2016a\a.phil\4.SparseUnmixing\MoneyInTheCode\nonConvexSULR\data\';%D:\Desktop\Papers\Others\SULRTen\SULoRa\data
    %saveas(ii,[path,filenameRoc]);
     print(figure(1),'-dpng', '-r300',[path,filenameRoc]);



% figure
% color1=[221,106,79]/255;
% color2=[0.286274509803922 0.580392156862745 0.768627450980392];%[134,192,125]/255;

% sre1=[4.11903589159810,5.07205095152670,10.5224355571433,15.0692618826506,10.6375829362698,4.26305974460700];
% sre2=[11.6949689326523,11.8067657412137,12.5215893792513,15.0692618826506,12.5995377824702,4.40613203897912]
% bound=20;
% sre1=[10.7559304195207,22.6275763684327,26.2989903099270,19.7281462936720,9.58057358130576,3.69705050047528];
% sre2=[24.4941166066957,24.9193427998280,26.2989903099270,23.1203375297322,11.4481000942130,3.90204325975059];
% bound=30;
% sre1=[18.1982964853735,36.6336945013535,33.0091639223799,19.1541400592877,9.15537049571538,3.66244173380081];
% sre2=[36.4300040417168,36.6336945013535,33.4595092238125,21.5947090645814,8.26362166284736,3.31321195958593];
% bound=40;
%------------------------------------------
% sre1=[6.68131673270853,7.52602359013855,11.8413125953554,18.0412689568474,15.7864340116609,6.34620304111818];
% sre2=[13.5464523720464,13.6457012437532,14.9685430938023,18.0412689568474,13.1326973127995,4.90727678697823];
% bound=20;
% sre1=[11.2819110515481,15.5894574321928,26.9860080824802,26.0952506511240,15.3992710658351,6.10728360373957];
% sre2=[25.8618169263614,26.2117427256146,26.9860080824802,22.6027941942319,10.6519214972590,3.62051811001536];
% bound=30;



%  sre1=[12.6349025668365,18.3472045006630,35.3047595935527,27.1447219148830,14.5611135917773,6.00545083653817];
%  sre2=[35.0357766589080,35.3047595935527,33.9481877802454,24.0699121562315,10.8286959951689,3.62973829682983];
% bound=40;
% yyaxis left;
% h=plot(sre1,'-o','LineWidth',2,'color',color1) ;set(h,'MarkerFaceColor',get(h,'color'));
% h=ylabel('SRE','FontSize', 14);set(h, 'FontName', 'Calibri', 'FontWeight', 'bold')
% set(gca,'YColor',color1);
% ylim([0,bound]);
% yyaxis right;
% h=plot(sre2,'-s','LineWidth',2,'color',color2) ;set(h,'MarkerFaceColor',get(h,'color'));set(gca,'FontSize',12);
% h=ylabel('SRE','FontSize', 14);
% set(h, 'FontName', 'Calibri', 'FontWeight', 'bold')
% set(gca,'YColor',color2);
% ylim([0,bound]);
% xlim([1, 6]);
% set(gcf,'unit','normalized','position',[0.4,0.4,0.20,0.25]);
% set(gca,'xtick',[1 2 3 4 5 6],'xticklabel',{'1e-5','1e-4','1e-3','1e-2','1e-1','1'},'looseInset',[0 0 0 0],'FontSize',10)
% xlabel('Parameters','FontName', 'Calibri', 'FontWeight', 'bold','FontSize', 14);
% h=legend('\alpha','\beta','location','northwest');
% set(h, 'FontSize', 14)
% grid on
