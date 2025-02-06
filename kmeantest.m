
clear
s = rng;
if exist('parfor') == 0
disp('Please install Parallel computing toolbox')
return;
end
if exist('gurobi')==0
disp('Please install Gurobi')
return;
end

if exist('mexEMD')==0
disp('mexEMD missing. Please add the src directory and subfolders to the path')
return;
end

if exist('cl_run_kernel')==0
disp('cl_run_kernel missing. Please add the src directory and subfolders to the path and compile. You can set a global variable forceCPU = true to run on CPU, but this is slower');
%Uncomment the rows below to skip GPU implementation, 
%global forceCPU;
%forceCPU = true;
return;
end

markers = 'x.os';

%Running the different problems in the paper in Section 5.2

if true
    load ruspini.mat;
    t = tic;
    X = (ruspini-mean(ruspini))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    Xn = diag(1./normx)*U'*X;

    [y,g]=GlobalKmean(Xn,3,1,1,tic,normx);
    ruspuniobj = sum(diag(X'*X))-sum(diag((g'*X')*(X*g))./sum(g)')
    ruspunitime = toc(t)
end


if true
    load iris.mat;
    t = tic;
    X = (A-mean(A))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    normx2 = 1;
    Xn = diag(normx2./normx)*U'*X;
    [y,g]=GlobalKmean(Xn,3,1,1,tic,(normx./normx2));
    irisobj = sum(diag(X'*X))-sum(diag(g'*X'*X*g)./sum(g)')
    iristime = toc(t)
end

if true
    load rice.mat;
    t = tic;
    X = (rice-mean(rice))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    normx2 = 1;
    Xn = diag(normx2./normx)*U'*X;
    [y,g]=GlobalKmean(Xn,2,1,1,tic,(normx./normx2));
    riceobj = sum(diag(X'*X))-sum(diag(g'*X'*X*g)./sum(g)')
    ricetime = toc(t)
end

if true
    n = 666;
    k = 3;
    load gr666.mat;
    t = tic;
    X = (gr666v-mean(gr666v))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    normx2 = 1;
    Xn = diag((normx2./normx))*U'*X;
    [y1,g1]=GlobalKmean(Xn,k,1,1,tic,(normx./normx2));
    gr666obj = sum(diag(X'*X))-sum(diag(g1'*X'*X*g1)./sum(g1)')
    gr666time = toc(t)
end

if true
    load ecoli.mat;
    t = tic;
    X = (ecolidata-mean(ecolidata))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    normx2 = 1;
    Xn = diag(normx2./normx)*U'*X;
    [y,g]=GlobalKmean(Xn,2,1,1,tic,normx./normx2);
    ecoliobj = sum(diag(X'*X))-sum(diag(g'*X'*X*g)./sum(g)')
    ecolitime = toc(t)
end


if true
    load wine.mat;
    t = tic;
    X = (wines-mean(wines))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));


    normx = normx/normx(1);

    normx2 = 1;

    Xn = diag(normx2./normx)*U'*X;

    [y,g]=GlobalKmean(Xn,2,1,1,tic,(normx./normx2));
    wineobj = sum(diag(X'*X))-sum(diag(g'*X'*X*g)./sum(g)')
    winetime = toc(t)
end

if true
    load realestate.mat;
    t = tic;
    X = (realestate-mean(realestate))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    normx2 = 1;
    Xn = diag(normx2./normx)*U'*X;
    [y,g]=GlobalKmean(Xn,2,1,1,tic,(normx./normx2));
    realestateobj = sum(diag(X'*X))-sum(diag(g'*X'*X*g)./sum(g)')
    realestatetime = toc(t)
end


if true
    load ai4i.mat;
    ai4i = ai4i(2:end,:);
    t = tic;
    X = (ai4i-mean(ai4i))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    normx2 = 1;
    Xn = diag(normx2./normx)*U'*X;
    [y,g]=GlobalKmean(Xn,2,1,1,tic,(normx./normx2));
    ai4iobj = sum(diag(X'*X))-sum(diag(g'*X'*X*g)./sum(g)')
    ai4itime = toc(t)
end


if true
    load liver.mat;
    t = tic;
    X = (liver-mean(liver))';
    [U,S,V] = svd(X);
    normx = diag(S(1:size(X,1),1:size(X,1)));
    normx = normx/normx(1);
    normx2 = 1;
    Xn = diag(normx2./normx)*U'*X;
    [y,g]=GlobalKmean(Xn,2,1,1,tic,(normx./normx2));
    liverobj = sum(diag(X'*X))-sum(diag(g'*X'*X*g)./sum(g)')
    livertime = toc(t)
end



% Synthetic tests in the Performance study Section 5.1
NSe = [50 500 5000];
m =2; %Dimensions. Change to what you want
k =3; % Clusters. Change to what you want
sigma =1; %Standard deviation of normal distribution.

for nk = 1:length(NSe)
    for iter = 1:10
        n = NSe(nk);

        
        X=rand(m,n);


        global totalbranches
        totalbranches = 0;


        clear g;
        clear g2;



        ns = [ceil(n*0.35) floor(0.32*n) floor(0.33*n)];

        c0 = [ones(ceil(n*0.35),1);2*ones(floor(0.32*n),1);3*ones(floor(0.33*n),1)];



        X = [randn(m,ns(1))*sigma (randn(m,ns(2))*sigma+[0;2;zeros(m-2,1)]) (randn(m,ns(3))*sigma+[2;0;zeros(m-2,1)])];

        X = X-mean(X,2);

        [U,S,V] = svds(X,size(X,1));
        normx = diag(S(1:size(X,1),1:size(X,1)));
        normx = normx/normx(1);
        Xn = diag(1./normx)*U'*X;




        t = tic;
        [yn,gn]=GlobalKmean(Xn,k,1,1,t,normx);
        tid = toc(t)

        sortorder = (sum(X)+1)*gn;
        [~,sortI] = sort(sortorder);
        gn = gn(:,sortI);

        yn = (X*gn)./sum(gn);

        objn = sum(pdist2(X',(normx.*yn)',"squaredeuclidean").*gn,"all")


    end
end
%load('resultsn10000d2k3.mat')
disp([ num2str(mean(res.totaliterations)) '&' num2str(mean(res.qptreeconstraints)) '&' num2str(mean(res.qptreesize)) '&' num2str(mean(res.totalbranches)) '&$' num2str(mean(res.tid)) '\pm' num2str( std(res.tid)) '$']);


