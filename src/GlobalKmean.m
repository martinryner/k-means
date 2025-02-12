function [Y0,LBGamma] = GlobalKmean(X,k,type,n0,timestart,normx)
delete backupqptree.mat
delete 'qptreebranch_*'
global totalbranches;
global LBthreshold;
LBthreshold = inf;
totalbranches = 0;
if type == 0
    [Y0,LBGamma] = GlobalKmeanfixed(X,k);
end

if type == 1



    Xnorm = 1;
    X = X/Xnorm;


    LBthreshold = inf;


    UB = inf;
    total_iterations= 0;
    qptreeconstraints=0;
    qptreesize = 0;

    localoptimums = [];



    %Take 6 at the time

    %totalactives = 1:(k*(size(X,1)+1)+1);
    totalactives = 1:((k-1)*(size(X,1)+1)+1);
    atthetime =ceil(length(totalactives));%floor(totalactives);%456451;%34551;%do all at once

    for sindex = 1:atthetime:length(totalactives)%1:k*size(X,1)+1 %There are n*k+1 extreme points in the starting simplex.
        delete backupqptree.mat
        delete 'qptreebranch_*'
        indices = sindex:(sindex+atthetime-1);
        indices = indices(indices<=length(totalactives));
        startactives= totalactives(indices)

        [Y02,LBGamma2,UB2,total_iterations2,qptreeconstraints2,qptreesize2] = GlobalKmeanNonFixednested(X,k,n0,[],[],timestart,localoptimums,[],Xnorm,startactives,normx);
        localoptimums = [localoptimums LBGamma2(:)];

        total_iterations= total_iterations+total_iterations2;
        qptreeconstraints=qptreeconstraints+qptreeconstraints2;
        qptreesize = qptreesize+qptreesize2;
        if UB2 < UB
            UB = UB2;
            Y0 = Y02;
            LBGamma = LBGamma2;
            if LBthreshold > UB
                LBthreshold = UB;
            end
        end

    end
    if exist(['resultsn' num2str(size(X,2)) 'd' num2str(size(X,1)) 'k' num2str(k) '.mat'],'FILE')==2
        load(['resultsn' num2str(size(X,2)) 'd' num2str(size(X,1)) 'k' num2str(k) '.mat'],'res')
    end

    if ~exist('res','var') || isempty(res)
        res.totaliterations=[];
        res.qptreeconstraints=[];
        res.qptreesize=[]
        res.tid = []
        res.totalbranches  =[];
    end

    res.totaliterations(end+1)=total_iterations;
    res.qptreeconstraints(end+1)=qptreeconstraints;
    res.qptreesize(end+1)=qptreesize;
    res.tid(end+1) = toc(timestart);
    res.totalbranches(end+1)  = totalbranches

    save(['resultsn' num2str(size(X,2)) 'd' num2str(size(X,1)) 'k' num2str(k) '.mat'],'res');



end
end

%A branch and bound wrapper around the algorithm handling the queue of
%branches.
function [Y0,LBGamma,UB,total_iterations,qptreeconstraints,qptreesize] = GlobalKmeanNonFixednested(X,k,n0,oldbranches,branchesb,timestart,localoptimums,oldbits,Xnorm,startactive,normx)

global LBthreshold;

Y0= [];
LBGamma=[];
UB=[];
total_iterations=[];
qptreeconstraints=[];
qptreesize=[];

[Y02,LBGamma2,UB2,total_iterations2,qptreeconstraints2,qptreesize2,newbranches,newbranchesb] = GlobalKmeanNonFixed(X,k,n0,oldbranches,branchesb,timestart,localoptimums,oldbits,Xnorm,startactive,normx);


UB = inf;
total_iterations= 0;
qptreeconstraints=0;
qptreesize = 0;
if ~isempty(LBGamma2)
    localoptimums = [localoptimums LBGamma2(:)];

    total_iterations= total_iterations+total_iterations2;
    qptreeconstraints=qptreeconstraints+qptreeconstraints2;
    qptreesize = qptreesize+qptreesize2;
    if UB2 < UB
        UB = UB2;
        Y0 = Y02;
        LBGamma = LBGamma2;
        if LBthreshold > UB
            LBthreshold = UB;
        end

    end

end



if ~isempty(newbranches)


    %Taking from the end of the list
    %n = size(newbranches,1)

    for n = size(newbranches,1):-1:1
        bit = [ones(n-1,1);-1];
        [Y02,LBGamma2,UB2,total_iterations2,qptreeconstraints2,qptreesize2] = GlobalKmeanNonFixednested(X,k,n0,[oldbranches;bit.*newbranches(1:n,:)],[branchesb bit'.*newbranchesb(1:n)],timestart,localoptimums,[oldbits;bit],Xnorm,startactive,normx);
        if ~isempty(LBGamma2)
            localoptimums = [localoptimums LBGamma2(:)];

            total_iterations= total_iterations+total_iterations2;
            qptreeconstraints=qptreeconstraints+qptreeconstraints2;
            qptreesize = qptreesize+qptreesize2;
            if UB2 < UB
                UB = UB2;
                Y0 = Y02;
                LBGamma = LBGamma2;
                if LBthreshold > UB
                    LBthreshold = UB;
                end

            end
        end
    end


end
end

%A version with fixed evenly large clusters
function [Y0,LBGamma,UB] = GlobalKmeanfixed(X,k)
c0 = sum(X.^2,"all");
max_iterations = 1000;

m = size(X,1);
n = size(X,2);
%Upper and lower bounds
e1 = zeros(1,k);
e1(1) = 1;
ubs = zeros(size(X,1),1);
lbs = zeros(size(X,1),1);




for kk = 1:m
    [~,G] = mexEMD(ones(n,1),(n/k)*ones(k,1),-X(kk,:)'*e1/sum(abs(X(kk,:)),"all") );
    tmp = X*G*k/n;
    ubs(kk) = tmp(kk,1);
    [~,G] = mexEMD(ones(n,1),(n/k)*ones(k,1),X(kk,:)'*e1/sum(abs(X(kk,:)),"all") );
    tmp = X*G*k/n;
    lbs(kk) = tmp(kk,1);
end

ub = ubs*ones(1,k);
lb = lbs*ones(1,k);

Ac = [];
bc = [];
LB = -inf;
UB = inf;
qptree = [];
LBGamma =[];



LBGammas = {};

qptree = bb4_kmean([],[],ub(:),lb(:),qptree);

disp('Adding Z1 = X1')
for kk1=1:m
    Zsum = zeros(m,k);
    Zsum(kk1,:) = 1;

    Acnew = Zsum(:)';
    nrm = norm(Acnew);
    Ac(end+1,:) = Acnew/nrm;
    bc(end+1) = sum(X(kk1,:))/nrm;

    qptree = bb4_kmean(Ac(end,:),bc(end)',ub(:),lb(:),qptree);


    Ac(end+1,:) = -Acnew/nrm;
    bc(end+1) = -sum(X(kk1,:))/nrm;
    qptree = bb4_kmean(Ac(end,:),bc(end)',ub(:),lb(:),qptree);

end





for kk = 1:max_iterations


    qptree = bb4_kmean(Ac(end,:),bc(end)',ub(:),lb(:),qptree);


    [UB,tin] = max(qptree.objs);
    x = qptree.EPs(:,tin);

    Y = reshape(x,size(X,1),k);
    Z = X'*Y;

    [~,Gammastar] = mexEMD(ones(n,1),(n/k)*ones(k,1),-Z/sum(abs(Z),"all") );
    y0 = X*Gammastar*k/n;




    %Make constraint

    y0 = y0(:);
    xn = x/norm(x);
    nrm = norm(xn);
    Ac(end+1,:) = xn/nrm;
    bc(end+1) = y0'*xn/nrm;
    if sum(y0.^2,"all")> LB
        LB = max(LB,sum(y0.^2,"all"))
        LBGamma = Gammastar;
    end

    tightconstr = false;
    if tightconstr

        [ANnew,Gamma,basis] =  TightconstraintsFixedsize(m,n,k,-Z(:)/sum(abs(Z),"all"),X,Ac,3);
        for tk = 1:length(ANnew)
            tmp = ANnew{tk};
            xn = tmp(:)/norm(tmp(:));
            if  min(sum(abs(Ac-xn'),2)) > 1000*sqrt(eps)

                Zt = X'*tmp;

                [~,Gammastart] = mexEMD(ones(n,1),(n/k)*ones(k,1),-Zt/sum(abs(Zt),"all") );


                y0 = X*Gammastart*k/n;

                %Make constraint

                y0 = y0(:);
                nrm = norm(xn);
                Ac(end+1,:) = xn/nrm;
                bc(end+1) = y0'*xn/nrm;


                qptree = bb4_kmean(Ac(end,:),bc(end)',ub(:),lb(:),qptree);


                if sum(y0.^2,"all")> LB
                    LB = max(LB,sum(y0.^2,"all"))
                    LBGamma = Gammastart;
                end

            else
                disp('already exist')
            end
        end

    end






    KMLB = c0-UB;
    KMUB = c0-LB;

    %gap = (UB-LB)

    gap = (KMUB-KMLB)/KMUB;
    vv = log10(length(qptree.objs));
    [KMLB KMUB gap vv]


    %gap = UB-LB
    if gap < 1e-6
        break;
    end


end

total_iterations = kk

Y0 = reshape(y0,size(X,1),k);
end





function [Y0,LBGamma,KMUB,total_iterations,qptreeconstraints,qptreesize,newbranches,newbranchesb] = GlobalKmeanNonFixed(X,k,n0,branches,branchesb,timestart,localoptimums,bits,Xnorm,startactive,normx)
global LBthreshold;
Y0=[];
LBGamma =[];
KMUB=[];
total_iterations=[];
qptreeconstraints=[];
qptreesize=[];
newbranches=[];
newbranchesb=[];

newbranches = [];
newbranchesb = [];





plotit = false;
extremepointmethod = true;
max_iterations = 4000;
constraintstolerance = 1e-9;
gaptolerance =1e-4;

Runwithassignments =true;


c0 = sum((diag(normx)*X).^2,"all");
m = size(X,1);

n = size(X,2);

mathcalX = [X;ones(1,n)];


%An ordering of the barycenters
orderb =0; %Set to inf if not using ordering constraints
orderA =[ones(1,size(X,1)) 1];
orderAX =orderA*mathcalX;


%Handle boundaries, the means must be inside the convex hull (which is
%approximated to a box) of X
Acn = [];
%e_i^T XG e_j <= max(X(i,:))*1^T G e_j
%kron(e_j^T,e_i^TX)vec(gamma) <= kron(e_j^T,max(X(i,:))*1^T) vec(gamma)
for k2 = 1:k
    for kk = 1:m
        e1 = zeros(1,k);
        e1(k2) = 1;
        Acn(end+1,:)= kron(e1,[1 -max(X(kk,:))]*[X(kk,:); ones(1,n)]);
    end
end
for k2 = 1:k
    for kk = 1:m
        e1 = zeros(1,k);
        e1(k2) = 1;
        Acn(end+1,:)= -kron(e1,[1 -min(X(kk,:))]*[X(kk,:); ones(1,n)]);
    end
end
bcn = zeros(2*m*k,1);

Z1b = sum(mathcalX,2);

Ac = [];
bc = [];
LB = -inf;
UB = inf;
qptree = [];
LBGamma =[];


if extremepointmethod

    if isfile('backupqptree.mat')


        load backupqptree.mat;

        lb = qptree.lb;
        zsimplex = qptree.zsimplex;
        Ac = qptree.A;
        bc = qptree.b;
        bc = reshape(bc,1,[]);%In this module, its a row vector

    else

        %Sloppy simplex
        if false
            %Bounding simplex
            for k2 = 1:k
                for kk = 1:size(X,1)+1
                    Der = zeros(size(X,1)+1,k);
                    Der(kk,k2) = 1;
                    Gl = LP2_Voronoi(mathcalX'*Der,n0,n,k,orderAX,orderb,Ac,bc,X);
                    %Gl = LP2_Voronoi(mathcalX'*Der,n0,n,k,[],orderb,Ac,bc,X);
                    tmpl = mathcalX(kk,:)*Gl;
                    lb(kk,k2) = tmpl(k2);
                end
            end

            Gu = LP2_Voronoi(-mathcalX'*ones(m+1,k),n0,n,k,orderAX,orderb,Ac,bc,X);

            zsimplex = mathcalX*Gu;
            zsimplex = zsimplex(:);

            qptree = bb4_kmean2_simplex([],[],lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,true,normx);

        else
            %Trimmed simplex

            %The flat subspace due to Z1 = X1
            v0 = kron(ones(1,k),eye(size(X,1)+1));
            v0 = v0/sqrt(k);
            %The non-flat subspace
            v1 = null(v0)';

            A = rand(size(v1,1));
            [ut,st,vt] = svd(A);

            V = [v0;v1];
            Vs.v0 = v0;
            Vs.v1 = v1;

            %Now hatz = V z, so the minimum in a direction d in V space is
            % min <d,hatz> = min < d,Vz> = min <d,V*mathcalX*Gamma>
            % = min<mathcalX'*V'*d,Gamma>
            %Remember, V is a tensor the V further down is the vectorized tensor


            %low bound (the first d+1 is fixed)
            v0d = zeros(m+1,k);
            v0d (end,:) = n/k;
            for nn = 1:size(V,1)
                [nn size(V,1)]
                zz = reshape(V(nn,:),m+1,k);


                Gut = LP2_Voronoi(mathcalX'*reshape(V(nn,:),m+1,k),n0,n,k,orderAX,orderb,Ac,bc,X);
                zsimplex2 = mathcalX*Gut;
                lb(nn,1) =V(nn,:)*zsimplex2(:);%The lower boound in the V space


                centr=X*Gut./sum(Gut);
                cost=pdist2((diag(normx)*X)',(diag(normx)*centr)','squaredeuclidean');
                lbcandidate = c0-sum(cost.*Gut,"all");


                if lbcandidate> LB

                    LB = max(LB,lbcandidate);

                    KMUB =Xnorm^2*( c0-LB);
                    LBGamma = Gut;
                    Y0 = (X*LBGamma)./sum(LBGamma);
                end





            end

            lb =V'*lb; %The constraint is so that  <V(nn,:),z-zsimplex2> >=0, i.e. V z = lb in this case, z = V^-1 lb = V^T*lb

            normal = sum(v1);
            Gu2 = LP2_Voronoi(-mathcalX'*reshape(normal,m+1,k),n0,n,k,orderAX,orderb,Ac,bc,X);

            centr=X*Gu2./sum(Gu2);
            cost=pdist2((diag(normx)*X)',(diag(normx)*centr)','squaredeuclidean');
            lbcandidate = c0-sum(cost.*Gu2,"all");


            if lbcandidate> LB

                LB = max(LB,lbcandidate);

                KMUB =Xnorm^2*( c0-LB);
                LBGamma = Gu2;
                Y0 = (X*LBGamma)./sum(LBGamma);
            end




            zsimplex = mathcalX*Gu2;
            zsimplex = zsimplex(:);



            qptree = bb4_kmean2_simplex([],[],lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx,true,Vs);

        end

        %max n points and sum_j (z_j) = sum_i x_i
        disp('Adding Z1 = X1')
        if false
            for kk1=1:m+1
                Zsum = zeros(m+1,k);
                Zsum(kk1,:) = 1;
                %Zsum = [zeros(m,k);ones(1,k)];
                Acnew = Zsum(:)';
                nrm = norm(Acnew);
                Ac(end+1,:) = Acnew/nrm;
                bc(end+1) = sum(mathcalX(kk1,:))/nrm;


                qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);


                Ac(end+1,:) = -Acnew/nrm;
                bc(end+1) = -sum(mathcalX(kk1,:))/nrm;

                qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);

            end
        else
            disp('SKIPPING Z1 = Z1. This is a test. Shouldnt affect at all since the subspace is already removed')
        end



        %Add the min and max boundaries for the z values bounded by the
        %convex hull of X, here approximated by a box

        disp('DZ >= n0')
        %positive DZ >= n0

        for kk1 = 1:k
            %  [kk1 max(qptree.objs)]
            Acnew = zeros(m+1,k);
            Acnew(end,kk1) = -1;
            bcnew = -n0;
            Ac(end+1,:) = Acnew(:)';
            bc(end+1) = bcnew;

            qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);
        end

        if true




            disp('very largest values constraint')
            for kk1 = 1:m
                for kk2 = 1:k
                    %   [kk1 kk2 max(qptree.objs)]
                    ej = zeros(1,m);
                    ej(kk1) = 1;
                    Acnew = [zeros(1,(kk2-1)*(m+1)) ej -max(X(kk1,:)) zeros(1,(k-kk2)*(m+1))];
                    Ac(end+1,:) = Acnew;
                    bc(end+1) = 0;

                    qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);


                    Acnew = -[zeros(1,(kk2-1)*(m+1)) ej -min(X(kk1,:)) zeros(1,(k-kk2)*(m+1))];
                    Ac(end+1,:) = Acnew;
                    bc(end+1) = 0;

                    qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);

                end
            end
        else
            disp('not adding largest values constraints')
        end





        if true && ~isinf(orderb)
            disp('Adding ordering constraint')
            %Add a^T z_j <= a^T z_{j+1}
            for kk2 = 1:k-1

                Acnew = [ zeros(1,(size(X,1)+1)*(kk2-1)) orderA -orderA zeros(1,(size(X,1)+1)*(k-1-kk2))   ];
                nrm = norm(Acnew);
                Ac(end+1,:) = Acnew/nrm;

                bc(end+1) = orderb;%0;%sqrt(size(X,2))/nrm;
                qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx,true);


            end
        else
            disp('NOT adding ordering constraint')
        end

        % File does not exist.
        save backupqptree.mat qptree;
    end



    totbits = [bits(1:end-1)];
    totbits = (totbits+1)/2;



    if isfile(['qptreebranch_' char(totbits'+48) '.mat'])
        filename = ['qptreebranch_' char(totbits'+48) '.mat'];
        load(filename,'qptree');
        disp(['Loading branch ' filename]);
        if ~isempty(branches)
            Acnew = branches(end,:);
            nrm = norm(Acnew);
            Acnew = Acnew/nrm;
            Ac(end+1,:) = Acnew;
            bc(end+1) = branchesb(end)/nrm;%+sqrt(eps);
            qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);
        end
        disp(['branch has' num2str(length(qptree.objs)) ' points']);

        Ac = qptree.A;
        bc = qptree.b;
        bc = reshape(bc,1,[]);%In this module, its a row vector
    else



        disp('Adding branching constraints')


        for kk2 = 1:size(branches,1)
            Acnew = branches(kk2,:);
            nrm = norm(Acnew);
            Acnew = Acnew/nrm;
            Ac(end+1,:) = Acnew;
            bc(end+1) = branchesb(kk2)/nrm;%+sqrt(eps);
            qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);
        end
    end






else

    %Add a^T z_j <= a^T z_{j+1}
    for kk2 = 1:k-1
        for t1 = 1:size(orderA,1)
            Acnew = [ zeros(1,(size(X,1)+1)*(kk2-1)) orderA(t1,:) -orderA(t1,:) zeros(1,(size(X,1)+1)*(k-1-kk2))   ];
            nrm = norm(Acnew)
            Ac(end+1,:) = Acnew/nrm;
            bc(end+1) = orderb/nrm;
        end
    end

    %Add the min and max boundaries for the z values bounded by the
    %convex hull of X, here approximated by a box
    for kk1 = 1:m
        for kk2 = 1:k
            ej = zeros(1,m);
            ej(kk1) = 1;
            Acnew = [zeros(1,(kk2-1)*(m+1)) ej -max(X(kk1,:)) zeros(1,(k-kk2)*(m+1))];
            nrm = norm(Acnew);
            Ac(end+1,:) = Acnew/nrm;
            bc(end+1) = 0;

            Acnew = -[zeros(1,(kk2-1)*(m+1)) ej -min(X(kk1,:)) zeros(1,(k-kk2)*(m+1))];
            nrm = norm(Acnew);
            Ac(end+1,:) = Acnew/nrm;
            bc(end+1) = 0;
        end
    end

    %max n points and sum_j (z_j) = sum_i x_i

    for kk1=1:m+1
        Zsum = zeros(m+1,k);
        Zsum(kk1,:) = 1;
        %Zsum = [zeros(m,k);ones(1,k)];
        Acnew = Zsum(:)';
        nrm = norm(Acnew)
        Ac(end+1,:) = Acnew/nrm;
        bc(end+1) = sum(mathcalX(kk1,:))/nrm;

        Ac(end+1,:) = -Acnew/nrm;
        bc(end+1) = -sum(mathcalX(kk1,:))/nrm;
    end

end

lastkkk = 2;
gap = 1;
KMUB = inf;
total_iterations = 0;
qptreeconstraints=size(qptree.A,1);
qptreesize=length(qptree.objs);




for kk = 1:max_iterations


    if extremepointmethod




        act= find(qptree.active);
        [UB,tin] = max(qptree.objs(act));
        x = qptree.EPs(:,act(tin));

        AA = [qptree.Alb;qptree.A];
        bb = [qptree.blb;qptree.b];
        base = [qptree.lbbases(:,act(tin));qptree.bases(:,act(tin))];
        xt = AA(base,:)\bb(base);

        if norm(xt-x)/norm(x) > sqrt(eps)
            disp(['interpolation and eq differs ' num2str(norm(xt-x))])
        end

        if isempty(x)
            return;
        end
        tmp = reshape(x,size(X,1)+1,k);

    else

        %As an alternative, the relaxed concave problem is solved by
        %Gurobi.
        gtic = tic;
        [xtest,qtest] = QPn([qptree.Alb;Ac],[qptree.blb;bc(:)],-100000000*ones(k*(m+1),1),100000000*ones(k*(m+1),1),n,m,k);
        if isempty(xtest)
            return;
        end
        tmp = xtest;
        x = tmp(:);
        UB = qtest;
        gtime =  toc(gtic)
    end


    Y = tmp(1:size(X,1),:);
    S = tmp(end,:);
    mathcalD = [2*(diag(normx.^2)*Y)./S;-diag(Y'*diag(normx.^2)*Y)'./S.^2];


    %Project mathcalD so that it is tangent to Z1 = mathcalX 1
    mathcalD = mathcalD-sum(mathcalD,2);



    Cb = mathcalX'*mathcalD;

    forcelocal = false;


    %least squares
    leastsquares  =false;
    LSGamma = [];
    if leastsquares && gap <0.1
        LSGamma = QPdirect(mathcalX,x,n0,n,m,k);



        if ~isempty(LSGamma)


            YY = diag(normx)*X*LSGamma;
            lbcandidate = sum(diag(YY'*YY)./(sum(LSGamma)'));
            if lbcandidate> LB
                forcelocal = true;
                LB = max(LB,lbcandidate);
                KMUB =Xnorm^2*( c0-LB);
                LBGamma = LSGamma;
                Y0 = (X*LBGamma)./sum(LBGamma);
            end

            fgamma = X*LSGamma;
            fgamma = sum(diag(fgamma'*fgamma)./(sum(LSGamma)'));

            Acnew = x-vec(mathcalX*LSGamma);

            bcnew = vec(mathcalX*LSGamma)'*(x-vec(mathcalX*LSGamma));
            bcnew = bcnew;

            nrm = norm(Acnew);
            Acnew = Acnew/nrm;
            bcnew = bcnew/nrm;
            mv = min(sum(abs([Ac bc']- [Acnew' bcnew]),2));
            if mv > constraintstolerance

                Ac(end+1,:) = Acnew(:)';
                bc(end+1) = bcnew;%Cb(:)'*Gammastar(:);
                qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,kk,Z1b,startactive,normx);
                forcelocal = true;

            else
                disp('constraint already exist')
            end

        end
    end




    centroid = (Y./S);


    if Runwithassignments
        mig =double( Assignment(X,centroid,normx));


        if  mathcalD(:)'*x-Cb(:)'*mig(:) < 0
            disp('oops')
        end

        if  sum(diff(sum(mathcalX*mig)')<-orderb)>0 || sum(sum(mig) == 0) >0 || mathcalD(:)'*x<Cb(:)'*mig(:)
            [Gammastartestg,basisg,~]  = LP2_Voronoi(-Cb,n0,n,k,orderAX,orderb,Ac,bc,X);
            shouldbezero = norm(Gammastartestg(:)-mig(:));

        else

            Gammastartestg = mig;
            basisg = find(mig(:)>0);
            basisg = [basisg;n*k+[1:k]'];%Why? Its recalculated further down before its needed.
            %find basisg

        end
    else
        [Gammastartestg,basisg,~]  = LP2_Voronoi(-Cb,n0,n,k,orderAX,orderb,Ac,bc,X);
    end




    YY = diag(normx)*X*Gammastartestg;
    lbcandidatet = sum(YY(1:end,:).^2./sum(Gammastartestg),"all");
    lbcandidateold = sum(diag(YY'*YY)./(sum(Gammastartestg)'));

    centr=X*Gammastartestg./sum(Gammastartestg);
    cost=pdist2((diag(normx)*X)',(diag(normx)*centr)','squaredeuclidean');
    lbcandidate = c0-sum(cost.*Gammastartestg,"all");
    if abs((lbcandidate-lbcandidateold)/abs(lbcandidate))> sqrt(eps)
        disp([ 'numerical instability in old UB calculation ' num2str(abs(lbcandidate-lbcandidateold))] );
    end

    if lbcandidate> LB
        forcelocal = true;
        LB = max(LB,lbcandidate);

        KMUB =Xnorm^2*( c0-LB);
        LBGamma = Gammastartestg;
        Y0 = (X*LBGamma)./sum(LBGamma);
    end



    Yt = (X*Gammastartestg)./sum(Gammastartestg);
    if abs(lbcandidate-LB)/LB < 0.01 && norm((Yt(:)-Y0(:)))/norm(Yt(:)) >0.2
        disp(['Close to optimal but far away ' num2str(Xnorm^2*( c0-lbcandidate))]);
    end

    %begin move
    Cb =  mathcalX'*mathcalD;
    Acnew = mathcalD(:)';
    Gammastartest = Gammastartestg;
    bcnew =  Cb(:)'*Gammastartestg(:);
    %Make constraint
    nrm = norm(Acnew);

    Acnew = Acnew/nrm;
    bcnew = bcnew/nrm;

    if Acnew*x < bcnew
        %Not good
        disp('wasnt able to cut off. Aborting');

        Acnew*x - bcnew
        qptreesize = length(qptree.objs)
        qptreeconstraints = size(Ac,1)
        Gammastar =  Gammastartest;
        break;

    else
        Gammastar =  Gammastartest;
    end





    mv = min(sum(abs([Ac bc']- [Acnew bcnew]),2));

    Ac(end+1,:) = Acnew(:)';
    bc(end+1) = bcnew;%Cb(:)'*Gammastar(:);

    if max(qptree.objs) < LB
        disp('all extremepoints are objective value infeasible')
    end

    qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,kk,Z1b,startactive,normx);
    if isempty(qptree.objs)
        disp('all extremepoints are objective value infeasible.')
        break;
    end
    if max(qptree.objs) < LB
        disp('all extremepoints are objective value infeasible')
    end



    if  extremepointmethod


        %Local optimum
        findlocaloptimum =false;%true;
        if    (forcelocal  ||( findlocaloptimum && gap < 0.01))


            %Local kmeans

            Gammalocal = Gammastartestg;
            zlocal = mathcalX*Gammastartestg;

            for nnn = 1:200


                Yl = zlocal(1:m,:);
                Sl = zlocal(end,:);


                localder = [2*diag(normx.^2)*Yl./Sl;-diag(Yl'*diag(normx.^2)*Yl)'./Sl.^2];
                Cbl =  mathcalX'*localder;

                if Runwithassignments
                    GammalocalnA = Assignment(X,(Yl./Sl),normx);
                    if    sum(diff(sum(mathcalX*GammalocalnA)')<-orderb)>0 || sum(sum(GammalocalnA) == 0) >0
                        [Gammalocaln,basisg,~]  = LP2_Voronoi(-Cbl,n0,n,k,orderAX,orderb,Acn,bcn,X);
                        shouldbezero = norm(Gammalocaln(:)-GammalocalnA(:));

                    else
                        Gammalocaln= GammalocalnA;
                    end
                else
                    [Gammalocaln,basisg,~]  = LP2_Voronoi(-Cbl,n0,n,k,orderAX,orderb,Acn,bcn,X);
                end

                zlocaln = mathcalX*Gammalocaln;
                if norm(Gammalocal(:)-Gammalocaln(:)) < sqrt(eps)
                    break;
                end
                zlocal = zlocaln;
                Gammalocal = Gammalocaln;
            end



            centr=X*Gammalocal./sum(Gammalocal);
            cost=pdist2((diag(normx)*X)',(diag(normx)*centr)','squaredeuclidean');
            lbcandidate =c0- sum(cost.*Gammalocal,"all");

            if lbcandidate> LB

                LB = max(LB,lbcandidate);
                KMUB =Xnorm^2*( c0-LB);
                LBGamma = Gammalocal;
                Y0 = (X*LBGamma)./sum(LBGamma);
            end

            YY =X*Gammalocal;
            Sn = sum(Gammalocal);
            mathcalDn = [2*diag(normx.^2)*YY./Sn;-diag(YY'*diag(normx.^2)*YY)'./Sn.^2];


            Acnew = mathcalDn(:);
            bcnew = mathcalDn(:)'*vec(mathcalX*Gammalocal);
            nrm = norm(Acnew);
            Acnew = Acnew/nrm;
            bcnew = bcnew/nrm;
            mv = min(sum(abs([Ac bc']- [Acnew' bcnew]),2));
            if mv > constraintstolerance

                Ac(end+1,:) = Acnew(:)';
                bc(end+1) = bcnew;%Cb(:)'*Gammastar(:);
                qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,kk,Z1b,startactive,normx);




            end


        end






    end




    if isempty(qptree.objs)
        disp('all extremepoints are objective value infeasible.')
        break;
    end




    %Pruning Z(end,:) must be integer
    pruneZend = false;
    if pruneZend
        maxvals = (max(qptree.EPs((1:k)*(m+1),:)'));
        minvals = (min(qptree.EPs((1:k)*(m+1),:)'));
        for pk = 1:k

            if maxvals(pk)> sqrt(eps)+floor(maxvals(pk))
                Acnew = zeros(m+1,k);
                Acnew(m+1,pk) = 1;
                bcnew = floor(maxvals(pk));
                Ac(end+1,:) = Acnew(:)';
                bc(end+1) = bcnew;
                qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,kk,Z1b,startactive,normx);
            end
        end
        for pk = 1:k
            if minvals(pk)+sqrt(eps)< ceil(minvals(pk))
                Acnew = zeros(m+1,k);
                Acnew(m+1,pk) = -1;
                bcnew = -ceil(minvals(pk));
                Ac(end+1,:) = Acnew(:)';
                bc(end+1) = bcnew;
                qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,kk,Z1b,startactive,normx);
            end
        end
    end



    KMLB =Xnorm^2*( c0-UB);
    KMUB =Xnorm^2*( c0-LB);





    gap = (KMUB-KMLB)/KMUB;


    if mod(kk,1) == 0
        vv = 0;
        if extremepointmethod
            vv = log10(length(qptree.objs));
        end
        disp(['LB: ' num2str(KMLB) ', UB:' num2str(KMUB) ', Gap: ' num2str(gap) ', EPs: ' num2str(length(qptree.objs)) ', active: ' num2str(sum(qptree.active)) ', elapsed time:' num2str(toc(timestart)) ' globalUB: ' num2str(LBthreshold)]);

    end
    if gap < gaptolerance
        break;
    end
    if toc(timestart) > 100000
        disp('Exceeding time limit')
        break;
    end

    if LBthreshold < KMLB
        disp('Lower bound higher than threshold')
        break;
    end

    gap2 = 1;
    if ~isinf(LBthreshold )
        gap2 = (LBthreshold-KMLB)/LBthreshold;
    end


    if sum(qptree.objs>LB) == 0
        disp('No extreme point left with value lower than the Upper bound')
        break;
    end


    sizethres = 5e5;%5e5;
    divisions = 0;

    while (length(qptree.objs) > sizethres) && gap > 0.2 && gap2 > 0.2 && divisions < 5 && max(qptree.objs) > LB

        divisions = divisions+1;
        disp(['Giving up. Dividing the branch. Size:' num2str(length(qptree.objs))]);




        origbranch = false;
        if origbranch

            EPs = qptree.EPs(:,qptree.objs> LB);
            objsEPs = qptree.objs(qptree.objs>LB);

            centroids = reshape(EPs,m+1,k,[]);
            centroids(1:end-1,:,:) = centroids(1:end-1,:,:)./centroids(end,:,:);
            centroids= centroids(1:end-1,:,:);
            for kk = 1:k
                kkcentroids = squeeze(centroids(:,kk,:));
                var(kk) = sum((kkcentroids-mean(kkcentroids,2)).^2,"all")/size(kkcentroids,2);
            end
            [~,kkmax] = max(var);


            kkcentroids = squeeze(centroids(:,kkmax,:));
            %centr = sum(kkcentroids.*objsEPs,2)/sum(objsEPs);

            centr = mean(kkcentroids,2);

            [~,maxi] = max(objsEPs);
            cmax = kkcentroids(:,maxi);

            %Testing to only include the maximum
            centr = (4*cmax + centr)/5;
            u = cmax-centr;
            %end test


            [u,s,v] = svds(kkcentroids-centr,1);


            newbranch = zeros(m+1,k);
            newbranch(:,kkmax) = [u; -u'*centr];



            newbranchb = 0;

        else
            [~, maxi] = max(qptree.objs);
            tmp = qptree.EPs(:,maxi);
            tmp = reshape(tmp,size(X,1)+1,k);

            Y = tmp(1:size(X,1),:);
            S = tmp(end,:);
            mathcalD = [2*(diag(normx.^2)*Y)./S;-diag(Y'*diag(normx.^2)*Y)'./S.^2];


            %Project mathcalD so that it is tangent to Z1 = mathcalX 1
            mathcalD = mathcalD-sum(mathcalD,2);
            newbranch = mathcalD(:);
            newbranchb=0;


            found = false;
            tryk = 1;

            [u,s,v] = svds(qptree.EPs-mean(qptree.EPs,2),size(qptree.EPs,1));


            while ~found

                if tryk == 1



                    newbranch = qptree.EPs(:,maxi);
                    newbranchb=0;


                else
                    newbranch = randn((size(X,1)+1)*k,1);
                    newbranch = u*s*newbranch;

                    newbranchb=0;



                end


                %Check sign of fcurrent optimum
                if newbranch(:)'*qptree.EPs(:,maxi) > newbranchb
                    newbranch = -newbranch;
                    newbranchb = -newbranchb;
                end



                nrm = norm(newbranch);
                newbranch = newbranch/nrm;
                newbranchb = newbranchb/nrm;

                testbranches(:,tryk) = newbranch(:);
                newbranch = newbranch(:)';

                minb = min(newbranch(:)'*qptree.EPs);
                maxb = max(newbranch(:)'*qptree.EPs);
                %How many points will the new branch get?
                vecc= [];
                newbranchtimesqptreeEPs =newbranch(:)'*qptree.EPs;

                for  btestk=1:50

                    btest = (50-btestk)/50*minb+(btestk)/50*maxb;

                    infeasiblesbit = (newbranchtimesqptreeEPs  > btest + sqrt(eps));
                    feasiblebit = (~infeasiblesbit) ;


                    feasibleadjacents= qptree.closematrix(:,infeasiblesbit)&feasiblebit';
                    netmoreextremepoints = sum(feasibleadjacents,"all")-sum( infeasiblesbit,"all" );



                    %  infeasiblesbit = (newbranchtimesqptreeEPs < btest +sqrt(eps));
                    % feasiblebit = (~infeasiblesbit) ;


                    %Flipping it instead
                    feasibleadjacents= qptree.closematrix(:,feasiblebit)&infeasiblesbit';
                    netmoreextremepoints2 = sum(feasibleadjacents,"all")-sum( feasiblebit,"all" );


                    vecc(btestk,4) = sum( feasiblebit,"all" );
                    vecc(btestk,2) = netmoreextremepoints;
                    vecc(btestk,3) = netmoreextremepoints2;
                    vecc(btestk,1) = btest;
                end


                idxpos = find(vecc(:,3) < 0 & vecc(:,2) < 0);
                if ~isempty(idxpos)
                    found = true;
                    %Multiobjective
                    [~,bestidx] = max(vecc(idxpos,3).*vecc(idxpos,2));
                    newbranchb = vecc(idxpos(bestidx),1);

                  %  figure(6951);
                  %  plot(vecc(:,2),vecc(:,3));
                  % drawnow;


                end
                %Backup plan.
                idxposback = find(vecc(:,2) < 0);
                [min3val(tryk), min3pos] = min(vecc(idxposback,2)-vecc(idxposback,3));
                min3pos = idxposback(min3pos);
                min3b(tryk) = vecc(min3pos,1);





                if ~found && tryk > 10
                    found = true;
                    [~,bestback] =min(min3val);
                    newbranchb = min3b(bestback);
                    newbranch =   testbranches(:,bestback);
                end

                tryk = tryk+1;

            end

            infeasiblesbit = (newbranch(:)'*qptree.EPs> newbranchb + sqrt(eps));
            feasiblebit = (~infeasiblesbit) ;


            feasibleadjacents= qptree.closematrix(:,infeasiblesbit)&feasiblebit';
            netmoreextremepoints = sum(feasibleadjacents,"all")-sum( infeasiblesbit,"all" );



            infeasiblesbit = (newbranch(:)'*qptree.EPs < newbranchb +sqrt(eps));
            feasiblebit = (~infeasiblesbit) ;


            feasibleadjacents= qptree.closematrix(:,infeasiblesbit)&feasiblebit';
            netmoreextremepoints2 = sum(feasibleadjacents,"all")-sum( infeasiblesbit,"all" );

            disp(['Found a cut which renders' num2str(netmoreextremepoints) ' and ' num2str(netmoreextremepoints2) ' net more points in its branches']);


        end






        newbranches(end+1,:) = newbranch(:)';
        newbranchesb(end+1) =newbranchb;

        %Add to this also
        branches(end+1,:) = newbranch(:)';
        branchesb(end+1) = newbranchb;

        Acnew = branches(end,:);




        nrm = norm(Acnew);
        Acnew = Acnew/nrm;
        Ac(end+1,:) = Acnew;
        bc(end+1) = newbranchb/nrm;


        %Save the old one
        totbits = [bits;ones(size(newbranches,1)-1,1)];
        totbits = (totbits+1)/2;

        %branchv = totbits.*(2.^(0:(length(totbits)-1))');
        filename = ['qptreebranch_' char(totbits'+48) '.mat'];
        save(filename, 'qptree');
        disp(['Saving branch ' filename]);


        global totalbranches;
        totalbranches = totalbranches+1;

        qptree = bb4_kmean2_simplex(Ac(end,:),bc(end)',lb(:),zsimplex,qptree,m,k,LB,0,Z1b,startactive,normx);




    end
end
if kk == max_iterations
    disp('warning: maximum iterations ')
end

qptreesize = length(qptree.objs)
qptreeconstraints = size(Ac,1)
total_iterations = kk
Y0 = Y0*Xnorm;






end








function [G,basis,result] = LP2_Voronoi(C,n0,n,k,orderAX,orderb,At,bt,X,S)
d = size(X,1);
%orderAX = [];%Skipping this
clear model;
%n = size(C,1);
%k = size(C,1);
At = [];
bt = [];
if ~isempty(At)
    At = At*kron(eye(k),[X;ones(1,n)]);
    bt = bt(:);
end


model.obj = C(:);

%G = ones(n,k);
%vec(eye(n)*G*ones(k,1)) = ones(n,1);
%kron(ones(k,1)',eye(n))* vec(G) = ones(n,1)
%vec(ones(1,n)*G*eye(k)) >= n0*ones(1,k);
%kron(eye(k),ones(1,n))* vec(G) >=n0*ones(k,1)
if  ~isempty(orderAX )
    Al = [];

    for kk = 1:k-1
        for t1 = 1:size(orderAX,1)
            e1 = zeros(k,1);
            e2 = e1;
            e1(kk) = 1;e2(kk+1) = 1;
            %Al = [Al;-kron(e1',ones(1,n))+kron(e2',ones(1,n))];
            Al = [Al;+kron(e1',orderAX(t1,:))-kron(e2',orderAX(t1,:))];
        end
    end


    orderbs = orderb*ones((k-1)*size(orderAX,1),1);




    if exist('S','var')

        model.A = sparse([kron(eye(k),ones(1,n));kron(eye(k),ones(1,n));Al;At;kron(ones(k,1)',eye(n))]);
        model.rhs = [floor(S(:));ceil(S(:));orderbs;bt(:);ones(n,1)];

        model.sense = cast(['>'*ones(k,1);'<'*ones(k,1);'<'*ones((k-1)*size(orderAX,1),1);'<'*ones(size(At,1),1)  ;'='*ones(n,1)],'char');


        %      model.A = sparse([kron(eye(k),ones(1,n));kron(eye(k),ones(1,n));kron(ones(k,1)',eye(n))]);
        %
        %     model.rhs = [floor(S(:));ceil(S(:));ones(n,1)];
        %     model.sense = cast(['>'*ones(k,1);'<'*ones(k,1);;'='*ones(n,1)],'char');
        % model.rhs = [(S(:));orderbs;bt(:);ones(n,1)];
        % model.sense = cast(['='*ones(k,1);'<'*ones((k-1)*size(orderAX,1),1);'<'*ones(size(At,1),1)  ;'='*ones(n,1)],'char');
    else
        %Original with all bells and whisles
        model.A = sparse([kron(speye(k),ones(1,n));Al;At;kron(ones(k,1)',speye(n))]);
        model.rhs = [n0*ones(k,1);orderbs;bt(:);ones(n,1)];

        model.sense = cast(['>'*ones(k,1);'<'*ones((k-1)*size(orderAX,1),1);'<'*ones(size(At,1),1)  ;'='*ones(n,1)],'char');

        %Removed At, bt
        %     model.A = sparse([kron(eye(k),ones(1,n));Al;kron(ones(k,1)',eye(n))]);
        % model.rhs = [n0*ones(k,1);orderbs;ones(n,1)];

        %model.sense = cast(['>'*ones(k,1);'<'*ones((k-1)*size(orderAX,1),1);'='*ones(n,1)],'char');

        %Removed also sorting
        %       model.A = sparse([kron(eye(k),ones(1,n));kron(ones(k,1)',eye(n))]);
        %   model.rhs = [n0*ones(k,1);ones(n,1)];

        %  model.sense = cast(['>'*ones(k,1);'='*ones(n,1)],'char');
    end

    model.lb = [zeros(k*n,1)];
    model.ub = [10*ones(k*n,1)];
    %model.vtype =  cast('b'*ones(k*n,1),'char');
    voronoi =false;
    if voronoi

        %   extending the objective and existing constraints
        model.obj = [model.obj; zeros(d*(k-1)*k,1)];
        model.A = [model.A spalloc(size(model.A,1), d*(k-1)*k,0)];


        Av = [];
        bv = [];
        %M = max(abs(X(:)))*d;
        M = n*100;;
        for jj = 1:k
            for mm = 1:k-1
                mm2 = mm;
                if mm>= jj
                    mm2 = mm+1;
                end

                leftzeros  = k*(k-1)*d-((jj-1)*(k-1)+(mm-1))*d-d;

                Anew = [(1+M)*kron(sparse(1,jj,1,1,k,1),speye(n)) sparse(n,((jj-1)*(k-1)+(mm-1))*d) -sparse(X)' sparse(n,leftzeros)];
                bnew = ones(n,1)*M;

                model.A = [model.A; Anew];
                model.rhs = [model.rhs;bnew];
                model.sense = [model.sense; cast('<'*ones(n,1),'char') ];

                Anew = [(1-M)*kron(sparse(1,mm2,1,1,k,1),speye(n)) sparse(n,((jj-1)*(k-1)+(mm-1))*d) -sparse(X)' sparse(n,leftzeros)];
                bnew = -ones(n,1)*M;
                model.A = [model.A; Anew];
                model.rhs = [model.rhs;bnew];
                model.sense = [model.sense; cast('>'*ones(n,1),'char') ];

            end
        end

        % model.vtype =  cast(['b'*ones(k*n,1);'c'*ones(d*k*(k-1),1)],'char');
        model.lb = [zeros(k*n,1);-1000000*ones(d*k*(k-1),1)];
        model.ub = [10*ones(k*n,1);1000000*ones(d*k*(k-1),1)];
    end


else
    model.A = sparse([kron(eye(k),ones(1,n));At;kron(ones(k,1)',eye(n))]);
    model.rhs = [n0*ones(k,1);bt(:);ones(n,1)];
    model.sense = cast(['>'*ones(k,1);'<'*ones(size(At,1),1);'='*ones(n,1)  ],'char');
    model.lb = [zeros(k*n,1)];
    model.ub = [10*ones(k*n,1)];


end


model.modelsense= 'min';
%model.vtype = cast(['B'*ones(n*k,1)],'char');
params.outputflag = 0;
params.OptimalityTol = 1e-9;
params.FeasibilityTol = 1e-9;
% if exist('interioronly','var')
% if interioronly
% params.method = 2;
% end
% end
result = gurobi(model, params);



if strcmp(result.status,'OPTIMAL')==1
    G = reshape(result.x(1:(k*n)),n,k);


    if true
        vbasis = find(result.vbasis == 0);
        cbasis = find(result.cbasis == 0);
        cbasis = cbasis+n*k;
        basis = [vbasis; cbasis];
    else
        basis = [];
    end
else
    G = [];
    basis = [];
end


% notbinary = (G>0 & G < 1);
% if sum(notbinary,"all") > 0
% disp('not binary')
% end


tsum = sum(diff(sum([X; ones(1,size(X,2))]*G)) < -orderb-sqrt(eps));
if tsum > 0
    disp('asdf')
end

%basis = find(G(:)> 0);

end







function [z,alphasum] = QPn(A,b,lb,ub,n,m,k)


clear model;

model.obj = [zeros(k*(m+1),1); ones(k,1)];
for j = 1:k
    model.quadcon(j).Qc = sparse([(m+1)*j],[(m+1)*k+j],[1],(m+2)*k,(m+2)*k);
    model.quadcon(j).Qc = model.quadcon(j).Qc-sparse((1:m)+(m+1)*(j-1),(1:m)+(m+1)*(j-1),ones(m,1),(m+2)*k,(m+2)*k);
    model.quadcon(j).rhs = 0;
    model.quadcon(j).q = zeros((m+2)*k,1);
    model.quadcon(j).sense = '<';
end

for ii = 1:m
    for j = 1:k-1
        model.quadcon(end+1).Qc = sparse([ ii+(j-1)*(m+1);(m+1)+(j)*(m+1)],[(m+1)+(j)*(m+1); ii+(j-1)*(m+1)],[1;-1],(m+2)*k,(m+2)*k);
        model.quadcon(end).rhs = 0;
        model.quadcon(end).q = zeros((m+2)*k,1);
        model.quadcon(end).sense = '<';
        %   ii+(j-1)*(m+1)
        %   (m+1)+(j)*(m+1)

        %  (m+1)+(j)*(m+1)
        % ii+(j-1)*(m+1)



    end
end
params.MIPGap =1e-2;
params.MIPGapAbs =1e-2;

model.Q = sparse(zeros((m+2)*k));

if size(A,1) > 0
    model.A = sparse([A zeros(size(A,1),k)]);
    model.rhs = b;
    model.sense = cast(['<'*ones(size(A,1),1)],'char');
else
    model.A = sparse(zeros(1,k*(m+2)));
    model.rhs = 0;
    model.sense = cast(['<'*ones(1,1)],'char');

end
%vtype  =  cast(['C'*ones(m,1); 'I'],'char');
%model.vtype = cast([kron(ones(k,1),vtype); 'C'*ones(k,1)],'char');
model.lb = [lb(:);-1000000*ones(k,1)];
model.ub = [ub(:);1000000*ones(k,1)];
model.modelsense= 'max';

params.NonConvex =2;
params.OptimalityTol = 1e-9;
params.FeasibilityTol = 1e-9;

params.outputflag = 0;
result = gurobi(model, params);
if strcmp(result.status,'OPTIMAL')==1
    z = reshape(result.x(1:((m+1)*k)),(m+1),k);
    alphasum= sum(result.x(((m+1)*k+1):end));
else
    z = [];
    alphasum = 0;
end


end





function Gamma = QPdirect(X,x0,n0,n,m,k)
clear model;

%if strcmp(minmax,'max')

%A1 = [kron(speye(k),ones(1,n)) zeros(k,(m+1)*k) ]; %>= n0
A1 = [kron(speye(k),ones(1,n)) spalloc(k,(m+1)*k,0) ]; %>= n0
%A2 = [kron(ones(1,k),speye(n)) zeros(n,(m+1)*k)]; %=1
A2 = [kron(ones(1,k),speye(n)) spalloc(n,(m+1)*k,0)]; %=1
A3 = [kron(speye(k),X) -speye((m+1)*k)];


model.A = sparse([A1;A2;A3]);
model.rhs = [n0*ones(k,1); ones(n,1); zeros(k*(m+1),1)];
model.sense= cast(['>'*ones(k,1);'='*ones(n+k*(m+1),1)],'char');

model.obj = [zeros(n*k,1);-2*x0];
qv = [zeros(n*k,1); ones(k*(m+1),1)];
model.Q = sparse(1:(n*k+k*(m+1)),1:(n*k+k*(m+1)),qv);

%params.NonConvex = 2;
model.lb = [zeros(n*k,1);-1000000000*ones(k*(m+1),1)];
model.ub = [ones(n*k,1);1000000000*ones(k*(m+1),1)];



model.modelsense= 'min';
params.outputflag = 0;
params.MIPGap =1e-2;
%params.MIPGapAbs =1e-20;
params.PoolSearchMode=1;
params.PoolSolutions = 50;
%params.OutputFlag = 1;
params.OptimalityTol = 1e-9;
params.FeasibilityTol = 1e-9;

result = gurobi(model, params);



if strcmp(result.status,'INFEASIBLE') ==1 || strcmp(result.status,'NUMERIC') ==1
    % disp("empty set")
    x = [];
    Gamma = [];
    obj = [];
    y = [];
    s = [];
    primalslack = [];
    opt = false;
    pool = [];
else
    x =  reshape(result.x(1:n*k),n,k);
    Gamma = reshape(x,n,k);

    obj = result.objval;
    opt = true;
    %    pool = result.pool;
    pool = [];

end
end


function mig = Assignment(X,centroid,normx)

Dcentroid= pdist2((diag(normx)*X)',(diag(normx)*centroid)',"squaredeuclidean");
mig = min(Dcentroid,[],2)==Dcentroid;

conflicting = [];
empties = find(sum(mig) == 0);

if ~isempty(empties)

    iter = 1;
    for k1 = 1:size(Dcentroid,2)

        for kk = 1:size(Dcentroid,2)
            em = kk;
            %The cost to flip
            additionalcost = Dcentroid( :,em)-sum(Dcentroid.*mig,2);
            additionalcost(conflicting) = inf;
            [addcost(kk,k1),newindex(kk,k1)] = min(additionalcost);

        end
        %Conflict?
        uniques = unique(newindex(:,iter));
        for kk = 1:length(uniques)
            if sum(newindex(:,iter) == uniques(kk)) > 1
                conflicting(end+1) = uniques(kk);
            end
        end

    end


    K = size(Dcentroid,2);
    %Find the least cost permutation for the empty
    [~,perm] =mexEMD(ones(K,1),ones(K,1),addcost);
    %[newindex perm]
    for k1 = 1:K %iteration
        idx = find(perm(:,k1)==1);
        mig(newindex(idx,k1),:) = 0;
        mig(newindex(idx,k1),idx) = 1;
    end

end
mig = double(mig);

%sums = sum([X;ones(1,size(X,2))]*mig);
%[~,si] = sort(sums);
%mig = mig(:,)

end
function y = vec(x)

y = x(:);
end





