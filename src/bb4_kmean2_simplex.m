function [qptree] = bb4_kmean2_simplex(Aend,bend,lb,zsimplex,qptree,m,km,lowbound,removebelows,Z1bound,startactive,normx, orderSimplex,Simplexsubspace)


if ~exist('qptree','var') || isempty(qptree)

    qptree.EPs = [];
    qptree.lbbases = [];
    qptree.bases = [];
    qptree.objs = [];
    qptree.active = [];
    qptree.uses = [];
    qptree.closematrix = [];
    qptree.lb = [];
    qptree.zsimplex = [];
    qptree.A = [];
    qptree.b = [];

end

if isempty(qptree.EPs) && size(qptree.A,1) > 0
    return;
end

if isempty(qptree.EPs)
    [qptree.EPs,qptree.lbbases,qptree.objs,qptree.Alb,qptree.blb] = SimplexEps(lb,zsimplex,m,km,Z1bound,orderSimplex,Simplexsubspace,normx);
    qptree.active = cast(spalloc(size(qptree.objs,1),size(qptree.objs,2),100),'logical');
    qptree.active(startactive) = true;
    qptree.bases =spalloc(1,size(qptree.lbbases,2),0);
    qptree.lb =lb;
    qptree.zsimplex = zsimplex;
    qptree.uses = zeros(1,size(qptree.objs,2));

    qptree.bases = cast(qptree.bases,'logical');
    tmpvector = [qptree.lbbases;qptree.bases];

    qptree.closematrix = sparse((tmpvector'*tmpvector==(length(lb)-1)));

end



oldmaxobj = max(qptree.objs);
%Which extreme points are no longer feasible?
if ~isempty(Aend)
   

    qptree.A = [qptree.A;Aend];
    qptree.b = [qptree.b;bend];


    tval =sqrt(eps);% 1*sqrt(eps);%1000*sqrt(eps);



    bs = Aend*qptree.EPs;
    [~, maxindex] = max(qptree.objs);


    infeasiblesbit = (Aend*qptree.EPs-bend >tval);
    infeasibles= find (infeasiblesbit);


    if size(qptree.bases,1) < size(qptree.A,1)
        qptree.bases(size(qptree.A,1),1) = 0;
    end

    feasiblebit = (~infeasiblesbit) ;

% [qptree.objs;qptree.objs > lowbound; feasiblebit]
%     sum(qptree.objs(feasiblebit) > lowbound)
% if max(qptree.objs(feasiblebit)) < lowbound
%     disp('After this, all points are lower than the lower bound.')
% end


    qptree.uses(feasiblebit) = qptree.uses(feasiblebit)+1;

    %closematrix contains all points that are adjacent by size(A,2)-1
    %constraints. feasiblebit is a binary vector with feasible points. The
    %hadamardproduct is the feasible adjacent points to every point.

    feasibleadjacents= qptree.closematrix(:,infeasiblesbit)&feasiblebit';%for every column k denoting an extreme point, the 1s in column k denotes a feasible adjacent extreme point to point k


    %We reduce the indexing to that we only have the feasible vectors left,
    %the otherones will be removed soon anyway



    reducedfeasibleadjacents =qptree.closematrix(:,feasiblebit)';


    closematrix1columns = cell(length(infeasibles),1);

    newlbbases = cell(1,length(infeasibles));
    newbases = cell(1,length(infeasibles));
    newEPs = cell(1,length(infeasibles));
    newobjs = cell(1,length(infeasibles));
    newactives= cell(1,length(infeasibles));

    producedn = 0;

    mergers = [];
    onesmove = 0;
    problemswithls = 0;
    for tk = 1:length(infeasibles)
        k = infeasibles(tk);



        lookup = feasibleadjacents(:,tk);%was k but with an unreduced feasibleadjacents


        if sum(lookup) > 0

            if sum(lookup) == 1
                %disp('moveit instead')
                onesmove = onesmove+1;
            end

            lbbit= qptree.lbbases(:,k);
            Abit = qptree.bases(:,k);



            %Do them all at once
            Ls = (bend-Aend*qptree.EPs(:,lookup))./(Aend*(qptree.EPs(:,k)-qptree.EPs(:,lookup)));

            if sum(Ls < -sqrt(eps) | Ls > 1) > 0
                problemswithls = problemswithls+1;
            end

            Ls(Ls < 0) = 0;
            Ls(Ls > 1) = 1;

            if ~isempty(Ls)

                newEP = Ls.*qptree.EPs(:,k)+(1-Ls).*qptree.EPs(:,lookup);


                newlbbase = lbbit&qptree.lbbases(:,lookup);

                newAbits = Abit&qptree.bases(:,lookup);
                newAbits(end,:) = 1;

                tmp = reshape(newEP,m+1,km,[]).*[normx;1];

                newobj = squeeze(sum(tmp(1:m,:,:).^2./tmp(m+1,:,:),[1 2]))';

                newEPs{tk} =newEP;
                newlbbases{tk} = newlbbase;
                newbases{tk} =newAbits;
                newobjs{tk} = newobj;

                newactives{tk} = sparse(ones(length(newobj),1),1:length(newobj),qptree.active(k) | qptree.active(lookup));


                closematrix1columns{tk} = find(reducedfeasibleadjacents(:,k))';

                producedn = producedn + length(newobj);
            end

        end
    end


    if problemswithls > 0
        disp(['problemswithls ' num2str(problemswithls)]);
    end
    clear reducedfeasibleadjacents;


    %Remove the infeasibles
    qptree.EPs(:,infeasibles) = [];
    qptree.lbbases(:,infeasibles) = [];
    qptree.bases(:,infeasibles) = [];
    qptree.objs(infeasibles) = [];
    qptree.active(infeasibles) = [];
    qptree.uses(infeasibles) = [];
    qptree.closematrix(:,infeasibles) = [];
    qptree.closematrix(infeasibles,:) = [];


    origsz = length(qptree.objs);



    %collect closematrix1, adjacency from existing points to the new points
    closematrix1columnsv = [closematrix1columns{:}];



    %Save them in sequential order. Pre allocation to save time





    ZSM = 2000;


    newlbbases= newlbbases(~cellfun('isempty',newlbbases));
    newbases= newbases(~cellfun('isempty',newbases));
    newEPs= newEPs(~cellfun('isempty',newEPs));
    newobjs=newobjs(~cellfun('isempty',newobjs));
    newactives = newactives(~cellfun('isempty',newactives));

    lbbasesN = logical(cell2mat(newlbbases));
    basesN = logical(cell2mat(newbases));
    EPsN = cell2mat(newEPs);
    objsN = cell2mat(newobjs);
    activesN = cell2mat(newactives);

    %disp(['added ' num2str(length(objsN)) ' points' ]);

    qptree.lbbases = [ qptree.lbbases  lbbasesN];
    qptree.bases = [qptree.bases basesN];
    qptree.EPs = [qptree.EPs EPsN];
    qptree.objs = [qptree.objs objsN];
    qptree.active = [qptree.active activesN];

    if length(qptree.objs)> 0
        qptree.uses(length(qptree.objs)) = 0;
    end


    clear newlbbases;
    clear newbases;
    clear newEPs;
    clear newobjs;


    %pad with zeros if not the correct size
    if size(qptree.bases,2) <length(qptree.objs)
        qptree.bases(1,length(qptree.objs)) = false;
    end

    if size(qptree.lbbases,2) <length(qptree.objs)
        qptree.lbbases(1,length(qptree.objs)) = false;
    end




    tmpvector1 = [lbbasesN;basesN];
    clear lbbasesN;
    clear basesN;

    %Reducing the bits so that only check the ones that are interesting
    rbit = (sum(tmpvector1,2)>0);



    tmpvector1 = tmpvector1(rbit,:);


    %The original slow one
    % closematrix2o =((tmpvector1'*tmpvector1==(size(A,2)-1)));

global forceCPU;
if isempty(forceCPU)
    forceCPU = false;
end
    rungpu = true;
    if rungpu && size(tmpvector1,2) >2000 && ~forceCPU


        SZM = 20000;%20000

        tid1 = tic;
        t1k = 0;
        asg = cell(ceil(size(tmpvector1,2)/SZM));
        bsg = cell(ceil(size(tmpvector1,2)/SZM));

        if size(tmpvector1,2) > 0
            for t1k = 1:ceil(size(tmpvector1,2)/SZM)
                k1 = 1+(t1k-1)*SZM;


                z1 = min(size(tmpvector1,2),(k1+SZM-1));
                t2k = 0;
                t2kmax = ceil(k1/SZM);
                parfor t2k = 1:t2kmax
                    k2 = 1+(t2k-1)*SZM;

                    z2 = min(size(tmpvector1,2),(k2+SZM-1));
                    [c1g,c2g,~] = find(fastNequalGPU(tmpvector1(:,k1:z1),tmpvector1(:,k2:z2)));
                    asg{t1k,t2k} = c1g(:)'+k1-1;
                    bsg{t1k,t2k} = c2g(:)'+k2-1;
                end
            end
        end

        tafg = [asg{:}];
        tbfg = [bsg{:}];

        cm2a = [tafg  tbfg];
        cm2b = [tbfg  tafg];

        tid1 = toc(tid1);

    else
        %Same thing without gpu


        SZM = 10000;
        as = cell(ceil(size(tmpvector1,2)/SZM));
        bs = cell(ceil(size(tmpvector1,2)/SZM));

        t1k = 0;

        ast = cell(ceil(size(tmpvector1,2)/SZM),1);

        for t1 = 1:SZM:size(tmpvector1,2)
            t1k =t1k+1;
            len1 = SZM;
            if t1+SZM > size(tmpvector1,2)
                len1 = size(tmpvector1,2)-t1;
            end
            t2k = 0;
            ast{t1k} = tmpvector1(:,t1:t1+len1);

            for t2 = 1:SZM:t1
                t2k = t2k+1;
                len2 = SZM;
                if t2+SZM > size(tmpvector1,2)
                    len2 = size(tmpvector1,2)-t2;
                end


                [ta, tb, ~] =find( (ast{t1k}' * ast{t2k} ==(size(qptree.A,2)-1)  ) ); %-2 since I removed the last one that all of them have


                as{t1k,t2k} = ta(:)'+t1-1;
                bs{t1k,t2k} = tb(:)'+t2-1;


            end
        end
        taf = [as{:}];
        tbf = [bs{:}];



        cm2a = [taf  tbf];
        cm2b = [tbf  taf];

    end


    %The old one build on the old sparse matrix.
    [closematrixrows,closematrixcols,~] = find(qptree.closematrix);
    %Old to new
    cm1b = 1:length(closematrix1columnsv);
    cm1a = closematrix1columnsv;
    %new to new was made above

    % When I've removed properly I can use this instead
    closematrixrows = [closematrixrows' cm1a (cm1b+origsz ) (cm2a+origsz )];
    closematrixcols = [closematrixcols' (cm1b+origsz) cm1a (cm2b+origsz )];
    %Dont calculate it
    qptree.closematrix = sparse(closematrixrows,closematrixcols,true);

% disp('after')
%     [qptree.objs;qptree.objs > lowbound]
%     sum(qptree.objs > lowbound)
% if max(qptree.objs) < lowbound
%     disp('After this, all points are lower than the lower bound.')
% end

end




end



function [EPs,lbbases,ubbases,objs] = BoxEps(lb,ub,m,km)
EPs = zeros(length(lb),2^length(lb));
%lbbases = zeros(length(lb),2^length(lb));
%ubbases = zeros(length(lb),2^length(lb));

lbbases = spalloc(length(lb),2^length(lb),100);
ubbases = spalloc(length(lb),2^length(lb),100);


for k = 0:2^length(lb)-1


    bit= bitget(k,1:length(lb))';
    EPs(:,k+1) = (1-bit).*lb+bit.*ub;
    lbbases(:,k+1) = sparse(1-bit);
    ubbases(:,k+1) = sparse(bit);
end

%tmp = EPs.^2;

tmp = reshape(EPs,m+1,km,[]);
objs = squeeze(sum(tmp(1:m,:,:).^2./tmp(m+1,:,:),[1 2]))';



%tmp(end,:) = [];
%objs = sum(tmp(1:end,:),1);
lbbases = cast(lbbases,'logical');
ubbases = cast(ubbases,'logical');

%Just pack them together, its pointless to keep them separate
lbbases = [lbbases;ubbases];
ubbases =[];%spalloc(1,size(lbbases,2),0);

[s,I] = sort(objs,'descend');
lbbases = lbbases(:,I);
EPs = EPs(:,I);
objs = objs(I);


end




function [EPs,lbbases,objs,Alb,blb] = SimplexEps(lb,zsimplex,m,km,Z1bound,orderSimplex,SimplexSubspace,normx)

if false%sloppy simplex

    lbbases = spalloc(km*(m+1)+1,km*(m+1),100);
    EPs = zeros(length(lb),km*(m+1));
    EPs(:,1) = lb;
    lbbases(:,1) = sparse([ones(km*(m+1),1); 0]);

    Alb = -eye(km*(m+1));
    blb = -lb;
    Alb = [Alb;ones(1,km*(m+1))];
    blb(end+1) = sum(zsimplex);

    for k = 1:km*(m+1)

        %ep = lb;
        %ep(k) = 0;
        %ep(k) = sum(zsimplex)-sum(ep); %här korsar simplexen axeln och planet <ep-zsimplex,1> = 0

        ep = lb;

        ep(k) = ep(k)+sum(zsimplex)-sum(lb); %här korsar simplexen axeln och planet <(ep-lb)-(zsimplex-lb),1> = 0


        EPs(:,k+1) =ep;
        bb = ones(km*(m+1)+1,1);
        bb(k) = 0;
        lbbases(:,k+1) =  sparse(bb);

        %bit= bitget(k,1:length(lb))';
        %EPs(:,k+1) = (1-bit).*lb+bit.*ub;
        %lbbases(:,k+1) = sparse(1-bit);
        %ubbases(:,k+1) = sparse(bit);
    end


    tmp = reshape(EPs,m+1,km,[]);
    objs = squeeze(sum(tmp(1:m,:,:).^2./tmp(m+1,:,:),[1 2]))';

else %trimmed simplex
    %The flat subspace
    v0 = SimplexSubspace.v0;
    v1 = SimplexSubspace.v1;

    V = [v0;v1];
    normal = sum(v1);

    %The constraints are
    %<z-lb,V(i,:)> >= 0
    %<z-zsimplex,normal> <= 0

    %Make sure null(v0) is orientended in the right direction (has the right sign)

    negs = find(V*(zsimplex-lb) < -sqrt(eps));
    V(negs,:) = -V(negs,:);
    if normal*(lb-zsimplex) > 0
        normal = -normal;
    end


    Alb = [-V;normal];
    blb = [-V*lb;normal*zsimplex];

    %The lower bound base
    lbbases(:,1) = sparse([ones(km*(m+1),1); 0]);
    %v0 is always active and then we alter the normal
    EPs(:,1) = Alb(find(lbbases(:,1)),:)\blb(find(lbbases(:,1)));

    for k = 1:size(v1,1)
        bbase = [ ones(size(v1,1)+1,1)];
        bbase(k) = 0;
        bbase =cast( [ones(size(v0,1),1);bbase],'logical');

        ep = Alb(find(bbase),:)\blb(find(bbase));
        EPs(:,end+1) = ep;
        lbbases(:,end+1) = bbase;

    end


    if   orderSimplex



        for kk = 1:size(EPs,2)
            [~,sortI(kk,:)] = sort(sum(reshape(EPs(:,kk),m+1,km)));
        end


        dd = pdist2(sortI,sortI,"squaredeuclidean");
        sums = sum(dd> 0);
        [~,II] = min(sums);
        sortI = sortI(II,:);




        %Rotate the simplex so that the ordering is correct
        [~,sortI] = sort(sum(reshape(EPs(:,1),m+1,km)));

        indexes = reshape(1:(m+1)*km,m+1,km);
        indexes = indexes(:,sortI);
        indexes = indexes(:);
        EPs = EPs(indexes,:);
        Alb = Alb(:,indexes);
    end


    tmp = reshape(EPs,m+1,km,[]).*[normx;1];
    objs = squeeze(sum(tmp(1:m,:,:).^2./tmp(m+1,:,:),[1 2]))';



end




lbbases = cast(lbbases,'logical');

[~,I] = sort(objs,'descend');
lbbases = lbbases(:,I);
EPs = EPs(:,I);
objs = objs(I);


end

