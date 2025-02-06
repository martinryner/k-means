function C = fastNequalGPU(vec1,vec2)
num_rows=size(vec1,2);
num_cols=size(vec2,2);
nn = ceil(size(vec1,1)/32);
offs = 0;
for n1=1:nn
len = 32;
if offs+len > size(vec1,1)
    len = size(vec1,1)-offs;
end
    uvec1(n1,:) = uint32(full(bit2int(vec1(offs+ (1:len),:),len)));
    uvec2(n1,:) = uint32(full(bit2int(vec2(offs+ (1:len),:),len)));
   offs = offs+32;
end

num_i=size(uvec1,1);
C(num_rows,num_cols) =false;
Range=[num_rows,num_cols,1];
settings = sprintf('-DNR=%d -DNC=%d -DNI=%d -DREAL=double',num_rows,num_cols,num_i);


[run_time]=cl_run_kernel(1,'./src/mul_kernel.cl',settings,'MM',Range,0,uvec1,uvec2,C,[1 1 2]);



