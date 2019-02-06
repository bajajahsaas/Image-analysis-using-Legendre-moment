%FUNCTION FOR TRANSLATION INVARIANT OF LEGENDRE MOMENTS
function s = Legendre_trans_fn(p,q,a,b,img)
 %p,q are order; a,b are translating factor

img=imtranslate(img,[a b]);
[nx,ny]=size(img);
img = double(img);
sum_F = sum(img(:));
sum_Fx = 0;
sum_Fy = 0;
for i = 0 : nx - 1
    for j = 0 : ny - 1
        x=2*i/(nx-1)-1;
        y=2*j/(ny-1)-1;
        sum_Fx = sum_Fx + x * img(i + 1,j + 1);
        sum_Fy = sum_Fy + y * img(i + 1,j + 1);
    end 
end

%%% 1st order moment
x_0 = double(sum_Fx / sum_F);
y_0 = double(sum_Fy / sum_F);

for x = 0 : nx - 1
    for n = 0 : p
        if(n == 0)
            M_LX(n + 1,x + 1) = 1; % p0(x)
        elseif(n == 1)
            M_LX(n + 1,x + 1) = 2*x/(nx-1)-1-x_0 ; %p1(x)
        else
            A = (2*n-1)*( 2*x/(nx-1)-1-x_0)/n;
            B = (n-1)/n;
            M_LX(n + 1,x + 1) = M_LX(n,x + 1) * A - M_LX(n-1,x + 1) * B;
% legendre polynomial pp(x)
        end 
    end
end

for y = 0 : ny - 1
    for n = 0 : q
        if(n == 0)
            M_LY(n + 1,y + 1) = 1;
        elseif(n == 1)
            M_LY(n + 1,y + 1) = 2*y/(ny-1)-1-y_0;
        else
            A = (2*n-1)*( 2*y/(ny-1)-1-y_0)/n;
            B = (n-1)/n;
            M_LY(n + 1,y + 1) = M_LY(n,y + 1) * A - M_LY(n-1,y + 1) * B;
% legendre polynomial pq(y)
        end 
    end
end

T=0;
for i = 1 : nx
    for j = 1 : ny
        T=T+(M_LX(p+1,i)*M_LY(q+1,j)*img(i,j));
    end 
end

s=T* ((2*p+1)*(2*q+1))/(nx*ny);