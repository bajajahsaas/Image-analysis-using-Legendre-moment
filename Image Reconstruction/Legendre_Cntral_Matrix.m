function s = Legendre_Cntral_Matrix(img,p,q)%%%%%% translation invariant LM
                                                    %%% max -- order
[nx,ny]=size(img);
img = double(img);

for x = 0 : nx - 1
    for n = 0 : p
        if(n == 0)
            M_LX(n + 1,x + 1) = 1; % p0(x)
        elseif(n == 1)
            M_LX(n + 1,x + 1) = 2*x/(nx-1)-1 ; %p1(x)
        else 
			A = (2*n-1)*( 2*x/(nx-1)-1)/n;
			B = (n-1)/n;
			M_LX(n + 1,x + 1) = M_LX(n,x + 1) * A - M_LX(n-1,x + 1) * B;% legendre polynomial pp(x)
        end
    end
end

for y = 0 : ny - 1
    for n = 0 : q
        if(n == 0)
            M_LY(n + 1,y + 1) = 1;
        elseif(n == 1)
            M_LY(n + 1,y + 1) = 2*y/(ny-1)-1;
        else 
			A = (2*n-1)*( 2*y/(ny-1)-1)/n;
			B = (n-1)/n;
			M_LY(n + 1,y + 1) = M_LY(n,y + 1) * A - M_LY(n-1,y + 1) * B;% legendre polynomial pq(y)
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
