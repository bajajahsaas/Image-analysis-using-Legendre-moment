function T1=inLegendre_Cntral_Matrix_c(img)%%%%%% translation invariant LM
                                                    %%% max -- order
[nx,ny]=size(img);
img = double(img);

for x = 0 : nx - 1
    for n = 0 : nx-1
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
    for n = 0 : ny-1
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

%storing legendre moments
for m = 0 : nx-1
    for n= 0 : ny-1
        L(m+1,n+1)=0;
        for i = 1 : nx 
            for j = 1 : ny
                    L(m+1,n+1)=L(m+1,n+1)+(M_LX(m+1,i)*M_LY(n+1,j)*img(i,j));
            end
        end
        L(m+1,n+1)=L(m+1,n+1)* ((2*m+1)*(2*n+1))/(nx*ny);
    end
end




for x = 0 : nx-1
    for y = 0 : ny-1
        T1(x+1,y+1)=0;
        for i = 0 : nx-1
            for j = 0 : ny-1
                    T1(x+1,y+1)=T1(x+1,y+1)+(M_LX(i+1,x+1)*M_LY(j+1,y+1)*L(i+1,j+1));
            end
        end
    end
end
