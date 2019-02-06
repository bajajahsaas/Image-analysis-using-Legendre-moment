%function to calculate features of an image
function d = feat(img,p,q)
%p,q are translating factor
a = Legendre_trans_fn(2,0,p,q,img);
b = Legendre_trans_fn(0,2,p,q,img);
c = Legendre_trans_fn(2,1,p,q,img);
e = Legendre_trans_fn(1,2,p,q,img);
f = Legendre_trans_fn(3,0,p,q,img);
g = Legendre_trans_fn(1,1,p,q,img);
d = [a b c e f g];
end