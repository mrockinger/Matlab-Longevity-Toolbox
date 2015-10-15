function y=HP4_Shape(x,age,scal)
% implements Heligman and Polard version 1
A = x(1)*scal(1);
B = x(2)*scal(2);
C = x(3)*scal(3);
D = x(4)*scal(4);
E = x(5)*scal(5);
F = x(6)*scal(6);
G = x(7)*scal(7);
H = x(8)*scal(8);
K = x(9)*scal(9);

age = age(:);

XX = G.*(H.^(age.^K));

y = A.^( (age+B).^C ) + D.*exp(-E.*(log(age)-log(F)).^2) + XX./(1+XX);
