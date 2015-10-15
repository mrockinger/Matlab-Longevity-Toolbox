function d=HP2(x,age,qx,scal)
% implements Heligman and Polard version 3

age=age(:);
qx=qx(:);

RHS = HP2_Shape(x,age,scal);

%d = sum( (RHS./qx - 1.0).^2 );

d = sum( (log(RHS)-log(qx)).^2 );