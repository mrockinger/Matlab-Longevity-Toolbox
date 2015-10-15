function d=HP1(x,age,qx,scal)
% implements Heligman and Polard version 1

age=age(:);
qx = qx(:);
px = 1-qx;

LHS = qx./px;

qxmodel = HP1_Shape(x,age,scal);

RHS = qxmodel./(1-qxmodel);

d = sum( ( log(LHS)-log(RHS) ).^2 );