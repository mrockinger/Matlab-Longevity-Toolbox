function   [k1tild,k2tild,k3tild,gztild]=stdize_CBD3(k1,k2,k3,gz,age,yearv)
% standardizes the parameters as to have essentiall no predictability in the cohort effects
k1=k1(:); k2=k2(:); k3=k3(:); gz=gz(:);


Nc=length(gz);

% prepare for OLS regression of gz on constant, cohortnbr  and cohortnbr^2
y = gz;
cohv = (1:Nc)'; % cohort number
x = [ones(Nc,1),cohv,cohv.^2] ;

res=ols(y,x);
phi1=res.beta(1);
phi2=res.beta(2);
phi3=res.beta(3);

N=size(age,1);
T=length(yearv);
tv=(1:T)';

gztild = gz - phi1 - phi2*cohv - phi3*cohv.^2; % residual
plot(gztild)

k3tild = k3 + phi3;


xb     = mean(age);
sigx2  = var(age);

k2tild = k2 - phi2 + 2*phi3*(xb-tv);

k1tild = k1 + phi1 + phi2*(tv-xb) + phi3*(tv.^2 - 2*tv.*xb + xb^2+sigx2);

% these are meant to be row vectors
k1tild = k1tild';
k2tild = k2tild';
k3tild = k3tild';


