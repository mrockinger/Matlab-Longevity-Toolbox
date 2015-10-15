function   [k1tild,k2tild,gztild]=stdize_CBD2(k1,k2,gz,age,yearv)
% standardizes the parameters as to have essentiall no predictability in the cohort effects
k1=k1(:); k2=k2(:); gz=gz(:);

Nc=length(gz);

% prepare for OLS regression of gz on constant and time 
y=gz;
cohv=(1:Nc)'; % cohort number
x=[ones(Nc,1),cohv];
res=ols(y,x);
phi1=res.beta(1);
phi2=res.beta(2);

gztild=gz-phi1-phi2*cohv; % residual

k2tild = k2-phi2;

T=length(yearv);
k1tild = k1 + phi1 + phi2*(1:T)'-phi2*mean(age);

% these are meant to be row vectors
k1tild = k1tild';
k2tild = k2tild';


