clc;clear;
% lambda2 > lambda1
%boundary of convergence speed

%%% parameter %%%
lambda1=0.01;
lambda2=2;
%
wd_min = -1;
wd_max = 1;
vd_min = 0.05;
vd_max = 0.5 ;
%%% parameter %%%

A0 = [0,wd_min,0 ; -wd_min,0,vd_min ; 0,0,0];
Aalpha = [0,wd_max-wd_min,0 ; -(wd_max-wd_min),0,0 ; 0,0,0];
Abeta = [0,0,0 ; 0,0,vd_max-vd_min ; 0,0,0];
B = [1,0;0,0;0,1];

A1 = A0;
A2 = A0+Aalpha;
A3 = A0+Abeta;
A4 = A0+Aalpha+Abeta;

% LMI practice
% Initialization
setlmis([]) 

% Set matrix variables.
% parameter_name = lmivar(type,struct) 
% For more information about "lmivar", please refer to https://www.mathworks.com/help/robust/ref/lmivar.html. 
% Variable P is a positive definite matrix, P = P'. W = P^{-1} is also a
% symmetric matrix.
W = lmivar(1,[3 1]);
% K1 and K2 are 2x3 vectors.
% Y = K*W .
Y = lmivar(2,[2 3]);  
% Setup LMIs
% lmiterm(termID,A,B,flag)
% (A_k W+B_k Y)+(WA_k^T+Y^T B_k^T )+2λW

lmiterm([1 1 1 W], A1, 1, 's')  % A*W+W*A'
lmiterm([1 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([1 1 1 W], 2*lambda1,1)  % 2*lampda*W

lmiterm([2 1 1 W], A2, 1, 's')  % A*W+W*A'
lmiterm([2 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([2 1 1 W], 2*lambda1,1)  % 2*lampda*W

lmiterm([3 1 1 W], A3, 1, 's')  % A*W+W*A'
lmiterm([3 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([3 1 1 W], 2*lambda1,1)  % 2*lampda*W

lmiterm([4 1 1 W], A4, 1, 's')  % A*W+W*A'
lmiterm([4 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([4 1 1 W], 2*lambda1,1)  % 2*lampda*W

lmiterm([-5 1 1 W],1,1)
%
lmiterm([-6 1 1 W], A1, 1, 's')  % A*W+W*A'
lmiterm([-6 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([-6 1 1 W], 2*lambda2,1)  % 2*lampda*W

lmiterm([-7 1 1 W], A2, 1, 's')  % A*W+W*A'
lmiterm([-7 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([-7 1 1 W], 2*lambda2,1)  % 2*lampda*W

lmiterm([-8 1 1 W], A3, 1, 's')  % A*W+W*A'
lmiterm([-8 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([-8 1 1 W], 2*lambda2,1)  % 2*lampda*W

lmiterm([-9 1 1 W], A4, 1, 's')  % A*W+W*A'
lmiterm([-9 1 1 Y], B, 1, 's')   % B*Y+Y'*B'
lmiterm([-9 1 1 W], 2*lambda2,1)  % 2*lampda*W
LMISYS = getlmis;

[tmin,xfeas] = feasp(LMISYS);

W_ = dec2mat(LMISYS,xfeas,W)
Y_ = dec2mat(LMISYS,xfeas,Y)
K = Y_*inv(W_)
is_feasiable = tmin < 0