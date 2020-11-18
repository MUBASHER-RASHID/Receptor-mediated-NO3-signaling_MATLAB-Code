
function jj=jacob(NI,x1,x2,x3) % jacobian jj doesnot contain x3, so I have removed x3 in function

%%a1=0.04; a2=0.1; a3=0.1;b1=1.2; b2=0.3; b3=0.7; %(b3=0.6)
%k1=0.07; k11=0.01; k12=0.01; k2=0.01; h1=0.5; h2=1.2; 

%%
%syms NI x1 x2 x3;
%a1=0.02; a2=0.1; a3=0.1; b1=0; b2=0.1; b3=0.3; 
%k1 =0.15; k11=0.1; k12=0.1; k2 =0.3; h1 =0.5; h2 =2.6;

%%
%f1= (a1*NI*x1^h1/(k1+x1^h1))-(b1*x1*x3);
%f2= (a2*x1/(k11+x1))*1/(k1+x1^h1)-(b2*x2/(k12+x2));
%f3= (a3*((x1^h1/(k1+x1^h1))+1)*x2^h2/(k2+x2^h2))-b3*x3;
%f=[f1;f2;f3];
%x=[x1,x2,x3];
%jj=jacobian(f,x);

%%
%IST jacob is given by b3=0.6 and IIND by b3=0.7

%jj=[                                                        NI/(50*x1^(1/2)*(x1^(1/2) + 7/100)) - NI/(50*(x1^(1/2) + 7/100)^2) - (6*x3)/5,                                                                                                                                                               0, -(6*x1)/5;
 %1/(10*(x1^(1/2) + 7/100)*(x1 + 1/100)) - x1^(1/2)/(20*(x1^(1/2) + 7/100)^2*(x1 + 1/100)) - x1/(10*(x1^(1/2) + 7/100)*(x1 + 1/100)^2),                                                                                                                (3*x2)/(10*(x2 + 1/100)^2) - 3/(10*(x2 + 1/100)),         0;
  %                                -(x2^(6/5)*(1/(20*(x1^(1/2) + 7/100)^2) - 1/(20*x1^(1/2)*(x1^(1/2) + 7/100))))/(x2^(6/5) + 1/100), (6*x2^(1/5)*(x1^(1/2)/(10*(x1^(1/2) + 7/100)) + 1/10))/(5*(x2^(6/5) + 1/100)) - (6*x2^(7/5)*(x1^(1/2)/(10*(x1^(1/2) + 7/100)) + 1/10))/(5*(x2^(6/5) + 1/100)^2),      -3/5];

%%

%jj =[                                                       NI/(50*x1^(1/2)*(x1^(1/2) + 7/100)) - 1/(50*(x1^(1/2) + 7/100)^2) - (6*x3)/5,                                                                                                                                                               0, -(6*x1)/5;
% 1/(10*(x1^(1/2) + 7/100)*(x1 + 1/100)) - x1^(1/2)/(20*(x1^(1/2) + 7/100)^2*(x1 + 1/100)) - x1/(10*(x1^(1/2) + 7/100)*(x1 + 1/100)^2),                                                                                                                (3*x2)/(10*(x2 + 1/100)^2) - 3/(10*(x2 + 1/100)),         0;
%                                  -(x2^(6/5)*(1/(20*(x1^(1/2) + 7/100)^2) - 1/(20*x1^(1/2)*(x1^(1/2) + 7/100))))/(x2^(6/5) + 1/100), (6*x2^(1/5)*(x1^(1/2)/(10*(x1^(1/2) + 7/100)) + 1/10))/(5*(x2^(6/5) + 1/100)) - (6*x2^(7/5)*(x1^(1/2)/(10*(x1^(1/2) + 7/100)) + 1/10))/(5*(x2^(6/5) + 1/100)^2),     -7/10];

%% Below Jacobian is created from program #callNANADgate_modified6.m

jj = [                                                         NI/(100*x1^(1/2)*(x1^(1/2) + 3/20)) - NI/(100*(x1^(1/2) + 3/20)^2),                                                                                                                                                              0,     0;
 1/(10*(x1^(1/2) + 3/20)*(x1 + 1/10)) - x1^(1/2)/(20*(x1^(1/2) + 3/20)^2*(x1 + 1/10)) - x1/(10*(x1^(1/2) + 3/20)*(x1 + 1/10)^2),                                                                                                                       x2/(10*(x2 + 1/10)^2) - 1/(10*(x2 + 1/10)),     0;
                               -(x2^(13/5)*(1/(20*(x1^(1/2) + 3/20)^2) - 1/(20*x1^(1/2)*(x1^(1/2) + 3/20))))/(x2^(13/5) + 3/10), (13*x2^(8/5)*(x1^(1/2)/(10*(x1^(1/2) + 3/20)) + 1/10))/(5*(x2^(13/5) + 3/10)) - (13*x2^(21/5)*(x1^(1/2)/(10*(x1^(1/2) + 3/20)) + 1/10))/(5*(x2^(13/5) + 3/10)^2), -3/10];
 
return
