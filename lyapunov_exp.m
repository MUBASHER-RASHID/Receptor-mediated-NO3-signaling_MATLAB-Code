clear all
%close all
clc
%enter the steady state
ss= xlsread('new_lyapunov_data_1'); % xlsread('Temp_final.xls'); %
%syms x1 x2 x3 NI %primary variable
%a1=0.04; a2=0.1; a3=0.1;%parameters
%b1=1.2; b2=0.3; b3=0.6;
%k1=0.07; k11=0.01; k12=0.01; k2=0.01; 
%h1=0.5; h2=1.2; 
%NI=xlsread('NIT.xls');
%f1= (a1*NI*x1^h1)/(k1+x1^h1);
%f2= (a2*x1/(k11+x1))*1/(k1+x1^h1)-(b2*x2/(k12+x2));
%f3= (a3*((x1^h1/(k1+x1^h1))+1)*x2^h2/(k2+x2^h2))-b3*x3;
%f=[f1;f2;f3];
%x=[x1,x2,x3];
%J=jacobian(f,x)
%E = subs(J, [NI x1 x2 x3], [0.004 0.0017 0.0077 0.029]);
%jj=@(NI,x1,x2,x3) [J(1,:);J(2,:);J(3,:)];
n=length(ss(:,1));
for i=1:n
    % mm is the jacobian matrix evaluated at steady states
    mm=jacob(ss(i,1),ss(i,2),ss(i,3),ss(i,4));  % calls jacob function created in program "jacob"
    
    % Stability of the steady states is determined by the eigen values ...
    % of the jacobian matrix evaluated at steady states. if det(mm)>0 and
    % trace(mm)<0, then the steady states are stable.
    % Reference: Book by JJ Tyson_Biochemical Oscillators
    Eigen(i,:)=eig(mm);                    
    lamda(i,:)=max(real(eig(mm)));                                     
end
disp('calculated lamda: the local lypunov exp:')
lamda
%disp(['lyap_exp: ' num2str(lamda)])
%plot(ss(:,1), abs(lamda'), '.-')

