%This matlab code is written by Abhishek Raghuwanshi(13720222)for Elec341 project at UBC
%The matlab script below is for each question individually based on question number put together. 
% To run it you would have to run individual sections of this code as given by question numbers.
 
%Part 1

%q1 is in pdf

%q2

% Define numerical values
K = 5;
k1 = 0.1;
k2 = 0.1;

% Define the transfer function T(s)
numerator = K*[1,1]; 
denominator = [1,3+K*k1+K*k2,K*k1+2*K*k2+2];
T = tf(numerator, denominator);

% Plot the step response
figure;
step(T);
title('elec341 project q2');
xlabel('Time');
ylabel('Amplitude');
grid on;

%q3

% Define numerical values
K = 5;
k1 = 0.1;
k2 = 0.1;

% Adjust K to satisfy zero steady-state error condition
K = 2 / (1- k1 - 2*k2);

%when k1<0,k2>0the graph has 1 full oscillation, when k1=0.1,k2=0.1, half oscillation that evens out at 1.

% Define the transfer function T(s) with adjusted K
numerator = K*[1,1];
 denominator = [1,3+K*k1+K*k2,K*k1+2*K*k2+2];
T = tf(numerator, denominator);

% Plot the step response
figure;
step(T);
title('elec341 project q3');
xlabel('Time');
ylabel('Amplitude');
grid on;


%q4

% Define numerical values
K = 2.85;
k1 = 0.1;
k2 = 0.1;

% Adjust K to satisfy zero steady-state error condition
K = 2 / (1- k1 - 2*k2);

%when k1<0,k2>0the graph has 1 full oscillation, when k1=0.1,k2=0.1, half oscillation that evens out at 1.

% Define the transfer function T(s) with adjusted K
numerator = K*[1,1];
 denominator = [1,3+K*k1+K*k2,K*k1+2*K*k2+2];
T = tf(numerator, denominator);
disp('K=');
disp(K);
disp('k1=');
disp(k1);
disp('k2=');
disp(k2);

poles = pole(T);
poles
% Values will be printed below.

%
%Part 2

%Q5 on pdf
%q6

 K1 = 1;
G_num=K1;
G_den=[1/2600,1/26,1,0]
G=tf(G_num,G_den);

H_num=1;
H_den=[0.04,1];
H=tf(H_num,H_den);

%closed loop transfer function:G*H
GH= series(G,H);
%open loop transfer function: G/(1+G*H)
T=feedback(G,H);

%Print the transfer function with K=1
disp("Final open loop transfer function is");
T

%simulate and plot unit step response:
figure;
step(T);
title('elec341 project q6');
xlabel('Time in sec');
ylabel('Response');
grid on;

%q7 on simulink

%q8

K1 = 1;
G_num=K1;
G_den=[1/2600,1/26,1,0]
G=tf(G_num,G_den);

H_num=1;
H_den=[0.04,1];
H=tf(H_num,H_den);

%closed loop transfer function:G*H
GH= series(G,H);
%open loop transfer function: G/(1+G*H)
T=feedback(G,H);

%Print the transfer function with K=1
disp("Final open loop transfer function is");
T
%need root locus plot of close loop system(GH)
rlocus(GH);

%q9

K1 = 9.25;%arrived at this by trial and error
G_num=K1;
G_den=[1/2600,1/26,1,0]
G=tf(G_num,G_den);

H_num=1;
H_den=[0.04,1];
H=tf(H_num,H_den);

%open loop transfer function:G*H
GH= series(G,H);
%closed loop transfer function: G/(1+G*H)
T=feedback(G,H);

%Print the transfer function with K=1
disp("Final open loop transfer function is");
T
%need root locus plot of close loop system(GH)
rlocus(T);

%q10
K1 = 9.25;
G_num=K1;
G_den=[1/2600,1/26,1,0]
G=tf(G_num,G_den);

H_num=1;
H_den=[0.04,1];
H=tf(H_num,H_den);

%open loop transfer function:G*H
GH= series(G,H);
%closed loop transfer function: G/(1+G*H)
T=feedback(G,H);

%Print the transfer function with K=1
disp("Final transfer function at K=9.25 is");
T
%unit step response:
step(T);