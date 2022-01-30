landa=0.7;
% Pd=[zeros(5,4)  ones(5,1)]
Pd=[1/3 2/3 ;0 1]
% Pd=[0.3 0.4 0.2 0.1 ;0.2 0.3 0.5 0;0.1 0 0.8 0.1;0.4 0 0 0.6]
rd=[8;-5]
% rd=[20*ones(4,1);0]
% inv(eye(2,2)-landa*Pd)
value=inv(eye(2,2)-landa*Pd)*rd

% Pd^100
