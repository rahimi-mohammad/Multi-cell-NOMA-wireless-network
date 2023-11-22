% This script calculates the solution of the max min problem by exhastive
% search and compares it to the soltion provided by SVD method.
%% freeze
        % N1=2;   % No. of users
        % N=2;   % No. of IRS elements
        % h_d=zeros(N1,1);
        % h_r=zeros(N,N1);
        % g=zeros(N,1);
        % load('h_d.mat','h_d');
        % load('g.mat','g');
        % load('h_r.mat','h_r');
%% exhaustive search
% M=5;
% channel_gain=0;
% opt_theta=zeros(N,1);
% min_theta=zeros(N,1);
% min_channel_gain=0;
% for i=1:N1
%    theta=zeros(N,1);
%    theta(N)=1;
%     while sum(theta)~=0
%         for k=1:2*M
%             theta(N)=k;
% %             calculate channel gain
%             if abs(h_d(i)+transpose(conj(h_r(:,i)))*diag(exp(1j*theta*pi/(2*M)))*g)>channel_gain
%                opt_theta=theta; 
%                channel_gain=abs(h_d(i)+transpose(conj(h_r(:,i)))*diag(exp(1j*theta*pi/(2*M)))*g);
%             end
%             
%         end
%         theta;
%         for m=1:N-1
%             theta(N-m)=theta(N-m)+floor(theta(N+1-m)/(2*M));
%             theta(N+1-m)=rem(theta(N+1-m),2*M);
%         end
%             theta(1)=rem(theta(1),2*M);
%         theta
%     end 
% %     theta=theta*pi/10;
%     if (channel_gain<min_channel_gain || i==1)
%         disp('here')
%         i
%         min_channel_gain=channel_gain;
%         min_theta=opt_theta;
%     end
%     theta;
% end
% exp(1j*min_theta*pi/(2*M))
