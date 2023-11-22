%% NE in EPG and 2 strategy game
EPG_utility = zeros(ceil(1/q_t)+1 ,ceil(1/q_t)+1 , 2 , SV) ; 
utility = zeros(2 , 2 , 2 , SV) ; 
P = zeros(2 , SV , 1) ;
NE_final_rate = zeros(2  , SV) ;
NE_final_utility = zeros(2  , SV) ;
s1_star = zeros(SV , 1) ;
s2_star = zeros(SV , 1) ;
EPG_final_rate = zeros(2 , SV) ;
EPG_final_utility= zeros(2 , SV) ;

for i = 1 : SV
    utility(1 , 1 , 1 ,i) = rate(1  , 1 , 1 ,i) ; 
    utility(1 , 1 , 2 ,i) = rate(1 , 1 , 2 ,i) ; 
    utility(1 , 2 , 1 ,i) = rate(1 , end , 1 ,i) ;     
    utility(2 , 1 , 2 ,i) = rate(end , 1 , 2 ,i) ; 

    utility(1 , 2 , 2 ,i) = rate(1 , end , 2 ,i)-r * N_i ; 
    utility(2 , 1 , 1 ,i) = rate(end , 1 , 1 ,i)-r * N_i ; 

    utility(2 , 2 , 1 ,i) = rate(1 + floor(size(rate , 1)/2) ,1 + floor(size(rate , 2)/2) , 1 ,i)-r * N_i/2 ; 
    utility(2 , 2 , 2 ,i) = rate(1 + floor(size(rate , 1)/2) ,1 + floor(size(rate , 2)/2) , 2 ,i)-r * N_i/2 ;        
    P( :  ,i) = NE(utility( :  , : ,  :  , i)) 
    NE_final_rate(  :  , i) = P( 1 , i) * P( 2 , i) * rate(1  ,1  , :  , i) +...
                            P( 1 , i) * ( 1 - P(2 , i)) * rate( 1 , end , :  , i) + ...
                            (1 - P( 1 , i) ) * P( 2 , i) * rate( end , 1 ,  :  , i) + ...
                            ( 1 - P( 1 , i) ) * (1 - P( 2 , i) ) * rate( 1 + floor(size(rate , 1)/2) ,1 + floor(size(rate , 2)/2) ,  :  , i) ;  
    NE_final_utility(  : , i) = P( 1 , i) * P( 2 , i) * utility(1  ,1  , :  ,i) +...
                           P( 1 , i) * ( 1 - P(2 , i)) * utility( 1 , 2 , :  , i) + ...
                           (1 - P( 1 , i) ) * P( 2 , i) * utility( 2 , 1 ,  :  , i) + ...
                           ( 1 - P( 1 , i) ) * (1 - P( 2 , i) ) * utility( 2 , 2 ,  :  , i) ; 
                       
   % Random Resource allocation for IRS 
%     Random_final_rate(1 , i) = 0.5 * (rate(1 , end , 1 , i) + rate(2 , 1 , 1 , i)) ; 
%     Random_final_rate(2 , i) = 0.5 * (rate(2 ,1 ,2 , i) + rate(1 , 2 , 2 , i)) ; 
%     Random_final_utility(1 , i) = 0.5 * (utility(1 ,2 ,1 , i)+utility(2 ,1 ,1 , i)) ; 
%     Random_final_utility(2 , i) = 0.5 * (utility(2 ,1 ,2 , i)+utility(1 ,2 ,2 , i)) ; 
                       
    for s1 = 1 : ceil(1/q_t)+1
       for s2 = 1 : ceil(1/q_t)+1
           N_i1 = floor((s1-1) * N_i * q_t) ;     % No. of elements allocated to BS1
           N_i2 = floor((s2-1) * N_i * q_t) ;     % No. of elements allocated to BS2
           EPG_utility(s1 , s2 , 1 , i) = rate(s1 , s2 , 1 , i)-N_i1 * r ; 
           EPG_utility(s1 , s2 , 2 , i) = rate(s1 , s2 , 2 , i)-N_i2 * r ; 
       end            
    end
    %     SaveOut('EPG_utility.mat' , EPG_utility( :  ,  : , : ,i)) ; 
    F = sum(EPG_utility( : , : , : ,i),3) ; 
    maximum = max(max(F)) ; 
    [s1_star(i),s2_star(i)] = find(F == maximum) ;
    
    EPG_final_rate(  : , i) = rate(s1_star(i),s2_star(i), : ,i) ; 
    EPG_final_utility(  : , i) = EPG_utility(s1_star(i),s2_star(i), : ,i) ; 


toc
end
disp('EPG optimal strategy')
disp(s1_star' - 1)
disp(s2_star' - 1)
NE_final_utility = 8 * NE_final_utility ;
EPG_final_utility = 8 * EPG_final_utility ;
% saveas(gcf,['EPG rate for ', num2str(N_i),' IRS elements-N1= ', num2str(N1),'-N2=',num2str(N2),'.epsc'])
%% plot 2 Utility
figure
tiledlayout(1,2)
nexttile
% draw BS1
set(gca,'FontName','Times New Roman','FontSize',11,'LineWidth',1.5) ; 
plot(0,0,'^','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b','linewidth',1)
hold on
% draw BS2
plot(x2,y2,'^','MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','linewidth',1)
% draw IRS
plot(x_i,y_i,'gs','MarkerSize',30,'linewidth',1)
% draw BS1 users
plot(x(1),start,'o-','MarkerSize',60 * R1/10,'linewidth',0.1)
plot(x(1),d,'o-','MarkerSize',60 * R1/10,'linewidth',0.1)
% draw BS2 users
plot(x(N1+1),start,'o-','MarkerSize',60 * R1/10,'linewidth',0.1)
% configs
title('User Location')
xlabel('x')
ylabel('y')
xlim([-5 150]) ; 
ylim([-5 50]) ; 
legend('BS1','BS2','IRS','UE1','UE1','UE2') ; 
grid on
nexttile
% draw results
plot(start : step : (SV-1) * step+start , EPG_final_utility(1 , 1 : SV) , 'b-*' , 'linewidth' , 1)
hold on
% plot(start : step : (SV-1) * step+start , Random_final_utility(1 , 1 : i) , 'b' , 'linewidth' , 1)
plot(start : step : (SV-1) * step+start , EPG_final_utility(2 , 1 : SV) , 'r-*' ,  'linewidth' , 1)
% plot(start : step : (SV-1) * step+start , Random_final_utility(2 , 1 : i) , 'r' , 'linewidth' , 1)


plot(start : step : (SV-1) * step+start ,  NE_final_utility(1 , 1 : SV) , 'b-o' , 'linewidth' , 1)
plot(start : step : (SV-1) * step+start , NE_final_utility(2 , 1 : SV) , 'r-o' , 'linewidth' , 1)
legend('BS1 EPG' , 'BS2 EPG' ,  'BS1 NE' ,  'BS2 NE') ; 
% legend('BS1 EPG' , 'BS1 RANDOM'...
%      , 'BS2 EPG' , 'BS2 RANDOM') ; 
% configs
grid on
set(gcf , 'color' , 'w');   % Set the figure background to white
title(['N = ' , num2str(N_i) , ' , r = ' , num2str(r) ])
xlabel('y')
ylabel('Utility')
xlim([start max((SV-1) * step+start , start+step)]) ; 
ylim([0 10]) ; 
