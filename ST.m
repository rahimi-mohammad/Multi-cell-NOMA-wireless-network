% clc  
% tic
% N_iter = 1 ;
% [ ~ , rate] = PotentialGame(M_t, M_r, N1, N2, N_i, P_T, x, y, z, x1, y1, z1,...
%                         x2, y2, z2, x_i, y_i, z_i, alpha_d, alpha_r, noise_power, N_iter,...
%                         R, r, q_t, interference) ;
% close all
i = 18 ;
% rate = rate(:,:,:,i) ;
%%
r_0 = 0.1 ;
StepSize = 0.05 ;
r = r_0 ;
N_i_0 = zeros(2, 1) ;
N_i_star = zeros(2, 1) ;
epsilon  = 0.0002 ;
utility = zeros(ceil(1/q_t) + 1,ceil(1/q_t) + 1,2) ;
EPG_utility = zeros(ceil(1/q_t)+1,ceil(1/q_t)+1,2) ; 
s1_star = 0 ;
s2_star = 0 ;
N_i11 = [] ;
N_i22 = [] ;
u11 = [] ;
u22 = [] ;
%%
while 1 == 1
    r = r + StepSize * 1 ;
    for s1 = 1:ceil(1/q_t)+1
        for s2 = 1:ceil(1/q_t)+1
            N_i1 = floor((s1-1) * N_i * q_t) ;     % No. of elements allocated to BS1
            N_i2 = floor((s2-1) * N_i * q_t) ;     % No. of elements allocated to BS2
            EPG_utility(s1, s2, 1) = N1*rate(s1, s2, 1, i) - N_i1 * r ; 
            EPG_utility(s1, s2, 2) = N2*rate(s1, s2, 2, i) - N_i2 * r ; 
        end            
    end
    F = sum(EPG_utility,3) ; 
    maximum = max(max(F)) ; 
    [s1_star,s2_star] = find(F == maximum) ;
    N_i_0 = N_i_star ;
    N_i_star = [floor((s1_star-1) * N_i * q_t) ; floor((s2_star-1) * N_i * q_t) ] ;     % No. of elements allocated to BS2
    N_i11 = [ N_i11 ; N_i_star(1) ] ;
    N_i22 = [ N_i22 ; N_i_star(2) ] ;
    u11 = [ u11 ; EPG_utility(s1_star, s2_star, 1) ] ;
    u22 = [ u22 ; EPG_utility(s1_star, s2_star, 2) ] ;
    %        r_0 = r ;
    %        r = r - StepSize * sign( N_i/2 - sum(N_i_star ,"all"))
    %        r = r + StepSize * sign(1 + sum(N_i_0 ,"all") - sum(N_i_star ,"all"))
    %        if sign(1 + sum(N_i_0 ,"all") - sum(N_i_star ,"all")) < 0
    %             StepSize= 0.5 * StepSize ;
    %        end
    %        abs(r * sum(N_i_star ,"all") - r_0 * sum(N_i_0 ,"all"))
    %        if abs(r * sum(N_i_star ,"all") - r_0 * sum(N_i_0 ,"all")) < epsilon
    %            break ;
    %        end
            if sum( N_i_star , "all" ) == 0
                break ;
            end
end

RIS_utility=(r_0 + StepSize * 1 : StepSize : r + epsilon)'.*(N_i11 + N_i22);
maximum = max(max(RIS_utility)) ;
maximum = find(RIS_utility==maximum)
maximum = r_0 + StepSize * maximum;
for s1 = 1:ceil(1/q_t)+1
        for s2 = 1:ceil(1/q_t)+1
            N_i1 = floor((s1-1) * N_i * q_t) ;     % No. of elements allocated to BS1
            N_i2 = floor((s2-1) * N_i * q_t) ;     % No. of elements allocated to BS2
            EPG_utility(s1, s2, 1) = rate(s1, s2, 1, i) - N_i1 * maximum ; 
            EPG_utility(s1, s2, 2) = rate(s1, s2, 2, i) - N_i2 * maximum ; 
        end            
end
F = sum(EPG_utility,3) ; 
maximum = max(max(F)) ; 
[s1_star,s2_star] = find(F == maximum) ;
BS11(i)=EPG_utility(s1_star, s2_star, 1);
BS22(i)=EPG_utility(s1_star, s2_star, 2);

toc
figure
tiledlayout(1,2)
nexttile
plot( r_0 + StepSize * 1 : StepSize : r + epsilon , N_i11 , 'r')
hold on
grid on
plot( r_0 + StepSize * 1 : StepSize : r + epsilon, N_i22 , 'b' )
xlabel('RIS Price')
ylabel('Used RIS elements')
legend( 'BS1 ', 'BS2 ' )
% plot( r_0 : StepSize : r , (N_i11 + N_i22) , 'r')
size((r_0 + StepSize * 1 : StepSize : r + epsilon)')
size(N_i22 )
nexttile
plot( r_0 + StepSize * 1 : StepSize : r + epsilon , (r_0 + StepSize * 1 : StepSize : r + epsilon)'.*(N_i11 + N_i22) , 'g-o')
hold on
grid on
plot(r_0 + StepSize * 1 : StepSize : r + epsilon,  u11(1 : end),'b-^' )
plot(r_0 + StepSize * 1 : StepSize : r + epsilon,  u22(1 : end) ,'r-*')
xlabel('RIS Price')
ylabel('Utility')
legend( 'RIS Utility', 'BS1 utility', 'BS2 utility' )
saveas(gcf,['RIS utility and elements used for y=', num2str(i)],'epsc')
% figure
% plot(N_i11 , u11  )
