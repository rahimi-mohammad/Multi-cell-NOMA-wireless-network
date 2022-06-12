% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method identifies wether a two-player game is an Exact Potential Game.
% strategy
% Inputs:
    % utility            - Utility function of each BS,
% Outputs:
    % U                  - Acheivable Rate of the BS.
% ------------------------------------------------------------------------


function y=EPGInd(utility)
    Ns1=size(utility, 1)
    Ns2=size(utility, 2)
    epsilon=0.01;
    y='True';
    for s1=1:Ns1-1
        for s2=1:Ns2-1
            for t1=s1+1:Ns1
                for t2=s2+2:Ns2
                    if abs([1, -1, 1, -1]*[utility(s1,s2,2)-utility(s1,s2,1);...
                                           utility(t1,s2,2)-utility(t1,s2,1);...
                                           utility(t1,t2,2)-utility(t1,t2,1);...
                                           utility(s1,t2,2)-utility(s1,t2,1)]...
                                                                                )>epsilon
                       [1, 1, 1, 1]*[utility(s1,s2,2)-utility(s1,s2,1);...
                                           -utility(t1,s2,2)+utility(t1,s2,1);...
                                           utility(t1,t2,2)-utility(t1,t2,1);...
                                          - utility(s1,t2,2)+utility(s1,t2,1)]
                                           
                                       [utility(s1,s2,2)-utility(s1,s2,1);...
                                           utility(t1,s2,2)-utility(t1,s2,1);...
                                           utility(t1,t2,2)-utility(t1,t2,1);...
                                           utility(s1,t2,2)-utility(s1,t2,1)]
                                                                                
                        y='False';
 
                    end
                end
            end            
        end
    end
        

end