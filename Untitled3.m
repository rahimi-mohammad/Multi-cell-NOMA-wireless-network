for i=1:SV
    for s1=1:ceil(1/q_t)+1
        N_i1=floor((s1-1)*N_i*q_t)    % No. of elements allocated to BS1

       for s2=1:ceil(1/q_t)+1
        N_i2=floor((s2-1)*N_i*q_t)   % No. of elements allocated to BS2
       utility(s1,s2,1,i)=rate(s1,s2,1,i)-N_i1*r;
       utility(s1,s2,2,i)=rate(s1,s2,2,i)-N_i2*r;
       end            
    end
%     SaveOut('utility.mat',utility(:,:,:,i));
    F=sum(utility(:,:,:,i),3);
    maximum = max(max(F));
    [s1_star,s2_star]=find(F==maximum)
    NE_final_rate( :, i)=rate(s1_star,s2_star,:,i)
end
