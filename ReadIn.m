% ------------------------------------------------------------------------ 
%  mrahimi7755 - Sharif University of Technology, Iran
% ------------------------------------------------------------------------
% Rate.m - This method saves the the given input variable into file_name
% .mat file
% Inputs:
    %file_name           - Name of the existing file to which you want to add your data ,
    % new                - The variable you want to write in the first file.
% Outputs:
% ------------------------------------------------------------------------    

function var=ReadIn(file_name)

for i=1:SV
    sample=load(file_name);
    for j=1:2
        start=1+(i-1)*7;
        rate(:,:,j,i)=sample.var(start:start+5,1:end-1,j)
        
        
    end
    
    
end
end

    