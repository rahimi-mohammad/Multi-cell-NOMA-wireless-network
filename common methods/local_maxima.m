function maxima = local_maxima(vector , epsilon)
% This function finds the local maxima in a vector

% Initialize variables
maxima = [];
% Find the indices of all local maxima
if vector(1)== max(vector) 
    maxima = [maxima vector(1)] ;
end
for i = 2:length(vector)-1
    if vector(i) >= vector(i-1) && vector(i) >= vector(i+1)
        maxima = [maxima vector(i)] ;
    end
end
% display(['Local maxima ' num2str(maxima)])
for i = 1 : length(maxima)
    if max(maxima) <= (1 + epsilon) * maxima(i)
        maxima = maxima(i) ;
        break
    end
end
% Return the values of the local maxima
maxima = max(maxima);
% display(['Local maxima ' num2str(maxima)])
end