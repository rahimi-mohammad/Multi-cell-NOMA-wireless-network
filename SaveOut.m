    function Y= SaveOut(file_name,new)
%         file_name = 'SaveBest.mat';                                     % Name File
%         if strcmp(flag,'iter')
%             ibest = state.Best(end);
%             ibest = find(state.Score == ibest,1,'last');
%             bestx = state.Population(ibest,:);
        previous = load(file_name);
        var = [previous.var; new];                                % Read Previous Results, Append New Value
        save(file_name, 'var')                                      % Write ‘Best Individual’ To File
%         end
%         changed = true;                                                 % Necessary For Code, Use  Appropriate Value
    end