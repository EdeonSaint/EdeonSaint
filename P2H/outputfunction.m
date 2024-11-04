function [state, options,optchanged] = outputfunction(options,state,flag)
%displays the function eval value at each iteration
disp(state.FunEval);
optchanged = false;
switch flag
 case 'init'
        disp('Starting the algorithm');
    case {'iter','interrupt'}
        disp('Iterating ...')
    case 'done'
        disp('Performing final task');
end