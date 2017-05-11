function population = test_individual(individual)    %(population, inputdata)
Path = pwd;
model_name = fullfile(Path,'gait9dof18musc.osim');
timevector = [0.0 15.0];
folder_dll='F:\u0113530\AI RUNNING\OpenSimInstall\bin';
old_PATH = getenv('PATH');
setenv('PATH', [folder_dll ';' old_PATH]);

NN_info = [ 14    18   size(individual.nodegenes,2)-14-18-1    size(individual.connectiongenes,2) ]; % input - output - hiddden - connection

value = Integrate_Runner_NEAT_vFAST_PRINT(model_name,timevector, individual.nodegenes, individual.connectiongenes,NN_info);

fitness = value + 1; % To make sure we have a positive value


setenv('PATH', old_PATH);
end