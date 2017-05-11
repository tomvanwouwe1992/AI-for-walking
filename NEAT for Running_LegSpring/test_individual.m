function population = test_individual(individual)    %(population, inputdata)
Path = pwd;
model_name = fullfile(Path,'gait9dof18musc.osim');
timevector = [0.0 15.0];
folder_dll = fullfile(pwd, '/../OpenSimInstall/bin');
old_PATH = getenv('PATH');
setenv('PATH', [folder_dll ';' old_PATH]);

NN_info = [ NEAT_PARAMS.number_input_nodes    NEAT_PARAMS.number_output_nodes   size(individual.nodegenes,2)-NEAT_PARAMS.number_input_nodes-NEAT_PARAMS.number_output_nodes-1    size(individual.connectiongenes,2) ]; % input - output - hiddden - connection

value = Integrate_Runner_NEAT_vFAST_PRINT(model_name,timevector, individual.nodegenes, individual.connectiongenes,NN_info);

fitness = value + 1; % To make sure we have a positive value


setenv('PATH', old_PATH);
end