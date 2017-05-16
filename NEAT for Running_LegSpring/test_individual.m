function fitness = test_individual(individual)    %(population, inputdata)
NEAT_PARAMS.max_generation      = 5000;
NEAT_PARAMS.population_size     = 120;  
NEAT_PARAMS.number_input_nodes  = 15;  
NEAT_PARAMS.number_output_nodes = 18; 
Path = pwd;
model_name = fullfile(Path,'gait9dof18musc.osim');
timevector = [0.0 15.0];
folder_dll = fullfile(pwd, '/../OpenSimInstall/bin');
old_PATH = getenv('PATH');
setenv('PATH', [folder_dll ';' old_PATH]);
print  = 1;
NN_info = [ NEAT_PARAMS.number_input_nodes    NEAT_PARAMS.number_output_nodes   size(individual.nodegenes,2)-NEAT_PARAMS.number_input_nodes-NEAT_PARAMS.number_output_nodes-1    size(individual.connectiongenes,2) ]; % input - output - hiddden - connection

fitness = Integrate_Runner_NEAT_vFAST(model_name,timevector, individual.nodegenes, individual.connectiongenes,NN_info, print);

setenv('PATH', old_PATH);
end