%%%%%%%%%%%%%%%% XOR experiment file  (contains experiment, receives genom, decodes it, evaluates it and returns raw fitnesses) (function)

%% Neuro_Evolution_of_Augmenting_Topologies - NEAT
%% developed by Kenneth Stanley (kstanley@cs.utexas.edu) & Risto Miikkulainen (risto@cs.utexas.edu)
%% Coding by Christian Mayr (matlab_neat@web.de)

function population = runner_simulation(population)    %(population, inputdata)
Path = pwd;
model_name = fullfile(Path,'gait9dof18musc.osim');
timevector = [0.0 15.0];
folder_dll = fullfile(pwd, '/../OpenSimInstall/bin');
%folder_dll='F:\u0113530\AI RUNNING\OpenSimInstall\bin';
old_PATH = getenv('PATH');
setenv('PATH', [folder_dll ';' old_PATH]);
print = 0;
parfor j = 1:length(population)
    
    individual = population(j);
   
        
    NN_info = [ NEAT_PARAMS.number_input_nodes    NEAT_PARAMS.number_output_nodes   size(individual.nodegenes,2)-NEAT_PARAMS.number_input_nodes-NEAT_PARAMS.number_output_nodes-1    size(individual.connectiongenes,2) ]; % input - output - hiddden - connection
 
    value = Integrate_Runner_NEAT_vFAST(model_name,timevector, individual.nodegenes, individual.connectiongenes,NN_info,print);
    
    fitness = value + 1; % To make sure we have a positive value
    if fitness < 0
        fitness = 0.01;
    end  
    population(j).fitness = fitness;
end
setenv('PATH', old_PATH);
end
