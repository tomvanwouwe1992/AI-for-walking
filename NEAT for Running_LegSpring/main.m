% Tom Van Wouwe
global NEAT_PARAMS

NEAT_PARAMS.max_generation      = 5000;
NEAT_PARAMS.population_size     = 80;  
NEAT_PARAMS.number_input_nodes  = 15;  
NEAT_PARAMS.number_output_nodes = 18; 

compileCPP;

neat_main(NEAT_PARAMS);








