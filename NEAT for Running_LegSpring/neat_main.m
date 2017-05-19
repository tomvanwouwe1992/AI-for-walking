function [] = neat_main(NEAT_PARAMS)
%%%%%%%%%%%%%%%%%%% Main NEAT file  (calls all other functions) (script file)

%% Neuro_Evolution_of_Augmenting_Topologies - NEAT
%% developed by Kenneth Stanley (kstanley@cs.utexas.edu) & Risto Miikkulainen (risto@cs.utexas.edu)
%% Coding by Christian Mayr (matlab_neat@web.de)
%% Adapted by Tom Van Wouwe


%list of parameters

%PARAMTERS OF MAIN FILE

maxgeneration=NEAT_PARAMS.max_generation; % maximum number of generations for generational loop

% if set to 1, will load population, generation, innovation_record and species_record from neatsave.mat at start of algorithm, this allows us to continue from where we were at when a certain interruption would occur
% if set to 0, algorithm will start with initial population, new species record and new innovation record, at generation=1 (default option)
load_flag=1;       % Choose for making new initial (random) population
save_flag=1;       % if set to 1, will save population, generation, innovation_record and species_record to neatsave.mat at every generation (default option)

average_number_non_disabled_connections=[];      % Will hold the average number over a whole population
average_number_hidden_nodes=[];                  % Will hold the average number over a whole population
max_overall_fitness=[];                          % Maximal fitness over the whole population

%PARAMETERS TO START INITIAL POPULATION WITH
population_size=NEAT_PARAMS.population_size;   % The population size is a constant over the whole algorithm. I.e. the number of genes within the population
number_input_nodes=NEAT_PARAMS.number_input_nodes;  % For now, we opt for three inputs. We will work on different options to fill these in (e.g.: position, velocity and acceleration of the COM in the first dimension)
number_output_nodes=NEAT_PARAMS.number_output_nodes; % Our network translates these inputs into muscle excitations. For now we have a simple model with three muscles.

vector_connected_input_nodes=1:number_input_nodes;  % Initially, connect all input nodes to all output nodes. We start from a fully developed feedforward network!


%SPECIATION PARAMETERS
% The following structure will contain various information on single
% species ( a species being a part of the population that is relatively similar, based on "evolutionary distance")
% This data will be used for fitness sharing, reproduction, and for visualisation purposes
species_record(1).ID=0;                   % Consecutive species ID's
species_record(1).number_individuals=0;   % Number of individuals in species
% Matrix will be 4 rows by the current number of generations existent columns and will contain (from top to bottom):
% -number of generation,
% -mean raw fitness within the species,
% -max raw fitness within the species,
% -index of individual in population which has produced max raw fitness
species_record(1).generation_record=[];

% "The number of excess and disjoint genes between two pairs of genomes is
% a natural measure for similarity. Another measure can be the difference
% between the weights when comparing matching genes(neurons). By making a
% linear combination of these three we find a distance between two
% structures."
% The parameters specified below are coefficients to calculate this
% distance.
speciation.c1=1.0; %Speciation parameters as published by Ken Stanley
speciation.c2=1.0;
speciation.c3=0.4;
speciation.target_number_of_species = floor(population_size/15);
% "Another aspect of allowing speciation is to make sure species don't grow
% too large. We introduce fitness sharing, to penalize large species within
% the population. We can calculate the distance between a member of the
% population and each other member of the population. If this distance is
% above a treshold we choose the sharing factor to be zero otherwise it is
% the distance itself. By dividing the fitness function by the sum of all
% the sharing distances we get explicit shared fitness."
speciation.threshold=3.1;

%REPRODUCTION PARAMETERS
% Stagnation
% We try to estimate when a species is starting to stagnate and thus needs
% a refocus. If the maximum fitness does not improve more than the treshold
% value for the chosen number of generations (eg 15) than we don't allow
% further reproduction. The champion of the species will be retained
% unchanged.
stagnation.threshold=1e-2;        % Threshold to judge if a species is in stagnation (max fitness of species varies below threshold) this threshold is of course dependent on your fitness function, if you have a fitness function which has a large spread, you might want to increase this threshold
stagnation.number_generation=maxgeneration + 1;  % If max fitness of species has stayed within stagnation.threshold in the last stagnation.number_generation generations, all its fitnesses will be reduced to 0, so it will die out
% Computation is done the following way: the absolute difference between the average max fitness of the last stagnation.number_generation generations and the max fitness of each of these generations is computed and compared to stagnation.threshold.
% if it stays within this threshold for the indicated number of generations, the species is eliminated

% Refocus
% If maximum overall fitness of population doesn't change within threshold for this number of generations, only the top two species are allowed to reproduce
refocus.threshold=1e-2;
refocus.number_generation=maxgeneration + 1;

% Initial setup
initial.kill_percentage=0.4; % The percentage OF EACH SPECIES which will be eliminated (lowest performing individuals)
initial.number_for_kill=5;   % The above percentage for eliminating individuals will only be used in species which have more individuals than min_number_for_kill
% Please note that whatever the above settings, the code always ensures that at least 2 individuals are kept to be able to cross over, or at least the single individual in a species with only one individual
initial.number_copy=3;       % Species which have equal or greater than number_copy individuals will have their best individual copied unchanged into the next generation. These species are thought of to be very good and for that they want to keep their best individual unchanged into the next generation.

%Selection (ranking and stochastic universal sampling)
selection.pressure=2; % Number between 1.1 and 2.0, determines selective pressure towards most fit individual of species

% Crossover
% Now how does this crossover work exactly?? So every time we
% take to genomes together we produce a random value between 0 and 1. If
% this value is larger than the percentage we will make a cross-over of
% these two for the next generation.
crossover.percentage=0.80;                 % Percentage governs the way (crossover vs mutation) in which new population will be composed from old population.  exception: species with just one individual can only use mutation
% The next two parameters will have an influence on the specifics of the
% crossover. If crossover has been selected, the interspecies probability
% governs the intra/interspecies parent composition being used. WE STILL HAVE TO FIND OUT DETAILS OF THIS. The multipoint probability governs for the
% standard-crossover in which way matching connection genes are inherited randomly from both parents.
crossover.probability_interspecies=0.00;
crossover.probability_multipoint=0.6;     % In the (1-crossover.probability_multipoint) cases, weights of the new connection genes are the mean of the corresponding parent genes

% Mutation
% After cross-over we will mutate. 3 possible mutations: add node, add
% link, mutation of weights. (Simgoid curve steepness is constant)
mutation.probability_add_node=0.02;        % Chance of adding a node
mutation.probability_add_connection=0.03;  % Chance of adding a connection
mutation.probability_recurrency=0.0;        % If we are in add_connection_mutation, this governs if a recurrent connection is allowed. Note: this will only activate if the random connection is a recurrent one, otherwise the connection is simply accepted. If no possible non-recurrent connections exist for the current node genes, then for e.g. a probability of 0.1, 9 times out of 10 no connection is added.
mutation.probability_mutate_weight=0.90;    % Chance of mutating the weights
mutation.weight_cap=4;                      % Weights will be restricted from -mutation.weight_cap to mutation.weight_cap. We need to make sure we can cover our possible output range but don't go over it too much.
mutation.weight_range=0.50;                  % Random distribution with width mutation.weight_range, centered on 0. mutation range of 5 will give random distribution from -2.5 to 2.5
mutation.probability_gene_reenabled=0.25;   % Probability of a connection gene being reenabled in offspring if it was inherited disabled


%%%%%%%%%%%%%%%%main algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if load_flag==0
    
    % CREATE INITIAL POPULATION
    % Set-up the initial population (remember that this will be a fully
    % developed network without a hidden layer! -->  this is given by the
    % "vector connected input nodes". This function returns two structures:
    % population and innovation record. See in initial_population.m for
    % detailed explanation.
    [population,innovation_record]=initial_population(population_size, number_input_nodes, number_output_nodes, vector_connected_input_nodes);
    % Tom: This initial population is created based on the number of in and output nodes. These are provided by our model, the other input variables are algorithm parameters.
    
    % DO SPECIATION ON OUR INITIAL POPULATION (i.e. divide in species)
    number_connections=(length(vector_connected_input_nodes)+1)*number_output_nodes;
    
    %
    population(1).species=1;                                         % Store the first individual of our population in "species 1"
    
    % Create a species reference matrix (abbreviated, only weights, since there are no topology differences in initial population).
    % Based on this information (that differs between individuals) we can calculate evolutionary distance.
    % For our initial population, evolutionary distance is only established
    % by a difference in connection weights (topologies are identical and
    % minimal). So we only need to store the weights.
    matrix_reference_individuals=population(1).connectiongenes(4,:);
    species_record(1).ID=1;
    species_record(1).number_individuals=1;   % For now our first species holds only one individual.
    
    %Loop through rest of individuals and either assign to existing species
    %or create new species and USE FIRST INDIVIDUAL OF NEW SPECIES AS REFERENCE
    % We will only have one species in our initial population, the code will
    % make clear why!
    for index_individual=2:size(population,2)
        assigned_existing_species_flag=0;          % Put the flag to see whether individual belongs to existing species as false
        new_species_flag=0;                        % Put flag whether to create new species on false
        index_species=1;
        
        % We now loop through the existing species,
        % If it is detected that either the individual belongs to the checked species or there is need for a new species (there is now species left to compare to) we stop the loop
        while assigned_existing_species_flag==0 && new_species_flag==0
            % Compute evolutionary distance, with initial population, we only need connection weights.
            % In our initial population this distance is the speciation
            % weight multiplied with the average difference between the
            % connection weights
            distance=speciation.c3*sum(abs(population(index_individual).connectiongenes(4,:)-matrix_reference_individuals(index_species,:)))/number_connections;
            if distance<speciation.threshold %If within threshold, assign to the existing species. The speciation treshold is 3. We will always be within treshold (max distance is 0.4*2)
                population(index_individual).species=index_species;
                assigned_existing_species_flag=1;
                species_record(index_species).number_individuals=species_record(index_species).number_individuals+1;
            end
            index_species=index_species+1;
            if index_species>size(matrix_reference_individuals,1) && assigned_existing_species_flag==0 %Outside of species references, must be new species
                new_species_flag=1;
            end
        end
        if new_species_flag==1 %check for new species, if it is, update the species_record and use individual as reference for new species
            population(index_individual).species=index_species;
            matrix_reference_individuals=[matrix_reference_individuals;population(index_individual).connectiongenes(4,:)];
            species_record(index_species).ID=index_species;
            species_record(index_species).number_individuals=1; %if number individuals in a species is zero, that species is extinct
        end
    end
    generation=1;
else % start with saved version of evolution
    load 'neatsave200'
end

%%%
%%%
%%%  Generational Loop
%%%
%%%
flag_solution=0;
maximal_fitness = [];
balance = [];
balance_end = [];
effort = [];
integration_times = [];
% We will create generation after generation as long as the maximum number
% of generations is not reached and as long as we don't reach our stop
% criterion (flag_solution)

while generation<maxgeneration && flag_solution==0
     
    % First thing to do is to evaluate our current population --> calculate fitnesses of individuals and store them in population(:).fitness
    % IMPORTANT reproduction assumes an (ALL POSITIVE!!!) evaluation function where a higher value means better fitness (in other words, the algorithm is geared towards maximizing a fitness function which can only assume values between 0 and +Inf)
    disp('start simulation of population');
    population=runner_simulation(population); % This function will return our current population with corresponding (raw) fitnesses for each individual
    disp('simulation finished - reproduce --> maximal fitness');
    disp(max([population(:).fitness]));
%     for i = 1:size(population,2)
%         population(i).fitness = population(i).fitness + 10^5;
%     end
    maximal_fitness = [maximal_fitness max([population(:).fitness])];
    if save_flag==1 % Backup copies of current generation
        save 'neatsave' population generation innovation_record species_record
    end
    
    eval(['save neatsave', num2str(generation), ' population generation innovation_record species_record']);
    [~,index] = max([population(:).fitness]);
    
%     ForwardIntegrateUsingGivenNetwork(population(index), inputdata)
%     % Compute mean and max raw fitnesses in each species and store in species_record.generation_record
%     max_fitnesses_current_generation=zeros(1,size(species_record,2));    % Create a vector that will hold the max fitness for each species
%     
    % We loop over each species to assign their mean fitness values and to check for stagnation
    for index_species=1:size(species_record,2)
        if species_record(index_species).number_individuals>0
            % Check within the selected species for the maximum fitness and
            % store the individual index
            [max_fitness,index_individual_max]=max(([population(:).species]==index_species).*[population(:).fitness]);
            % Calculate the mean fitness of the selected species. 
            mean_fitness=sum(([population(:).species]==index_species).*[population(:).fitness])/species_record(index_species).number_individuals;
            % Compute stagnation vector (last stagnation.number_generation-1 max fitnesses plus current fitness
            % If the species already exists for a self defined number of
            % generations we want to check if the species is not stagnating
            % (i.e. improvement is not rising)
            if size(species_record(index_species).generation_record,2)>stagnation.number_generation-2
                % Stagnation vector holds the max fitnesses of a species of the
                % last (stagnation.number_generation) generations, and as well
                % the newest max fitness
                
                stagnation_vector=[species_record(index_species).generation_record(3,size(species_record(index_species).generation_record,2)-stagnation.number_generation+2:size(species_record(index_species).generation_record,2)),max_fitness];   % Stagnation vector holds
                if max_fitness > 50
                if sum(abs(stagnation_vector-mean(stagnation_vector))<stagnation.threshold)==stagnation.number_generation % If the max fitnesses for all generations in the stagnation vector are within a range of the treshold we set the mean fitness of the species very small
                    mean_fitness=0.01; %set mean fitness to small value to eliminate species (cannot be set to 0, if only one species is present, we would have divide by zero in fitness sharing. anyways, with only one species present, we have to keep it)
                end
                end
                
            end
            species_record(index_species).generation_record=[species_record(index_species).generation_record,[generation;mean_fitness;max_fitness;index_individual_max]];
            max_fitnesses_current_generation(1,index_species)=max_fitness;
        end
    end
    
    %check for refocus -- top species
    [~,index_top_species]=max(max_fitnesses_current_generation);
    % Check how much generations the top species already exists, if is more
    % than the number of refocus generations we check whether a refocus is necessary
    if size(species_record(index_top_species).generation_record,2)>refocus.number_generation
        index1=size(species_record(index_top_species).generation_record,2)-refocus.number_generation; % Index of the first generation number of the species that is within refocus horizon 
        index2=size(species_record(index_top_species).generation_record,2);                           % Index of last generation
        if sum(abs(species_record(index_top_species).generation_record(3,index1:index2)-mean(species_record(index_top_species).generation_record(3,index1:index2)))<refocus.threshold)==refocus.number_generation % If for all generations within the refocus horizon the top species does not change enough we want to refocus
            [~,vector_cull]=sort(-max_fitnesses_current_generation); % We sort species max fitnesses - vector_cull holds the indices in the right order
            vector_cull=vector_cull(1,3:sum(max_fitnesses_current_generation>0));  % only hold certain species (the two best!) -- we don't want to reproduce from inferior species
            for index_species=1:size(vector_cull,2)
                index_cull=vector_cull(1,index_species);
                if max_fitness > 50
                species_record(index_cull).generation_record(2,size(species_record(index_cull).generation_record,2))=0.01; % Put their fitness very low
                end
            end
        end
    end

    % After having computed all necessary information regarding our current population, it is time to update this population to a next generation.
    if flag_solution==0
        %call reproduction function with parameters, current population and species record, returns new population, new species record and new innovation record
        [population,species_record,innovation_record]=reproduce(population, species_record, innovation_record, initial, selection, crossover, mutation, speciation, generation, population_size);
        
    end
    %increment generational counter
    generation=generation+1;
    
    if generation > 5000
        flag_solution=1;
    end
end
end