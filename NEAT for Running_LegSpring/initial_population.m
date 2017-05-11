%%%%%%%%%%%%%%%%%%%%%%% generate minimal initial population 

%% Neuro_Evolution_of_Augmenting_Topologies - NEAT 
%% developed by Kenneth Stanley (kstanley@cs.utexas.edu) & Risto Miikkulainen (risto@cs.utexas.edu)
%% Coding by Christian Mayr (matlab_neat@web.de)

function [population,innovation_record]=initial_population(number_individuals,number_input_nodes,number_output_nodes,vector_connected_input_nodes)

% ABOUT OUR POPULATION AND ITS ENCODING
% We encode our NN in NEAT in a very specific manner (that may look redundant), which will appear very useful for the implementation.
% We encode our NN via Node -and Connectiongenes.

% NODEGENES
% For each node (column) we store the following information in a cell (row)
% - Node ID
% - Node type (input, output, hidden, bias)
% - Node input state (enabled/disabled) - for a bias node always enabled (1)
% - Node output state (enabled/disabled) - for a bias node always enabled (1)
% The number of columns tells us immediately a lot about the complexity of
% our topology since it tells us the number of nodes (= input + output +
% hidden + bias(1)). For our current problem this is 7!

% CONNECTIONGENES
% For each connection (column) we store the following information in a cell (row)
% - Innovation number
% - Connection from

% - Connection to
% - Weight
% - Enabled?? y/n

% For each member in our population we also store the fitness and the
% species ID to which it belongs!


% ABOUT INNOVATION
% An innovation can be two things --> either a new link or a new neuron! A
% link comes in between two nodes, a new node comes in between two nodes as
% well and inheritely asks for the introduction of two new links (that's
% three innovations at once!).
% --> Innovation concerns with topology! (mutation of weights, cross over is not innovation)
% "Innovation_record" tracks innovations in a 5 rows by (number of innovations) columns matrix. Every column holds the information of a specific innovation.
% The first row contains innovation number, this serves as identification
% The second row tells which node the innovation connects from --> connect_from_node 
% The third row tells which node the innovation connects to --> connect_to_node 
% The fourth row holds the node (neuron) number if a new node is introduced as mutation.
% In the fifth row the generation at which the innovation occured is stored(generation is assumed to be zero for the innovations in the initial population)
% We only need to store connection genes as innovation, as a new node (innovation) always goes with new connections!
   
% We have a fully developed feedforward network with no hidden layers so we
% can easily calculate our number of connections. Note that we add one node
% to the input nodes (bias node)
number_connections=(length(vector_connected_input_nodes)+1)*number_output_nodes; 

% This will initialize the "connection from" row of our connectiongenes!
vector_connection_from=rep([vector_connected_input_nodes,number_input_nodes+1],[1 number_output_nodes]); % [1 2 3 4 1 2 3 4 1 2 3 4]

%Initializing the next row (connection to is a (little) bit more
%difficult). In this way it is done in a very general way, so that if we
%change our problem (and input/output structure we don't have to rewrite
%all of this)
vector_connection_to=[];
for index_output_node=(number_input_nodes+2):(number_input_nodes+1+number_output_nodes)   % outputnode index starts at number of input nodes + 1 (biasnode) + 1 (just to take the next free number)!
   vector_connection_to=[vector_connection_to,index_output_node*ones(1,length(vector_connected_input_nodes)+1)]; % [ 5 5 5 5  6 6 6 6 7 7 7 7]
end
connection_matrix=[vector_connection_from;
                   vector_connection_to];
% All of the above is equal for all members of the initial population,
% because it was concerning topology. For the weights, and state of the
% connections there is some random number chosen for each individual in the initial population.      
% Loop over individuals

for index_individual=1:number_individuals 
   % We specify our node information
   population(index_individual).nodegenes=[1:(number_input_nodes+1+number_output_nodes);                  % Node IDs --> 1 by 1 to our number of nodes
                                           ones(1,number_input_nodes),4,2*ones(1,number_output_nodes);    % Type: first 3 are input (type 1), 4th is bias(type 4), 5th to last is output (type 2).
                                           zeros(1,number_input_nodes),1,zeros(1,number_output_nodes);    % input state all zero, except the bias node which is enabled
                                           zeros(1,number_input_nodes+1+number_output_nodes)];            % output state is similar
   % We specify connection information                                    
   population(index_individual).connectiongenes=[1:number_connections;                                    % Innovation number
                                                 connection_matrix;                                       % from - to nodes have been already worked out above based on the basic structure
                                                 rand(1,number_connections)*2-1;                          % Weights are random numbers that in between -1 to 1.
                                                 ones(1,number_connections)];                             % All connections are enabled
   population(index_individual).fitness=0;                                                                % Fitness is set to 0
   population(index_individual).species=0;                                                                % Species to which our individual/member belongs to is 0.
end

% Finally we store in our innovation record the following information
% Innovation number, from, to --> directly from the connectiongenes matrix
% The highest node ID is the highest node ID in the current population (no
% hidden nodes yet), the generation number is 0 as well.
innovation_record=[population(index_individual).connectiongenes(1:3,:);    % Take innovation number and from and to connections from the connectiongenes
                        zeros(size(population(index_individual).connectiongenes(1:2,:)))];    
innovation_record(4,size(innovation_record,2))=max(population(1).nodegenes(1,:)); % replace 4th row with highest node ID for initial population. WAAROM ENKEL LAATSTE INNOVATIE DIE HIGHEST NODE ID KRIJGT???