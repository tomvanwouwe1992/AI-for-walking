function z=ForwardIntegrationRunner(x,model_name,timevector, nodegenes, connectiongenes,NN_info,print)

    for i = 1:size(x,2)
       connectiongenes(4,i) = x(i);
    end
        
    z = 10-Integrate_Runner_NEAT_vFAST(model_name,timevector,nodegenes, connectiongenes,NN_info,print);

end