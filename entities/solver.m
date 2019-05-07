classdef solver
    methods(Static)
        
        function vec = get_node_index_to_equation_index_mapping(...
                number_of_nodes, nodes_with_direchlet_boundary_conditions)
            
            vec = -ones(number_of_nodes,1);
            cnt = 1;
            
            for k = 1 : number_of_nodes
                if ~ismember(k, nodes_with_direchlet_boundary_conditions)
                    vec(k) = cnt;
                    cnt = cnt + 1;
                end
            end
        end
    end
end