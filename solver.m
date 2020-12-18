
function [u, u_prime] = solver(n_ddl, elem, H, G )

    % Known values
    u_known = ones(1,n_ddl)';
    for i = 1 : length(elem)
        u_known(i) = elem(i).val;
    end

    % Construction of the linear system
    A = zeros(n_ddl);
    b = zeros(n_ddl, 1);

    for i = 1:length(elem)
        if elem(i).type == 1 %Dirichlet
            A(:, i ) = -G(:, i);
            b = b + -H( :, i ) * u_known(i) ;

        elseif elem(i).type == 0 %Neumann
            A(:, i ) = H( :, i );
            b = b + G(:,i) * u_known(i) ;
        end
    end

    % Solve for unknown
    x = mldivide(A,b); 

    % Reconstruction of values u and u_prime
    u = zeros(1,length(elem))'; 
    u_prime = zeros(1,length(elem))'; 
    for i = 1:length(elem)
        if elem(i).type == 1
            u(i)  = u_known(i) ;
            u_prime(i)  = x(i) ;

        elseif elem(i).type == 0
            u_prime(i)  = u_known(i) ;
            u(i)  = x(i) ;
        end
    end

end