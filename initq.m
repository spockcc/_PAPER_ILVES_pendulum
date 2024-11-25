function [q, flag, it, his, res]=initq(q0, param, tau, maxit)

% Finds a solution to the constraint equation 
%
%     g(q) = 0
%
% using Newton's method. There are m constraints and n atoms
% 
% INPUT:
%   q0      initial guess for q
%   param   structure describing the constraints
%             param.dim   the dimension of the vectorspace
%             param.m   number of constraints
%             param.n   number of atoms
%             param.dg  the Jacobian of the constraint function
%  tau      terminate if norm(g(q),'inf') <= tau
%  maxit    maximum number of Newton steps
%
% OUTPUT:
%   q      the final approximation of the solution
%   flag   if norm(g(q),'inf') <= tau, then flag = 1 else flag = 0.
%   it     the number of Newton steps completed
%   his    array such that his(:,j) is the jth approximation
%   res    res(j)=norm(g(q(:,j)),'inf')
% 
% MINIMAL WORKING EXAMPLE: initq_mwe1.m

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2022-02-28  Initial programming and testing
%   2022-08-03  Added support for finite dimensional vector spaces

% Isolate the dimension of the vector space
dim=param.dim;

% Isolate the number of atoms
n=param.n;

% Isolate the constraint function
g=param.g;

% Isolate the Jacobian of the constraint function
dg=param.dg;

% Allocate space for his and res
his=zeros(dim*n,maxit); 
res=zeros(maxit+1,1);

% Initialize q
q=q0;

% Evaluate the initial value 
aux=g(q);

% Initialize his and res
his(:,1)=q; res(1)=norm(aux,'inf');

% Assume failure
flag=0;

% Initialize the number of Newton steps
it=0;

if (res(1)<=tau)
    % The constraint are satisfied
    flag=1; i=1;
else
    % The constraints are not satisfied.
    for i=1:maxit
        % Do a Newton step using the least norm solution
        % TODO: Replace this with a calculation that uses A'=Q*R
        q=q-dg(q)'*((dg(q)*dg(q)')\aux); 
        % q=q-dg(q)\aux;
        his(:,i+1)=q;
        % Evaluate the new residual
        aux=g(q); res(i+1)=norm(aux,'inf');
        if (res(i+1)<=tau)
            % The constraints are satisfied
            flag=1; break;
        end
    end
end

% Remove trailing zeros by truncating the arrays
res=res(1:i+1); his=his(:,1:i+1);

% The number of Newton steps
it=i;