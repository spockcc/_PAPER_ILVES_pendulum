function [q, v, time, lambda, it, residual, flag, energy]=shake(q0, v0, param, a, b, N1, N2, solver, tol, maxit)

% Another implementation of the celebrated SHAKE algorithm
%
% This function attempts to complete N = N1*N2 timesteps. 
% Information is recorded every N2 timesteps.
% 
% INPUT:
%   q0      vector of length dim*n, the initial position of the n particles
%   v0      vector of length dim*n, the initial velocities of the n particles
%   param   a structure describing the simulation
%             param.dim   the dimension of the problem, say, dim = 2 or dim = 3
%             param.m     number of constraints
%             param.n     number of atoms
%             param.g     constraint function
%             param.dg    Jacobian of the constraint function
%             param.phi   the potential 
%             param.f     the forcefield, f = -dphi
%             param.mass  mass(i) is the mass of the ith particle
%   a, b    the simulations covers the time interval (a,b)
%   N1, N2  the total number of time steps is N1*N2
%           the state is saved every N2 time steps  
%   solver  a string that specifies the constraint solver
%              'newton'  Newton's method with the nonsymmetric Jacobian
%              'quasi'   Quasi-Newton method with symmetric approximation
%              
%   tol     tolerance used by Newton's method
%   maxit   maximum number of Newton step  
%
% OUTPUT:
%   q       array such that q(:,i) is the position at time t(i)
%   v       array such that v(:,i) is the velocity at time t(i)
%   t       vector of length at least N1+1 recording the time
%   lambda  array such that lambda(:,i) are the multipliers at time t(i)
%   it      the number of Newton steps needed at time t(i)
%   flag    an integer such that
%             flag=-i   if the solver failed at timestep i
%   energy  vector such that energy(i) is the total energy at time t(i)
%
% MINIMAL WORKING EXAMPLE: shake_mwe1

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2022-02-15 Initial programming and testing
%   2022-08-03 Added support for finite dimensional vector spaces
%   2023-03-06 Improved diagonal solve and updated inline documentation

% /////////////////////////////////////////////////////////////
% Extract relevant parameters and functions from the input
% /////////////////////////////////////////////////////////////

% Isolate the dimension of the vector space
dim=param.dim;

% Isolate the number of constraints
m=param.m;

% Isolate the number of atoms
n=param.n;

% Isolate the constraint function
g=param.g;

% Isolate the Jacobian of the constraint function
G=param.dg;

% Isolate the potential
phi=param.phi;

% Isolate the forcefield
f=param.f;

% Isolate the masses 
mass=param.mass;

% Construct the inverse of the diagonal mass matrix
% M=diag(kron(mass,ones(1,dim)));
iM=diag(kron(1./mass,ones(1,dim)));

% Define the kinetic energy
kin=@(v)kinetic_energy(dim, n, mass,v);

% ////////////////////////////////////////////////////////////
% Allocate space for the output
% ////////////////////////////////////////////////////////////

% Positions
q=zeros(dim*n,N1+1);

% Velocities
v=zeros(dim*n,N1+1);

% Time instances
time=linspace(a,b,N1+1);

% Lagrange multipliers
lambda=zeros(m,N1+1);

% Number of Newton steps
it=zeros(1,N1+1);

% Residual
residual=zeros(1,N1+1);

% Total energy
energy=zeros(1,N1+1);

% ///////////////////////////////////////////////////////////////////
% Compute an initial position that is compatible with the constraints
% ///////////////////////////////////////////////////////////////////

[q0, flag]=initq(q0, param, tol, maxit);

if (flag==0)
    fprintf("Unable to solve initial constraint equation\n");
    return;
end

% At this point we have g(q0) = 0 with the required tolerance.

% ///////////////////////////////////////////////////////////////////
% Record initial values
% ///////////////////////////////////////////////////////////////////

% Initial positions
q(:,1)=q0;

% Initial velocities
v(:,1)=v0;

% ///////////////////////////////////////////////////////////////////
% Temporal matters
% ///////////////////////////////////////////////////////////////////

% The duration of the simulations
T=b-a;

% Basic timestep
h=T/(N1*N2);

% ///////////////////////////////////////////////////////////////////
% Compute v(0.5) and q(1.0)
% ///////////////////////////////////////////////////////////////////

% These are the fundamental equations
%
% v(1/2) = v0 + 0.5*h*inv(M)*(f(q0) - G(q0)'*lambda(0))
%  q(1)  = q0 +     h*v(1/2)
%
% where we must find lambda(0) such that
%
%   g(q(1)) = 0
%
% We need v(1/2) and q(1) in order to contunue

% Assume failure to converge
flag=0;

% These function values will be recycled
f0=f(q0); G0=G(q0);

% Initialize the Lagrange multipliers
l=zeros(m,1);

% Compute initial value of vbar using l=0;
vbar=v0+0.5*h*(iM*f0);

% Compute initial value of qbar
qbar=q0+h*vbar;

% This matrix will be recycled
A=-0.5*h^2*(iM*G0');

% Evaluate constraint function and residual norm
gaux=g(qbar); res=norm(gaux,'inf');

if (res<=tol)
    flag=1; count=0;
else
    for count=1:maxit
        % Do a Newton step
        l=l-(G(qbar)*A)\gaux; 
        % Update vbar
        vbar=v0+0.5*h*(iM*(f0-G0'*l));
        % Update qbar
        qbar=q0+h*vbar;
        % Evaluate constraint function and residual
        gaux=g(qbar); res=norm(gaux,'inf');
        % Test of convergence
        if (res<=tol)
            flag=1;
            break;
        end
    end
end            
    
if (flag==0)
    fprintf('Unable to compute v(0.5)\n'); 
    return;
end
  
% Save number of Newton steps
it(1)=count;

% Save final residual
residual(1)=res;

% Save v(0.5) and q(1);
vold=vbar; qold=qbar;

% Save Lagrange multipliers
lambda(:,1)=l;

% Save initial energy
energy(1)=phi(q0)+kin(v0);

% //////////////////////////////////////////////////////////////////////// 
%   Main loop follows below
% ////////////////////////////////////////////////////////////////////////

% We are using a double loop rather than a single loop to collect the data
% needed for Richardson's techniques.
%
% There are two solvers implemented and while we could have test the choice
% of solver at every timestep, it is faster to have two separate blocks, 
% one for each solver.

if (strcmp(solver,'newton')==1)
    % Do the exact Newton method
    for i=1:N1
        for j=1:N2
            % Compute the one dimensional index k of the current step
            k=(i-1)*N2+j;
            
            % These are the fundamental equations
            %
            %  v(k+0.5) = v(k-0.5) + h*inv(M)*(f(q(k)) - G(q(k))'*lambda(k))
            %  q(k+1  ) = q(k    ) + hv(k+0.5)
            %
            % where we must find lambda(k) such that
            %
            %    g(q(k+1)) = 0
            %
            % We enter the inner loop with vold=v(k-0.5) and qold=v(k)
            % The objective is to compute v(k+0.5) and v(k+1)
            % We will record v(k) and q(k) every N2 steps
            %
            % When k = 1, we have vold=v(0.5) and qold=q(1)
            
            % Initialize Lagrange multiplier
            l=zeros(m,1);
            
            % These function values will be recycled
            G0=G(qold); f0=f(qold);
            
            % This matrix will be recycled
            A=-h^2*(iM*G0');
            
            % Compute initial value of vbar using l=0
            vbar=vold+h*(iM*f0);
            
            % Compute initial value of qbar
            qbar=qold+h*vbar;
            
            % Evaluate constraint function and residual
            gaux=g(qbar); res=max(abs(gaux)); 
            
            % Assume failure to converge!
            flag=0;
            
            if (res<=tol)
                flag=1; count=0;
            else
                for count=1:maxit
                    % Do a Newton step
                    l=l-(G(qbar)*A)\gaux;
                    % Update vbar
                    vbar=vold+h*(iM*(f0-G0'*l));
                    % Update qbar
                    qbar=qold+h*vbar;
                    % Compute residual and residual norm
                    gaux=g(qbar); res=max(abs(gaux));
                    % Test for convergence
                    if (res<=tol)
                        flag=1;
                        break;
                    end
                end
            end
            
            % Test for convergence
            if (flag==0)
                fprintf('Unable to compute Lagrange multipliers at step k = %d\n',k);
                flag=-k;
                return;
            end
            
            
            % ////////////////////////////////////////////////////////////////////////
            % At this point: qbar = q(n+1) and vbar = v(n+1/2)
            % ////////////////////////////////////////////////////////////////////////
            
            % We only need to record values every N2 timesteps
            if (j==N2)
                % Save the values q(n) and v(n)
                q(:,i+1)=qold; v(:,i+1)=(vold+vbar)/2;
                % Save the Langrange multipliers
                lambda(:,i+1)=l;
                % Save the total energy of the system
                energy(i+1)=phi(q(:,i+1))+kin(v(:,i+1));
                % Save number of Newton steps
                it(i+1)=count;
                % Save final residual
                residual(i+1)=res;               
            end
            
            % Prepare for next iteration
            qold=qbar;
            vold=vbar;
        end
    end    
else
    % Use a quasi-Newton solver    
    for i=1:N1
        for j=1:N2
            % Compute the one dimensional index k of the current step
            k=(i-1)*N2+j;
                        
            % These are the fundamental equations
            %
            %  v(k+0.5) = v(k-0.5) + h*inv(M)*(f(q(k)) - G(q(k))'*lambda(k))
            %  q(k+1  ) = q(k    ) + hv(k+0.5)
            %
            % where we must find lambda(k) such that
            %
            %    g(q(k+1)) = 0
            %
            % We enter the inner loop with vold=v(k-0.5) and qold=v(k)
            % The objective is to compute v(k+0.5) and v(k+1)
            % We will record v(k) and q(k) every N2 steps
            %
            % When k = 1, we have vold=v(0.5) and qold=q(1)
            
            % Initialize Lagrange multiplier
            l=zeros(m,1);
            
            % These function values will be recycled
            G0=G(qold); f0=f(qold);
            
            % This matrix will be recycled
            A=-h^2*G0*(iM*G0');
            
            % Compute initial value of vbar using l=0
            vbar=vold+h*(iM*f0);
            
            % Compute initial value of qbar
            qbar=qold+h*vbar;
            
            % Evaluate constraint function and residual
            gaux=g(qbar); res=max(abs(gaux));
                        
            % Assume failure to converge!
            flag=0;
            
            if (res<=tol)
                flag=1; count=0;
            else
                for count=1:maxit
                    % Do a quasi Newton step
                    l=l-A\gaux;
                    % Update vbar
                    vbar=vold+h*(iM*(f0-G0'*l));
                    % Update qbar
                    qbar=qold+h*vbar;
                    % Compute the residual and its norm
                    gaux=g(qbar); res=max(abs(gaux));
                    % Test for convergence
                    if (res<=tol)
                        flag=1;
                        break;
                    end
                end
            end
            
            % Test for convergence
            if (flag==0)
                fprintf('Unable to compute Lagrange multipliers at step k = %d\n',k);
                flag=-k;
                return;
            end
            
            
            % ////////////////////////////////////////////////////////////////////////
            % At this point: qbar = q(k+1) and vbar = v(k+1/2)
            % ////////////////////////////////////////////////////////////////////////
            
            % We only need to record values every N2 timesteps
            if (j==N2)
                % Save the values q(k) and v(k)
                q(:,i+1)=qold; v(:,i+1)=(vold+vbar)/2;
                % Save the Langrange multipliers
                lambda(:,i+1)=l;
                % Save the total energy of the system
                energy(i+1)=phi(q(:,i+1))+kin(v(:,i+1));
                % Save number of Newton steps
                it(i+1)=count;
                % Save final residual
                residual(i+1)=res;
            end
            
            % Prepare for next iteration
            qold=qbar;
            vold=vbar;
        end
    end
end
