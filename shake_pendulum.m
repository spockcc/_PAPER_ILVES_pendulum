% Demonstrates different aspects of the SHAKE algorithm

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%  2024-10-03 Initial programming and testing.
%  2024-10-09 Interesting experiment found:
%                tol = 1e-31
%                method = 'quasi'
%                maxit = 200
%                N2 in 16, 32, 64
%             For this experiment we have both decay and growth of the
%             total energy
%             


% Clear work space
clear;

% Test if the experiments have already been conducted

if ~isfile('shake_pendulum.mat')

    % Define the pendulum 
    pendulum1a;

    % Define the time interval
    time1a;

    % Brief description of the first set of experiment
    % Tolerance: tiny
    % Solver: Newton
    % Time step: 
    %   Constant in each experiment. 
    %   Reduced from one experiment to the next
    % Structure: exp1

    % Number of experiments
    num=3;

    % Define tolerance
    tol=1e-12;

    % Choose solver and number of iterations
    method='newton'; maxit=20;

    % Number of recorded states per period
    % It is customary to have at least 5 times steps per period-
    % We have many more ...
    nsp=10;

    % Set the number of samples
    N1=nsp*np;

    % Set the number of time steps between recorded states
    N2=10*2.^(0:num-1);

    % Preallocate space for the results
    r=zeros(dim,N1+1,num);
    v=zeros(dim,N1+1,num);
    t=zeros(1,N1+1,num);
    lambda=zeros(1,N1+1,num);
    it=zeros(1,N1+1,num);
    res=zeros(1,N1+1,num);
    flag=zeros(1,N1+1,num);
    energy=zeros(1,N1+1,num);

    % Conduct the experiments
    for i=1:num
        [r(:,:,i), v(:,:,i), t(:,:,i), lambda(:,:,i),...
            it(:,:,i), res(:,:,i), flag(:,:,i), energy(:,:,i)]=...
            shake(r0, v0, param, a, b, N1, N2(i), method, tol, maxit);
    end
    % Embed the results into a suitable data structure
    exp1.r=r; 
    exp1.v=v;
    exp1.t=t; 
    exp1.lambda=lambda;
    exp1.it=it; 
    exp1.res=res; 
    exp1.flag=flag; 
    exp1.energy=energy;
    exp1.num=num;
    exp1.tol=tol;
    exp1.N1=N1;
    exp1.N2=N2;


    % Tolerance: variable
    % Solver: Simplified Newton method = fixed Jacobian
    % Time step: fixed
    % Structure: exp2
   
    % Define tolerance
    q=2:1:13; tol=10.^(-q);

    % Number of experiments
    num=numel(q);

    % Choose solver and number of iterations
    method='quasi'; maxit=200;

    % Number of recorded states per period
    % It is customary to have at least 5 times steps per period-
    % We have many more ...
    nsp=10;

    % Set the number of samples
    N1=nsp*np;

    % Set the number of time steps between recorded states
    N2=10;

    % Preallocate space for the results
    r=zeros(dim,N1+1,num);
    v=zeros(dim,N1+1,num);
    t=zeros(1,N1+1,num);
    lambda=zeros(1,N1+1,num);
    it=zeros(1,N1+1,num);
    res=zeros(1,N1+1,num);
    flag=zeros(1,N1+1,num);
    energy=zeros(1,N1+1,num);

    % Conduct the experiments
    for i=1:num
        [r(:,:,i), v(:,:,i), t(:,:,i), lambda(:,:,i),...
            it(:,:,i), res(:,:,i), flag(:,:,i), energy(:,:,i)]=...
            shake(r0, v0, param, a, b, N1, N2, method, tol(i), maxit);
    end
    % Embed the results into a suitable data structure
    exp2.r=r; 
    exp2.v=v;
    exp2.t=t; 
    exp2.lambda=lambda;
    exp2.it=it; 
    exp2.res=res; 
    exp2.flag=flag; 
    exp2.energy=energy;
    exp2.num=num;
    exp2.tol=tol;




    % Tolerance: fixed
    % Solver: Newton and Simplified Newton
    % Time step: fixed
    % Structure: exp3
   
    % Define tolerance
    tol=10.^(-4);

    % Number of experiments
    num=2;

    % Choose solver and number of iterations
    method=["newton","quasi"]; maxit=200;

    % Number of recorded states per period
    % It is customary to have at least 5 times steps per period-
    % We have many more ...
    nsp=10;

    % Set the number of samples
    N1=nsp*np;

    % Set the number of time steps between recorded states
    N2=10;

    % Preallocate space for the results
    r=zeros(dim,N1+1,num);
    v=zeros(dim,N1+1,num);
    t=zeros(1,N1+1,num);
    lambda=zeros(1,N1+1,num);
    it=zeros(1,N1+1,num);
    res=zeros(1,N1+1,num);
    flag=zeros(1,N1+1,num);
    energy=zeros(1,N1+1,num);

    % Conduct the experiments
    for i=1:num
        [r(:,:,i), v(:,:,i), t(:,:,i), lambda(:,:,i),...
            it(:,:,i), res(:,:,i), flag(:,:,i), energy(:,:,i)]=...
            shake(r0, v0, param, a, b, N1, N2, method(i), tol, maxit);
    end

    % Embed the results into a suitable data structure
    exp3.r=r; 
    exp3.v=v;
    exp3.t=t; 
    exp3.lambda=lambda;
    exp3.it=it; 
    exp3.res=res; 
    exp3.flag=flag; 
    exp3.energy=energy;
    exp3.num=num;
    exp3.tol=tol;
    exp3.method=method;

    % Save all data
    save('shake_pendulum.mat','exp1','exp2','exp3');
else
    % Recover all data
    load shake_pendulum.mat

    % Generate figure from experiment 1
    f1=figure();
    title(['Tol = ',num2str(exp1.tol)]);
        
    for i=1:exp1.num
        % Extract data
        t=exp1.t(:,:,1); energy=exp1.energy(:,:,i); e0=energy(1,1);
        % Compute drift
        drift=energy-e0; 
        rel=drift/e0;
        
        % Plot the relative drift
        subplot(1,exp1.num,i); plot(t,log2(abs(rel))); grid; grid minor;
        title("N2 = ",num2str(exp1.N2(i)));
        xlabel('time (s)');
        ylabel('log2(abs((E(t) - E(0))/E(0)))');
        xlim([min(t) max(t)]); ylim([-45 -8]);
    end
    % f1.Position=[150 700 560 420];
    % saveas(f1,'good_pendulum.eps');
    print('good_pendulum','-dpdf','-bestfit');

    % Generate figure from experiment 2
    % f2=figure('DefaultAxesPosition', [0.0, 0.0, 1, 1]);
    f2=figure();
    tiledlayout(1,5, 'Padding', 'none', 'TileSpacing', 'compact'); 
    % Clear rel from memory
    clear rel;
    
    for i=1:exp2.num
       
        % Extract data
        t=exp2.t(:,:,1); energy=exp2.energy(:,:,i); e0=energy(1,1);
        % Compute drift
        drift=energy-e0; 
        rel(:,:,i)=drift/e0;
    end
    % Determine common upper and lower bounds for the relative drift
    aux=abs(rel); idx=aux>0; aux=aux(idx);
    p=floor(min(log2(aux))); q=ceil(max(log2(aux)));
    for i=1:5
        nexttile
        % Plot the relative drift
        % subplot(2,5,i); 
        plot(t,log2(rel(:,:,i)));  
        xlabel('time (s)'); 
        ylabel('log2(abs((E(t) - E(0))/E(0)))');
        grid; grid minor; ylim([p q]); 
        title('tol =',num2str(exp2.tol(i)));
        xlim([min(t), max(t)]);
       
    end
    f2=gcf();
    f2.Position=[682.5000  443.5000  2*560.0000  418.0000];
    % saveas(f2,'bad_pendulum.eps');
    print('bad_pendulum','-dpdf','-bestfit');
    % Generate figure from experiment 1
    f3=figure(); 
    
    title(['Tol = ',num2str(exp3.tol)]);
    % Clear rel from memory
    clear rel;
    
    for i=1:exp3.num
        % Extract data
        t=exp3.t(:,:,1); energy=exp3.energy(:,:,i); e0=energy(1,1);
        % Compute drift
        drift=energy-e0; 
        rel(:,:,i)=drift/e0;
    end
   
    for i=1:exp3.num
        % Plot the relative drift
        subplot(1,exp3.num,i); plot(t,sign(rel(:,:,i)),'*');  
        xlabel('time (s)'); 
        ylabel('sign((E(t) - E(0))/E(0))');
        grid; grid minor;
        title('method =',exp3.method(i));
        % legend(exp3.method(i),'Location','East');
        ylim([-1.1,1.1]);
    end
    % f3.Position=[150 100 560 420];
    % saveas(f3,'flip.eps');
    print('flip','-dpdf','-bestfit');
end