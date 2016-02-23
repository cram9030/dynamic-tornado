%global tunnelWing K M C Cdp S_ref C_ref B_ref alpha_0 centroid ActLoc rhom_inf airSpeed

usrfun = 'dynamicTwist2016_2_6';

% setup the snopt input parameters
t = [0:dt:tfin]';
x_init = [6*pi/180*sin(2*pi*5*t)]; %intialize to sine and cos at 5.4 hertz frequency

nnObj   = length(x_init); %number of input variables
%set lower and uper bounds of decision variables to be within physical
%limits
xlow = [-6*pi/180*ones(nnObj,1)];
xupp = [10*pi/180*ones(nnObj,1)]; 

% lower and upper bound of constraints
Flow =[      0; 
             weight-weight*0.01;
             0;
             0;
             (weight-weight*0.1)*ones(nnObj-1,1)]; 
         
Fupp =[      1e6;
             weight+weight*0.01;
             maxAcc;
             0;
             (weight+weight*0.1)*ones(nnObj-1,1)]; 

% Change the following default values if you want better performance
% Read SNOPT.pdf before making any changes to this part.
nnCon  = length(Flow); %number of constraints including the cost function
xmul   = zeros(nnObj,1); xstate = zeros(nnObj,1);
Fmul   = zeros(nnCon,1); Fstate = zeros(nnCon,1);
ObjAdd = 0; ObjRow = 1;
A  = [];  iAfun = [];  jAvar = []; 
[iGfun,jGvar] = find(ones(nnCon,nnObj));

%%%% set options with integer input assignments %%%%%%%%%%%%
%snseti   ( 'Major Iteration limit', 5000); % SNOPT userg guide pp. 62 
snseti   ( 'Derivative option', 0); % '0' some derivatives are unknown, '1': all derivatives are known.
snseti   ( 'Verify level', 3); % verify user provided derivative. '-1': no checking, '1:3' increasing verify accuracy
%snsetr   ( 'Major print level', 3);

% call snopt
solveopt = 1;
Runtime = cputime;
tic

warning('off','all')
[x_opt,F,xmul,Fmul,inform] = snopt( x_init,xlow,xupp,xmul,xstate, ...
				    Flow,Fupp,Fmul,Fstate, ...
				    @dynamicTwist2016_2_6,ObjAdd,ObjRow,A,iAfun,jAvar, ...
				    iGfun,jGvar);
Runtime = cputime-Runtime;
toc
%%%%%%%%%%%%%%%%%%%%%%
F(1)