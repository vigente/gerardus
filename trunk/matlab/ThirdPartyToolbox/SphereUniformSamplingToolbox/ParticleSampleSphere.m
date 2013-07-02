function [V,Tri,Ue_i,Ue]=ParticleSampleSphere(varargin)
% Create an approximately uniform triangular tessellation of the unit 
% sphere by minimizing generalized electrostatic potential energy 
% (aka Reisz s-energy) of the system of charged particles. Effectively, 
% this function produces a locally optimal solution to the problem that 
% involves finding a minimum Reisz s-energy configuration of N equal 
% charges confined to the surface of the unit sphere (s=1 corresponds to 
% the problem originally posed by J. J. Thomson). 
%
% SYNTAX:
% [V,Tri,Ue_i,Ue]=ParticleSampleSphere(option_i,value_i,option_j,value_j,...)
%
% OPTIONS:
%   - 'N'    : desired number of particles. Corresponding value of N must 
%              be a positive interger greater than 9 and less than 1001.
%              N=200 particles is the default setting. 
%   - 'Vo'   : particle positions used to initialize the search.
%              Corresponding value of Vo must be a N-by-3 array, where n is
%              the number of particles. N=10 is the lowest permissible 
%              number of particles. Initializations consisting of more than 
%              1E3 particles are admissible but may lead to unreasonably 
%              long optimization times.
%   - 's'    : Reisz s-energy parameter used to control the strength of
%              particle interactions. Corresponding value of s must be
%              a real number greater than zero. s=1 is the default setting.
%   - 'Etol' : covergence tolerance. Coresponding value of Etol must be
%              a real, positive number. Etol=1E-5 is the default setting.
%   - 'Nitr' : Maximum number of iterations. Corresponding value of Nitr
%              must be a non-negative interger. Nitr=1E3 is the default 
%              setting.
%
% OUTPUT:
%   - V     : N-by-3 array of vertex positions.
%   - Tri   : M-by-3 list of face-vertex connectivities. 
%   - Ue_i  : N-by-1 array of particle energies.
%   - Ue    : K-by-1 array of energy scores, where K-1 was the total number 
%             of iterations. Ue(1) corresponds to the energy of the initial
%             configuration.
%
% EXAMPLE: 
% -------------------------------------------------------------------------
% % Uniformly distribute 100 particles across the surface of the unit sphere
%
%  [V,Tri,~,Ue]=ParticleSampleSphere('N',100);
%  figure, plot(0:numel(Ue)-1,Ue,'.-')
%  set(get(gca,'Title'),'String','Optimization Progress','FontSize',20)
%  xlabel('Iteration #','FontSize',15)
%  ylabel('Reisz s-enrgy','FontSize',15)
%  TR=TriRep(Tri,V);
%  figure, h=trimesh(TR); set(h,'EdgeColor','b'), axis equal
% -------------------------------------------------------------------------
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: June.2012  
%

% Check the inputs
[V,s,Etol,Nitr]=VerifyInputArgs(varargin);
N=size(V,1);    % number of particles
clear varargin
    
% Compute geodesic distances between the points
DOT=V*V';       % dot product 
DOT(DOT<-1)=-1; DOT(DOT>1)=1;
GD=acos(DOT);   % geodesic distance

% Evaluate the energy functional 
GD(1:(N+1):end)=Inf; % set diagonal entries of GD to Inf
Ue_ij=1./((GD.^s)+eps);
Ue_i=sum(Ue_ij,2);
Ue=sum(Ue_i);

% Iteratively optimize the position of the particles along the negative 
% gradient of the energy functional using an adaptive Gauss-Seidel 
% update scheme -----------------------------------------------------------

fprintf('\nWait while particle positions are being optimized ...\n')
t0=clock;
i=0; 

a=ones(N,1);    % step sizes used during position updates
a_min=1E-14;    % minimum step size
a_max=0.1;      % maximum step size
dE=Inf;
while i<Nitr && dE>Etol
    
    i=i+1;
    
    % Sort the particles according to their energy contribution
    [~,idx_sort]=sort(Ue_i,'descend');
    
    % Update the position of individual particles
    for k=1:N
        
        j=idx_sort(k);
        
        idx_j=true(N,1); % particle indices, except the current one 
        idx_j(j)=false;
        
        DOTj=DOT(idx_j,j); 
        GDj=GD(idx_j,j); 
        
        % Compute the gradient for the j-th particle
        dVj=bsxfun(@times,s./(sqrt(1-DOTj.^2)+eps),V(idx_j,:));                
        dVj=bsxfun(@rdivide,dVj,(GDj.^(s+1))+eps);
        if i<5 % remove contributions from particles which are too close
            idx_fail=sum(dVj.^2,2)>(1E4)^2;
            dVj(idx_fail,:)=[];
            if isempty(dVj) % perturb initial positions
                V=V+randn(size(V))/1E5;
                V_L2=sqrt(sum(V.^2,2));
                V=bsxfun(@rdivide,V,V_L2);
                DOT=V*V';
                GD=acos(DOT);
                break
            end 
        end
        dVj=sum(dVj,1);
        
        % Only retain the tangential component of the gradient
        dVj_n=(dVj*V(j,:)')*V(j,:);
        dVj_t=dVj-dVj_n;        
        
        % Adaptively update position of the j-th particle
        m=0; 
        Uj_old=sum(Ue_ij(j,:));
        while true
        
            % Update the position
            Vj_new=V(j,:)-a(j)*dVj_t;
            
            % Constrain the point to surface of the sphere
            Vj_new=Vj_new/norm(Vj_new);
                     
            % Recompute the dot products and the geodesics
            DOTj=sum(bsxfun(@times,V,Vj_new),2);
            DOTj(DOTj<-1)=-1; DOTj(DOTj>1)=1;
            GDj=acos(DOTj);
            GDj(j)=Inf;
            
            Ue_ij_j=1./((GDj.^s)+eps);
            
            % Check if the system potential decreased
            if sum(Ue_ij_j)<=Uj_old
                
                V(j,:)=Vj_new;
                
                DOT(j,:)=DOTj';
                DOT(:,j)=DOTj;
                
                GD(j,:)=GDj';
                GD(:,j)=GDj;
                
                Ue_ij(j,:)=Ue_ij_j';
                Ue_ij(:,j)=Ue_ij_j;
                
                if m==1, a(j)=a(j)*2.5; end
                if a(j)>a_max, a(j)=a_max; end
                
                break
                
            else
                
                if a(j)>a_min
                    a(j)=a(j)/2.5;
                    if a(j)<a_min, a(j)=a_min; end
                else
                    break
                end
                
            end
            
        end
        
    end
    
    % Evaluate the total energy of the system 
    Ue_i=sum(Ue_ij,2);
    Ue(i+1)=sum(Ue_i);

    % Progress update
    if i==1
        fprintf('Iteration#\tEnergy Score\n')
        fprintf('%u\t        %.3f\t\n',0,Ue(1))
    else
        if mod(i,50)==0
        fprintf('%u\t        %.3f\t\n',i,Ue(end))
        end
    end
    
    % Change in energy
    if i>=10, dE=(Ue(i-2)-Ue(i+1))/10; end
 
    % Reset the step sizes
    if mod(i,20)==0, a=a_max*ones(N,1); end
    
end
clear DOT GD Ue_ij
if Nitr==0, fprintf('%u\t        %.3f\t\n',i,Ue(1)); end
if mod(i,50)~=0, fprintf('%u\t        %.3f\t\n',i,Ue(end)); end
fprintf('Optimization completed after %u iterations. Elapsed time: %5.1f sec\n',i,etime(clock,t0))

try
    Tri=fliplr(convhulln(V));
catch err
    msg=sprintf('Unable to triangulate the points. %s',err.message);
    disp(msg)
    Tri=[];
end


%==========================================================================
function [V,s,Etol,Nitr]=VerifyInputArgs(VarsIn)
% Make sure all supplied input argumets have valid format

% Default settings
V=RandSampleSphere(200,'stratified');
s=1; Etol=1E-5; Nitr=1E3;
if isempty(VarsIn), return; end

% First check that there is an even number of inputs
Narg=numel(VarsIn);
if mod(Narg,2)~=0
    error('Incorrect number of input arguments')
end

% Get the properties
FNo={'N','Vo','s','Etol','Nitr'};
flag=false(1,5); exit_flag=false;
for i=1:Narg/2
    
    % Make sure the input is a string
    str=VarsIn{2*(i-1)+1};
    if ~ischar(str)
        error('Input argument #%u is not valid',2*(i-1)+1)
    end
    
    % Get the value
    Val=VarsIn{2*i};
    
    % Match the string against the list of avaliable options  
    chk=strcmpi(str,FNo);
    id=find(chk,1);
    if isempty(id), id=0; end
    
    switch id
        case 1 % number of particles
            
            % Check if 'initialization' option has also been specified
            if flag(2) 
                error('Ambigious combination of options. Specify option ''%s'' or option ''%s'', but not both.',FNo{2},FNo{1})
            end
            
            % Check the format
            if sum(Val(:)<10)>0 || sum(Val(:)>1E3)>0 || ~isreal(Val) || ~isnumeric(Val) || numel(Val)~=1 
                error('Incorrect entry for the ''%s'' option. N must be a positive interger greater than 9 and less than 1001.',FNo{1})
            end
            V=RandSampleSphere(round(Val),'stratified'); %#ok<*NASGU>
            
        case 2 % initialization
        
            % Check if 'number' option has also been specified
            if flag(1) 
                error('Ambigious combination of options. Specify option ''%s'' or option ''%s'', but not both.',FNo{1},FNo{2})
            end
            
            % Check the format
            if ~isreal(Val) || ~isnumeric(Val) || size(Val,2)~=3 || size(Val,1)<10 
                error('Incorrect entry for the ''%s'' option. Vo must be a N-by-3 array, where N is the number of particles. N=10 is the lowest number allowed.',FNo{2})
            end
                        
            % Make sure the particles are constrained to the surface of the unit sphere
            V=Val;
            V_L2=sqrt(sum(V.^2,2));
            V=bsxfun(@rdivide,V,V_L2);
            clear Val V_L2 
            
            % Check if there are more than 1E3 particles
            if size(V,1)>1E3
                
                % Construct a 'yes'/'no' questdlg
                choice = questdlg('Default particle limit exceeded. Would you like to continue?', ...
                    'Particle Limit Exceeded','   YES   ','   NO   ','   NO   ');
                
                % Handle response
                if strcmpi(choice,'   NO   '), exit_flag=true; end
                
            end
            
        case 3 % s parameter
            
            % Check the format
            if Val<0 || ~isreal(Val) || ~isnumeric(Val) || numel(Val)~=1 
                error('Incorrect entry for the ''%s'' option. s must be a positive real number.',FNo{3})
            end
            s=Val;
            
        case 4 % energy tolerance parameter
        
            % Check the format
            if Val<0 || ~isreal(Val) || ~isnumeric(Val) || numel(Val)~=1 
                error('Incorrect entry for the ''%s'' option. Etol must be a positive real number.',FNo{4})
            end
            Etol=Val;
            
        case 5 % maximum number of iterations 
            
            % Check the format
            if Val<0 || ~isreal(Val) || ~isnumeric(Val) || numel(Val)~=1 
                error('Incorrect entry for the ''%s'' option. Nitr must be a non-negative integer.',FNo{5})
            end
            Nitr=Val;
            
        otherwise
            error('''%s'' is not a valid option',str)
    end
    flag(id)=true;
    
end

if exit_flag, Nitr=0; end %#ok<*UNRCH>

