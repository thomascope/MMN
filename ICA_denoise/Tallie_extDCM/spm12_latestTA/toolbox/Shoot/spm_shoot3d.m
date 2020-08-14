function varargout = spm_shoot3d(v0,prm,args, F)
% Geodesic shooting
% FORMAT [phi,Jphi,v1,theta,Jtheta] = spm_shoot3d(v0,prm,args, F)
% v0   - Initial velocity field n1*n2*n3*3 (single prec. float)
% prm  - Differential operator parameters
% prm  - 8 settings
%        - [1][2][3] Voxel sizes
%        - [4][5][6][7][8] Regularisation settings.
%          Regularisation uses the sum of
%          - [4] - absolute displacements
%          - [5] - laplacian
%          - [6] - bending energy
%          - [7] - linear elasticity mu
%          - [8] - linear elasticity lambda
% args - Integration parameters
%        - [1] Num time steps
% F    - optional Fourier transform of Green's function (saves a little time).
%
% phi  - Forward deformation field n1*n2*n3*3 (single prec. float)
% Jphi - Forward Jacobian tensors n1*n2*n3 (single prec. float)
% v1   - Final velocity field n1*n2*n3*3 (single prec. float)
% theta - Inverse deformation field n1*n2*n3*3 (single prec. float)
% Jtheta   - Inverse Jacobian tensors n1*n2*n3 (single prec. float)
%
% This code generates deformations and their Jacobian determinans from
% initial velocity fields by gedesic shooting.  See the work of Miller,
% Younes and others.
%
% LDDMM (Beg et al) uses the following evolution equation:
%     d\phi/dt = v_t(\phi_t)
% where a variational procedure is used to find the stationary solution
% for the time varying velocity field.
% In principle though, once the initial velocity is known, then the
% velocity at subsequent time points can be computed.  This requires
% initial momentum (m_0), computed (using differential operator L) by:
%     m_0 = L v_0
% Then (Ad_{\phi_t})^* m_0 is computed:
%     m_t = |d \phi_t| (d\phi_t)^T m_0(\phi_t)
% The velocity field at this time point is then obtained by using
% multigrid to solve:
%     v_t = L^{-1} m_t
%
% These equations can be found in:
% Younes (2007). "Jacobi fields in groups of diffeomorphisms and
% applications". Quarterly of Applied Mathematics, vol LXV,
% number 1, pages 113-134 (2007).
%
% Note that in practice, (Ad_{\phi_t})^* m_0 is computed differently,
% by multiplying the initial momentum by the inverse of the Jacobian
% matrices of the inverse warp, and pushing the values to their new
% location by the inverse warp (see the "pushg" code of shoot3).
% Multigrid is currently used to obtain v_t = L^{-1} m_t, but
% this could also be done by convolution with the Greens function
% K = L^{-1} (see e.g. Bro-Nielson).
%
%________________________________________________________
% (c) Wellcome Trust Centre for NeuroImaging (2009)

% John Ashburner
% $Id: spm_shoot3d.m 6008 2014-05-22 12:08:01Z john $

spm_diffeo('boundary',0); % Set boundary condition
args0 = [8 4 4];
if nargin<3,
    args = args0;
else
    if numel(args)<numel(args0),
        args = [args args0((numel(args)+1):end)];
    end
end
verb     = false;
N        = args(1);   % # Time steps
d        = size(v0);
d        = d(1:3);
vt       = v0;

if ~isfinite(N),
    % Number of time steps from an educated guess about how far to move
    N = double(floor(sqrt(max(max(max(v0(:,:,:,1).^2+v0(:,:,:,2).^2+v0(:,:,:,3).^2)))))+1);
end

if nargin<4,
    F = spm_shoot_greens('kernel',d,prm);
end

if verb, fprintf('N=%g:', N); end

m0       = spm_diffeo('vel2mom',v0,prm); % Initial momentum (m_0 = L v_0)
if verb, fprintf('\t%g', 0.5*sum(v0(:).*m0(:))/prod(d)); end

% Compute initial small deformation and Jacobian matrices from the velocity.
% The overall deformation and its Jacobian tensor field is generated by
% composing a series of small deformations.
[ phi, Jphi]     = spm_diffeo('smalldef', vt,1/N);

% If required, also compute the forward and possibly also its Jacobian
% tensor field. Note that the order of the compositions will be the
% reverse of those for building the forward deformation.
if nargout>=5,
    [theta,Jtheta] = spm_diffeo('smalldef', vt,-1/N);
elseif nargout>=4,
    theta          = spm_diffeo('smalldef', vt,-1/N);
end

for t=2:abs(N),
    mt             = spm_diffeo('pushg',m0,phi,Jphi);
    %vt            = spm_diffeo('mom2vel', mt, [prm 2 2],vt); % Solve via V-cycles
    vt             = spm_shoot_greens(mt,F,prm);
    if verb, fprintf('\t%g', 0.5*sum(vt(:).*mt(:))/prod(d)); end

    [  dp,   dJ]   = spm_diffeo('smalldef',  vt,1/N);       % Small deformation
    [ phi, Jphi]   = spm_diffeo('comp', dp, phi, dJ, Jphi); % Build up large def. from small defs

    clear dp dJ
%sk=2;plot(phi(:,1:sk:end,1),phi(:,1:sk:end,2),'k-',phi(1:sk:end,:,1)',phi(1:sk:end,:,2)','b-'); axis image; drawnow
%J=squeeze(Jphi);
%dj  = J(:,:,1,1).* J(:,:,2,2) -  J(:,:,1,2).* J(:,:,2,1);
%if any(~isfinite(dj)); crash; end
%imagesc(dj); axis image; drawnow

    % If required, build up forward deformation and its Jacobian tensor field from
    % small deformations
    if nargout>=5,
        [   dp,    dJ] = spm_diffeo('smalldef',  vt,-1/N);        % Small deformation
        [theta,Jtheta] = spm_diffeo('comp', theta,dp, Jtheta,dJ); % Build up large def. from small defs
        clear dp dJ
    elseif nargout>=4,
        dp             = spm_diffeo('smalldef',  vt,-1/N);        % Small deformation
        theta          = spm_diffeo('comp', theta, dp);           % Build up large def. from small defs
        clear dp
    end
    drawnow
end
if verb, fprintf('\n'); end


varargout{1} = phi;
varargout{2} = Jphi;
if nargout>=3,
    mt           = spm_diffeo('pushg',m0,phi,Jphi);
    varargout{3} = spm_shoot_greens(mt,F,prm);
end
if nargout>=4, varargout{4} = theta; end
if nargout>=5, varargout{5} = Jtheta;   end
%__________________________________________________________________________________


