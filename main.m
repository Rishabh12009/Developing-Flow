%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Developing Flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Geometry and grid
Lx = 0.1; Ly = 0.01; % Lx:lenght(m) Ly:Diameter(m)
nx = 100; ny = 50; % grid

% physical parameters
rho = 1000; %density
mu = 0.756*10^-3; %viscosity

% input pressure
p1 = 0.1;

%numerical parameters
relax_gssor = 1.2;
relax_p = 0.1; relax_u = 0.7; relax_v = 0.7;
num_max_simple_iter = 2;
tol = 1e-6;
maxiteru = 2; maxiterv = 2; maxiterp = 10000;

[xe, ye, xc, yc, dxe, dye, dxc, dyc] = set_grid(Lx, Ly, nx, ny, 1);


% set initial conditions
u = zeros(nx+1,ny);
v = zeros(nx,ny+1);
p = zeros(nx,ny);

ust = zeros(size(u)); unew = zeros(size(u));
vst = zeros(size(v)); vnew = zeros(size(v));
pst = zeros(size(p)); pnew = zeros(size(p));

% begin iterations
pst = p; p_pr = p*0;
timestep=0;
tend=1;

dt=1e-3;
U = p1*((2*Ly)^2)/(12*mu*Lx);
Re = rho*U*Ly/mu; nu = 1/Re;

for t = 0:dt:tend
    timestep =timestep+1;
for simple_iter = 1:num_max_simple_iter

  % enforce ghost values according to boundary conditions
  [u, v, p] = set_boundary_conditions(nx,ny,u,v,p,p1);

  pst = p; 

  % compute u-coefficients
  [aPu, aE, aW, aN, aS, b] = get_u_coeffs(nx,ny,u,v,pst,nu,dxc,dxe,dyc,dye);
  ust = solve_gssor_u(aPu, aE, aW, aN, aS, b, u, nx, ny, relax_gssor, maxiteru, tol);

  % compute v-coefficients
  [aPv, aE, aW, aN, aS, b] = get_v_coeffs(nx,ny,u,v,pst,nu,dxc,dxe,dyc,dye);
  vst = solve_gssor_v(aPv, aE, aW, aN, aS, b, v, nx, ny, relax_gssor, maxiterv, tol);
  
  % compute p-coefficients
  [aPp, aE, aW, aN, aS, b, du, dv] = get_p_coeffs(nx,ny,ust,vst,aPu,aPv,dxc,dyc,dt,rho);
  p_pr = solve_gssor_p(aPp, aE, aW, aN, aS, b, p_pr, nx, ny, relax_gssor, maxiterp, tol);
  % correct u, v and p
  pnew = pst + (relax_p*p_pr);
  unew(2:nx,:) = (1-relax_u)*u(2:nx,:) + relax_u*(ust(2:nx,:) + du(2:nx,:).*(p_pr(1:nx-1,:)-p_pr(2:nx,:)));
  vnew(:,2:ny) = (1-relax_v)*v(:,2:ny) + relax_v*(vst(:,2:ny) + dv(:,2:ny).*(p_pr(:,1:ny-1)-p_pr(:,2:ny)));

  [unew, vnew, pnew] = set_boundary_conditions(nx,ny,unew,vnew,pnew,p1);

  % check for convergence
  l2u(simple_iter) = sqrt(mean((u - unew).^2,'all'));
  l2v(simple_iter) = sqrt(mean((v - vnew).^2,'all'));
   
  resmax =  max([l2u(simple_iter) l2v(simple_iter)]);

  % prepare for next iteration
  u = unew; v = vnew; p = pnew;
end

 udata(:,:,timestep) = u;
 vdata(:,:,timestep) = v;
 pdata(:,:,timestep) = p;
 u = u + dt.*u;
 v = v + dt.*v;
 p = p + dt.*p;
 disp(strcat('timestep=',num2str(timestep,'%04d;'),' u: ', num2str(l2u(simple_iter),'%04d'), ' v: ', num2str(l2v(simple_iter),'%04d')))
 % if converged, exit the loop
  if(resmax < tol)
      disp('steady state')
      break;
   end
end

set_plot(xe,ye,udata,vdata,pdata);
