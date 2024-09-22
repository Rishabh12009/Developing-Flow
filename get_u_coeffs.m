function [aP, aE, aW, aN, aS, b] = get_u_coeffs(nx,ny,u,v,p,kdiff,dxc,dxe,dyc,dye)

Fw = zeros(size(u)); Fe = zeros(size(u)); Fn = zeros(size(u)); Fs = zeros(size(u));
Dw = zeros(size(u)); De = zeros(size(u)); Dn = zeros(size(u)); Ds = zeros(size(u));
zo = zeros(size(u)); Sp = zeros(size(u)); Su = zeros(size(u)); b  = zeros(size(u));

% geometric quantities
cellvol = (dxe.*dyc');
celldy = ones(size(dxe)).*dyc';
celldx = dxe.*ones(size(dyc'));

% source terms
b(2:nx,:) = Su(2:nx,:).*cellvol + (p(1:nx-1,:) - p(2:nx,:)).*celldy;


% convective and diffusive fluxes
Fw(2:nx,:) = 0.5*(u(1:nx-1,:     ) + u(2:nx,:)     ) .*celldy;   Dw(2:nx,:) = repmat(kdiff./dxc(1:nx-1) , 1, ny) .*celldy;
Fe(2:nx,:) = 0.5*(u(3:nx+1,:     ) + u(2:nx,:)     ) .*celldy;   De(2:nx,:) = repmat(kdiff./dxc(2:nx) ,   1, ny) .*celldy;
Fs(2:nx,:) = 0.5*(v(1:nx-1,1:ny  ) + v(2:nx,1:ny)  ) .*celldx;   Ds(2:nx,2:ny) = repmat(kdiff./dye', nx-1, 1).*celldx(:,2:ny);
Fn(2:nx,:) = 0.5*(v(1:nx-1,2:ny+1) + v(2:nx,2:ny+1)) .*celldx;   Dn(2:nx,1:ny-1) = repmat(kdiff./dye', nx-1, 1).*celldx(:,1:ny-1);

advec_scheme = 2;
if(advec_scheme==1) % upwind scheme
  aE = De + max(zo, -Fe);  aN = Dn + max(zo, -Fn);  
  aW = Dw + max(Fw,  zo);  aS = Ds + max(Fs,  zo);  
elseif(advec_scheme==2) % C-D scheme
  aE = De - 0.5*Fe;  aN = Dn - 0.5*Fn;  
  aW = Dw + 0.5*Fw;  aS = Ds + 0.5*Fs;  

aS(:, 1) = 0; Fs(:, 1) = 0;   % bottom boundary
aN(:,ny) = 0; Fn(:,ny) = 0;   % top boundary

% calculate aP
aP = aE + aW + aS + aN + (Fe-Fw)+(Fn-Fs);
aP(2:nx,:) = aP(2:nx,:) - Sp(2:nx,:).*cellvol;


end
