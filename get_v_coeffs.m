function [aP, aE, aW, aN, aS, b] = get_v_coeffs(nx,ny,u,v,p,kdiff,dxc,dxe,dyc,dye)

Fw = zeros(size(v)); Fe = zeros(size(v)); Fn = zeros(size(v)); Fs = zeros(size(v));
Dw = zeros(size(v)); De = zeros(size(v)); Dn = zeros(size(v)); Ds = zeros(size(v));
zo = zeros(size(v)); Sp = zeros(size(v)); Sv = zeros(size(v)); b  = zeros(size(v));

% geometric quantities
cellvol = (dxc.*dye');
celldy = ones(size(dxc)).*dye';
celldx = dxc.*ones(size(dye'));

% source terms
b(:,2:ny) = Sv(:,2:ny).*cellvol + (p(:,1:ny-1) - p(:,2:ny)).*celldx;

% convective and diffusive fluxes
Fw(:,2:ny) = 0.5*(u(1:nx,2:ny  ) + u(1:nx  ,1:ny-1)) .*celldy;   Dw(2:nx,2:ny) = repmat(kdiff./dxe , 1, ny-1) .* celldy(2:nx,:);
Fe(:,2:ny) = 0.5*(u(2:nx+1,2:ny) + u(2:nx+1,1:ny-1)) .*celldy;   De(1:nx-1,2:ny) = repmat(kdiff./dxe , 1, ny-1) .* celldy(1:nx-1,:);
Fs(:,2:ny) = 0.5*(v(:,1:ny-1   ) + v(:     ,2:ny  )) .*celldx;   Ds(:,2:ny) = repmat(kdiff./dyc(1:ny-1)', nx, 1) .* celldx;
Fn(:,2:ny) = 0.5*(v(:,3:ny+1   ) + v(:     ,2:ny  )) .*celldx;   Dn(:,2:ny) = repmat(kdiff./dyc(2:ny)'  , nx, 1) .* celldx;


advec_scheme = 2;
if(advec_scheme==1) % upwind scheme
  aE = De + max(zo, -Fe);  aN = Dn + max(zo, -Fn);  
  aW = Dw + max(Fw,  zo);  aS = Ds + max(Fs,  zo);  
elseif(advec_scheme==2) % C-D scheme
  aE = De - 0.5*Fe;  aN = Dn - 0.5*Fn;  
  aW = Dw + 0.5*Fw;  aS = Ds + 0.5*Fs;  

% calculate aP
aP = aE + aW + aS + aN + (Fe-Fw)+(Fn-Fs);
aP(:,2:ny) = aP(:,2:ny) - Sp(:,2:ny).*cellvol;

end
