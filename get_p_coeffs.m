function [aPp, aE, aW, aN, aS, b, du, dv] = get_p_coeffs(nx,ny,ust,vst,aPu,aPv,dxc,dyc,dt,rho)

  aPp = zeros(nx,ny);
  aE = zeros(nx,ny);  aW = zeros(nx,ny);
  aN = zeros(nx,ny);  aS = zeros(nx,ny);
  b  = zeros(nx,ny);
  du = zeros(size(ust));  dv = zeros(size(vst));

  atmp = repmat(dyc', nx,1);
  aE = atmp.^2./aPu(2:nx+1,:);   b = b - ust(2:nx+1,:).*atmp;
  aW = atmp.^2./aPu(1:nx,:);   b = b + ust(1:nx,:).*atmp;
  du(2:nx,:) = atmp(2:nx,:)./aPu(2:nx,:);

  atmp = repmat(dxc, 1, ny);
  aN = atmp.^2./aPv(:,2:ny+1);   b = b - vst(:,2:ny+1).*atmp;
  aS = atmp.^2./aPv(:,1:ny);   b = b + vst(:,1:ny).*atmp;
  dv(:,2:ny) = atmp(:,2:ny)./aPv(:,2:ny);
   
  aW(1,:) = 0;  aN(:,ny) = 0; aE(nx,:) =0;  aS(:,1) = 0;

  aPp = aE + aW + aN + aS + b;
  
  %outlet b.c.
  aS(nx,:)=0; aN(nx,:)=0; aW(nx,:)=0; aPp(nx,:)=1; b(nx,:)=0;

end
