function p = solve_gssor_p(aP, aE, aW, aN, aS, b, p, nx, ny, alp_relax, max_iter, tol)

  ppad = zeros(nx+2, ny+2);
  ppad(2:nx+1,2:ny+1) = p;

  for iter=1:max_iter

    ppad(1,:) =ppad(2,:);    ppad(nx+2,:) = ppad(nx+1,:);
    ppad(:,1) = ppad(:,2);    ppad(:,ny+2) = ppad(:,ny+1);

    ppadold = ppad;

    for i=2:nx
     for j=1:ny
       ipad = i+1; jpad=j+1;
       ppad(ipad,jpad) = (1-alp_relax)*ppadold(ipad,jpad) + alp_relax*(b(i,j) + aE(i,j)*ppad(ipad+1,jpad) + aW(i,j)*ppad(ipad-1,jpad) + aN(i,j)*ppad(ipad,jpad+1) + aS(i,j)*ppad(ipad,jpad-1))/(aP(i,j));

     end
    end

    % calculate errors and check for convergence
    l2diff = sqrt(sum((ppad-ppadold).^2, 'all')/((nx+2)*(ny+2)));
    norm_fac = mean(abs(ppad),'all');
    if(norm_fac < tol) 
        norm_fac = 1; 
    end

    if(l2diff/norm_fac < tol) 
        break; 
    end
  end

  p = ppad(2:nx+1,2:ny+1);

end
