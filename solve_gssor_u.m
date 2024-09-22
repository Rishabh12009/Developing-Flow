function u = solve_gssor_u(aP, aE, aW, aN, aS, b, u, nx, ny, alp_relax, max_iter, tol)
  for iter=1:max_iter
    uold = u;
    for i=2:nx
     for j=2:ny-1
       u(i,j) = (1-alp_relax)*uold(i,j) + alp_relax*(b(i,j) + aE(i,j)*u(i+1,j) + aW(i,j)*u(i-1,j) + aN(i,j)*u(i,j+1) + aS(i,j)*u(i,j-1))/aP(i,j);
     end
    end


    % calculate errors and check for convergence
    l2diff = sqrt(sum((u-uold).^2, 'all')/((nx+1)*ny));
    norm_fac = mean(abs(u),'all');
    if(norm_fac < tol) norm_fac = 1; end

    if(l2diff/norm_fac < tol) break; end
  end

end
