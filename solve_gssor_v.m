function v = solve_gssor_v(aP, aE, aW, aN, aS, b, v, nx, ny, alp_relax, max_iter, tol)
  for iter=1:max_iter
    vold = v;
    for i=2:nx-1
     for j=2:ny
       v(i,j) = (1-alp_relax)*vold(i,j) + alp_relax*(b(i,j) + aE(i,j)*v(i+1,j) + aW(i,j)*v(i-1,j) + aN(i,j)*v(i,j+1) + aS(i,j)*v(i,j-1))/aP(i,j);
     end
    end

    % calculate errors and check for convergence
    l2diff = sqrt(sum((v-vold).^2, 'all')/((ny+1)*nx));
    norm_fac = mean(abs(v),'all');
    if(norm_fac < tol) norm_fac = 1; end

    if(l2diff/norm_fac < tol) break; end
  end

end
