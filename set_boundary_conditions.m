function [u, v, p] = set_boundary_conditions(nx,ny,u,v,p,p1)

  % west: rigid wall
  u(1,:) = u(2,:);
  v(1,:) = v(2,:); 
  p(1,:) = p1;

  % east: rigid wall
  u(nx+1,:)= u(nx,:);                   
  v(nx,:) = v(nx-1,:); 
  p(nx,:) = 0.0;
  
  % south: rigid wall
  u(:,1) = 0.0;
  v(:,1) = 0.0;

  % north: rigid wall
  u(:,ny) = 0.0; 
  v(:,ny+1) = 0.0;

end
