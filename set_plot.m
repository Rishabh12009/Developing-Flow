function  set_plot(xe,ye,udata,vdata,pdata)

figure(1)
for l= 1:size(udata,3)
h=contourf(xe,ye,udata(2:end,:,l)', 'LineStyle', 'none');
set(gca,'fontsize',14), colorbar,xlabel('Lenght(m)'), ylabel('Diameter(m)'), title(strcat('u-velocity'))
drawnow
end

figure(2)
for l= 1:10:size(pdata,3)
h=contourf(xe,ye,pdata(:,:,l)', 'LineStyle', 'none');
set(gca,'fontsize',14), colorbar,xlabel('Lenght(m)'), ylabel('Diameter(m)'), title(strcat('pressure'))
drawnow
end

end

