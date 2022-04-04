function [nX] = comp_normal(X)

dX = [X(2:end,:);X(1,:)]-[X(end,:);X(1:end-1,:)];
vel=sqrt(dX(:,1).^2+dX(:,2).^2);
nX=[dX(:,2) -dX(:,1)]./[vel,vel];