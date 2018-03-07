D = load('output.txt');
x = D(:,1);
y = D(:,2);
qvec = 2:length(x)+1;
semilogy(qvec,x, '-o','LineWidth',1.2, 'MarkerSize',8);
set(gca,'fontsize',18);
xlabel('$q$','interpreter','latex')
ylabel('Log Error')
ylim([1e-16,100])
title('Error in $u_x$','Interpreter','latex')
print -depsc2 media/derivative_x_error
print -dpng media/derivative_x_error

close all
semilogy(qvec,y, '-o','LineWidth',1.2, 'MarkerSize',8);
set(gca,'fontsize',18);
xlabel('$q$','interpreter','latex')
ylabel('Log Error')
ylim([1e-16,100])
title('Error in $u_y$','Interpreter','latex')
print -depsc2 media/derivative_y_error
print -dpng media/derivative_y_error
exit