x = load('output.txt');
qvec = 2:length(x)+1;
semilogy(qvec,x, '-o','LineWidth',1.2, 'MarkerSize',8);
set(gca,'fontsize',18);
xlabel('$q$','interpreter','latex')
ylabel('Log Error')
ylim([1e-16,10])
title('Error in Projection','Interpreter','latex')
print -depsc2 media/coeff_error
print -dpng media/coeff_error
exit