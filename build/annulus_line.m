x = load('output.txt');
qvec = 2:length(x)+1;
semilogy(qvec,x, '-o','LineWidth',1.2, 'MarkerSize',10);
set(gca,'fontsize',18);
xlabel('Number of LGL Nodes')
ylabel('Log Error')
ylim([1e-16,1])
title('Error in Area of Half Annulus','Interpreter','latex')
print -depsc2 media/annulus_line
print -dpng media/annulus_line
exit