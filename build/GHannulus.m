x = load('output.txt');
qvec = 2:length(x)+1;
semilogy(qvec,x, '-o','LineWidth',1.2, 'MarkerSize',10);
set(gca,'fontsize',18);
xlabel('Number of LGL Nodes')
ylabel('Log Error')
title('Error in Area of Half Annulus','Interpreter','latex')
print -depsc2 media/GHannulus
print -dpng media/GHannulus
