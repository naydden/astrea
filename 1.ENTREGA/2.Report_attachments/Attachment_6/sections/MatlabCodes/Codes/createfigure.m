function createfigure(X1, YMatrix1,llegenda1,llegenda2,llegenda3,llegenda4, titol,eixX, eixY)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 30-Nov-2016 11:29:05

% Create figure
figure1 = figure('Name',titol, 'NumberTitle','off');
colormap('hsv');

% Create axes
axes1 = axes('Parent',figure1,'XMinorTick','on','CLim',[0 1]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',2);
set(plot1(1),'Color',[1 0 0],'DisplayName',llegenda1);
set(plot1(2),'Color',[1 0.694117665290833 0.39215686917305],...
    'DisplayName',llegenda2);
set(plot1(3),...
    'Color',[0.0431372560560703 0.517647087574005 0.780392169952393],...
    'DisplayName',llegenda3);

% Create title
title(titol,'Interpreter','latex','FontSize',20);

% Create xlabel
xlabel(eixX,'Interpreter','latex','FontSize',16);

% Create ylabel
ylabel(eixY,'Interpreter','latex','FontSize',16);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','Best');

