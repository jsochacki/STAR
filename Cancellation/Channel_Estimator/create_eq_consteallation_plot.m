function [figure1 axes1] = create_consteallation_plot(figure_in, axis_in, X1, Y1, X2, Y2, X3, Y3, SNR2, EVM2, SNR3, EVM3, xmin, xmax, ymin, ymax)

if ishandle(figure_in)
   figure1 = figure_in;
   axes1 = axis_in;
else
   figure1 = figure;
   axes1 = axes('Parent',figure1);
   hold(axes1,'on');
   ylabel('Q (V)');
   xlabel('I (V)');
   title('Receive Constellation');
   xlim(axes1,[xmin xmax]);
   ylim(axes1,[ymin ymax]);
   box(axes1,'on');
   set(axes1,'XGrid','on','YGrid','on');
   legend1 = legend(axes1,'show');
   %set(legend1,'Location','north');
   set(legend1,...
      'Position',[0.411317569996558 0.457936510898231 0.206081076223101 0.109730845899562]);
end
if ~isempty(X1)
   plot(axes1, X1,Y1,'DisplayName',num2str([SNR2, EVM2], 'Receive Data Without equalization SNR = %2.2f, EVM = %2.2f'),'Marker','o',...
      'LineStyle','none');
end
if ~isempty(X2)
   plot(axes1, X2,Y2,'DisplayName',num2str([SNR3, EVM3], 'Receive Data With equalization SNR = %2.2f, EVM = %2.2f'),'Marker','o',...
      'LineStyle','none');
end
if ~isempty(X3)
   plot(axes1, X3,Y3,'DisplayName','Noiseless Constellation Points','Marker','o',...
      'LineStyle','none',...
      'Color',[0 0 0]);
end