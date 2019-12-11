function [figure1 axes1] = create_gain_plot(figure_in, axis_in, X1, Y1, X2, Y2, X3, Y3, xmin, xmax, ymin, ymax)


if ishandle(figure_in)
   figure1 = figure_in;
   axes1 = axis_in;
else
   figure1 = figure;
   axes1 = axes('Parent',figure1);
   hold(axes1,'on');
   xlabel('Input Power (dBm)');
   title('System Gain');
   ylabel('Gain (dB)');
   xlim(axes1,[xmin xmax]);
   ylim(axes1,[ymin ymax]);
   box(axes1,'on');
   set(axes1,'FontName','Times New Roman','XGrid','on','XTick',...
      (xmin:1:xmax),...
      'XTickLabelRotation',45,'YGrid','on');
   legend1 = legend(axes1,'show');
   set(legend1,'Location','northwest','FontSize',10);
end
if ~isempty(X1)
   plot(axes1, X1,Y1,'DisplayName','System Without Pre-Distortion','Marker','.',...
      'LineStyle','none');
end
if ~isempty(X2)
   plot(axes1, X2,Y2,'DisplayName','System With Pre-Distortion','Marker','.',...
      'LineStyle','none');
end
if ~isempty(X3)
   plot(axes1, X3,Y3,'DisplayName','Linear Gain',...
      'Color',[0 0 0]);
end

