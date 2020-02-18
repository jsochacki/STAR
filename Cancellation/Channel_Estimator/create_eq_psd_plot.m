function [figure1 axes1] = create_psd_plot(figure_in, axis_in, X1, Y1, X2, Y2, X3, Y3, oversampling_rate, ymin, ymax)

if ishandle(figure_in)
   figure1 = figure_in;
   axes1 = axis_in;
else
   figure1 = figure;
   axes1 = axes('Parent',figure1);
   hold(axes1,'on');
   xlabel('F/Fs');
   title('Power Spectrum');
   ylabel('Power Spectral Density (dB)');
   xlim(axes1,[-oversampling_rate/2 oversampling_rate/2] ./ oversampling_rate);
   ylim(axes1,[ymin ymax]);
   box(axes1,'on');
   set(axes1,'FontName','Times New Roman','XGrid','on','XTick',...
      (-oversampling_rate/2:0.5:oversampling_rate/2) ./ oversampling_rate,...
      'YGrid','on');
   legend1 = legend(axes1,'show');
   set(legend1,'Location','south','FontSize',10);
end
if ~isempty(X1)
   plot(axes1, X1,Y1,'DisplayName','System Input');
end
if ~isempty(X2)
   plot(axes1, X2,Y2,'DisplayName','System Without equalization');
end
if ~isempty(X3)
   plot(axes1, X3,Y3,'DisplayName','System With equalization');
end
