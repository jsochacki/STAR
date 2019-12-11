function create_performance_plot(X1, Y1, Y2, Y3, xmin, xmax, ymin, ymax)

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
ylabel('Measurement Value');
xlabel('Feedback Path SNR (dB)');
title('Transmitter Measurements Vs Feedback Path SNR');
xlim(axes1,[xmin xmax]);
ylim(axes1,[ymin ymax]);
box(axes1,'on');
set(axes1,'XGrid','on','YGrid','on');
legend1 = legend(axes1,'show');
set(legend1,'Location','east');

plot(axes1, X1, Y1,'DisplayName','SNR Of Transmitter Signal (dB)');
plot(axes1, X1, Y2,'DisplayName','EVM Of Transmitter Signal (%)');
plot(axes1, X1, Y3,'DisplayName',['-NMSE Of Difference Between Transmitter' newline ' Signal And Desired Transmit Signal (dB)']);

end