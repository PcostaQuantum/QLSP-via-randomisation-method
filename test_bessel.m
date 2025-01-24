

Delta = 0.5; 

x = -100:0.05:100; 
y = zeros(size(x)); 

for k = 1:1:max(size(x))
   y(k) = pdf_JLPSS(x(k),Delta);
   y2(k) = pdf_RMopt(x(k),Delta);
end

y = y/y(floor(end/2)-1)*y2(floor(end/2)-1);
y((end+1)/2) = y((end+1)/2-1);
y = y/sum(y)*sum(y2);

plot(x,y,'LineWidth',1.5)
hold on
plot(x,y2,'LineWidth',1.5)
set(gca,'FontSize',18)
axis([-40,40,0,0.3])
legend('JLPSS','Optimal','Location','northeast')
xlabel('t')
ylabel('p(t)')


