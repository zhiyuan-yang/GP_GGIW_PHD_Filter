filter = {'GP','RM'};
%plot GP filter
for i=10:10:30
    figure(i/10);
    hold on;
    for j=1:2
        load(strcat('OSPA_',filter{j},'_',num2str(i)));
        if j == 1
            plot(0.1*(1:300), OSPA_GP,'b')
        else
            plot(0.1*(1:300), OSPA_RM,'r')
        end
    end
    hold off;
end