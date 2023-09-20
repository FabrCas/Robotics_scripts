function plot_it(fs,interval)
    clf
    if not(size(interval) == 2)
        error("intervall not properly set");
    end 
%     disp(length(fs));
    hold on
    for i=1:length(fs)
        fplot(fs(i),interval)
    end
    hold off 
end