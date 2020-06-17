%% Function to Unwrap data

function unWrapped = UnWrap(data,scale,correct)
unWrapped = data;
flip = (max(data)-min(data))/scale;

for i=1:length(unWrapped)-1
    if unWrapped(i+1) - unWrapped(i) > flip
        unWrapped(i+1:end) = unWrapped(i+1:end) - flip*scale;
    elseif unWrapped(i+1) - unWrapped(i) < - flip
        unWrapped(i+1:end) = unWrapped(i+1:end) + flip*scale;
    end
end

if correct
    h = figure;
    subplot(2,1,1);
    plot(data);
    set(gca,'FontSize',16);
    title('Original Data');
    axis tight;
    subplot(2,1,2);
    plot(unWrapped);
    set(gca,'FontSize',16);
    title('Unwrapped Data');
    axis tight;
    numBaddies = input('# of pts to change?:');
    if numBaddies > 0
        for bads = 1:numBaddies
            shiftPt = input('Pt to shift?');
            shiftVal = input('Up(+1) or Down(-1)?');
            unWrapped(shiftPt:end) = unWrapped(shiftPt:end) + shiftVal*flip*scale;
            subplot(2,1,2)
            plot(unWrapped);set(gca,'FontSize',16);title('Unwrapped Data');
            axis tight;
        end
        subplot(2,1,2)
        plot(unWrapped);set(gca,'FontSize',16);title('Unwrapped Data');
        axis tight;
    end
    pause;
    delete(h);
end

end