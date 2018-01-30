function display_spectrogram(Y,fs,tw,ov,the_title,MIN_C,MAX_C,step_time,phase)
%%%%%%%%%%%%%%%%% Simple method to display a spectrogram %%%%%%%%%%%%%%%%%%

if(~exist('phase','var'))
    phase=false;
end

if(exist('MIN_C','var'))
    Y(Y<MIN_C)=MIN_C;
    Y(Y>MAX_C)=MAX_C;
    Y(end,end)=MIN_C;
    Y(end,end-1)=MAX_C;
end


if(~exist('step_time','var'))
    step_time = 1;
end

fig=figure;clf(fig);
set(fig,'WindowStyle','docked');
imagesc(Y);

title(the_title,'FontSize',30);

if(phase)
    colormap('HSV');
end

secs = 0:step_time:floor(size(Y,2)*tw*ov/100000);
secs_ticks = secs*100000/(tw*ov);

set(gca,'XTick',secs_ticks);
set(gca,'XTickLabel',secs);
set(get(gca,'XLabel'),'String','Time (s)','FontSize',25)

set(gca,'Xlabel')

freqs = fs/2000:-2:0;
freqs_ticks = size(Y,1)-size(Y,1)/fs*2000*freqs;
freqs_ticks(1)=freqs_ticks(1)+1;

set(gca,'YTick',freqs_ticks);
set(gca,'YTickLabel',freqs); 
set(get(gca,'YLabel'),'String','Frequency (kHz)','FontSize',25)

set(gca,'FontSize',28);

end