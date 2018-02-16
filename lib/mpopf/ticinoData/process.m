cd /Users/Juraj/Dropbox/USI/PhD/publications/mpbench;

% load 15 min data from Swiss Grid
filename = 'dataTicino960.csv';
M = csvread(filename);

% compute hourly averages
sum = 0;
cnt = 0;
data = [];
for i = 1:length(M)
    sum = sum + M(i);
    cnt = cnt + 1;
    
    if (cnt == 4)
       data = [data sum/4];
       sum = 0;
       cnt = 0;
    end
end

% plot hourly data
plot(data)
dlmwrite('TI240hrs.dat',data');