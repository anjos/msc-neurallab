%% Hello emacs, this is -*- matlab -*-
%% $Id: analyze.m,v 1.1 2001/07/13 15:48:02 andre Exp $
%% André Rabello <Andre.Rabello@ufrj.br>

%% This Matlab script reads the data produced during network building
%% and testing and produce some plots for analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN ROUTINE AREA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = readdat(runfile, 15); % Read the run input;

%% Print the 'run' data. Such data is organized as follows:
%% data( 1,:) := current step
%% [TRAIN] data( 2,:) := best SP found for this run
%% [TRAIN] data( 3,:) := threshold value for the best SP found
%% [TRAIN] data( 4,:) := MSE for the current network
%% [TRAIN] data( 5,:) := efficiency for class 1 (usually electrons)
%% [TRAIN] data( 6,:) := efficiency for class 2 (usually jets)
%% [TRAIN] data( 7,:) := PE
%% [TRAIN] data( 8,:) := -1
%% [TEST]  data( 9,:) := best SP found for this run
%% [TEST]  data(10,:) := threshold value for the best SP found
%% [TEST]  data(11,:) := MSE for the current network
%% [TEST]  data(12,:) := efficiency for class 1 (usually electrons)
%% [TEST]  data(13,:) := efficiency for class 2 (usually jets)
%% [TEST]  data(14,:) := PE
%% [TEST]  data(15,:) := 1 in the case this net was saved. 0 otherwise
%%
%% I will plot the following stuff:
%% 1) MSE evolution for train  and test with the saved path;
%% 2) SP evolution for train and test with the saved path;
%% 3) PE evolution for train and test (dont need saved path);

%% MSE Plotting
subplot(1,1,1);
plot(data(1,:), data(4,:), 'g--',data(1,:), data(11,:), 'b-');
title('MSE evolution over training');
xlabel('Training Steps');
ylabel('MSE');
grid on;
legend('train set','test set');
print -depsc2 'mse-evolution.eps';

%% SP Plotting
subplot(1,1,1);
plot(data(1,:), data(2,:), 'g--',data(1,:), data(9,:), 'b-');
title('SP Evolution over training');
xlabel('Training Steps');
ylabel('SP');
grid on;
legend('train set','test set');
print -depsc2 'sp-evolution.eps';

%% PE Plotting
if data(7,1) ~= 0.0,
  subplot(2,1,1);
  plot(data(1,:), data(7,:), 'g-');
  grid;
  title('PE Evolution (TRAIN SET)');
  xlabel('Training Steps');
  ylabel('PE');
  subplot(2,1,2);
  plot(data(1,:), data(14,:), 'b-;test;');
  grid;
  title('PE Evolution (TEST SET)');
  xlabel('Training Steps');
  ylabel('PE');
  print -depsc2 'pe-evolution.eps';
end

%% Evolution off all three above for the test set
subplot(2,1,1);
plot(data(1,:), data(11,:), 'g-');
title('Comparison between MSE and SP quantities');
grid on;
xlabel('Training Steps');
ylabel('MSE');
subplot(2,1,2);
plot(data(1,:), data(9,:), 'b-');
grid on;
xlabel('Training Steps');
ylabel('SP');
print -depsc2 'evolution.eps';

%% The network output for testing and training
subplot(1,1,1);
train = readdat(train_out,2);
test  = readdat(test_out,2);
%% Find the last threshold where I saved data
tra = data(3,data(15,:)>0.5);
tsa = data(10,data(15,:)>0.5);
tra = tra(length(tra));
tsa = tsa(length(tsa));
yline = 0:0.2:1;
trline = tra * ones(size(yline)); clear tra;
tsline = tsa * ones(size(yline)); clear tsa;
grid on;
subplot(2,1,1);
[n1, x1] = myhist(train(1,train(2,:)>0),40,1);
bar(x1,n1, 'g-');
[n2, x2] = myhist(train(1,train(2,:)<0),40,1);
hold on;
bar(x2,n2, 'b-');
plot(trline,yline,'k-');
hold off;
title('Net Output (TRAIN SET)');
xlabel('Outputs');
ylabel('Probability');
legend('threshold','electron','jet');

subplot(2,1,2);
[n1, x1] = myhist(test(1,test(2,:)>0),40,1);
bar(x1,n1,'g-');
[n2, x2] = myhist(test(1,test(2,:)<0),40,1);
hold on;
bar(x2,n2,'b-');
plot(tsline,yline,'k-');
hold off;
title('Net Output (TEST SET)');
xlabel('Outputs');
ylabel('Probability');
legend('threshold','electron','jet');
print -depsc2 'net-output.eps';

%% The Relevance of train and test
subplot(1,1,1);
trrelev = readdat(train_rel, 1);
tsrelev = readdat(test_rel, 1);
subplot(1,2,1);
barh(trrelev,'g-');
title('Train Relevance');
ylabel('Inputs');
xlabel('Relevance');
grid on;
subplot(1,2,2);
barh(tsrelev,'b-');
title('Test Relevance');
ylabel('Inputs');
xlabel('Relevance');
grid on;
print -depsc2 'relevance.eps';


%% The Importance of train and test
trimp = readdat(train_imp, 1);
tsimp = readdat(test_imp, 1);
subplot(1,2,1);
barh(trimp,'g-');
title('Train Importance');
ylabel('Inputs');
xlabel('Importance');
grid on;
subplot(1,2,2);
barh(tsimp, 'b-');
title('Test Importance');
ylabel('Inputs');
xlabel('Importance');
grid on;
print -depsc2 'importance.eps';

%% A comparative test
%subplot(1,1,1);
%barh([tsrelev; tsimp]','g-','stacked');
%grid on;
%legend('Relevance','Importance');
%title('Relevance versus Importance (TEST SET)'); 
%ylabel('Inputs');
%xlabel('Values');
%print -depsc2 'test-relXimp.eps';

%% The output efficiency when compared to (a) given test(s)
subplot(1,1,1);
for i = 1:length(efflist)
  eff = readdat(efflist(i).name, 2);
  if i==1, plot(25*(1-eff(2,:)),eff(1,:), efflist(i).line);
  else
    hold on;
    plot(25*(1-eff(2,:)),eff(1,:), efflist(i).line);
    hold off;
  end
end
grid on;
legend(efflist.label,4); % all at once!
title('Characteristic Curves');
xlabel('Jet background (kHz)');
ylabel('Electron Efficiency');
axis([0 10 0.8 1]); % this is the part I really want to take a look
print -depsc2 'efficiency.eps';






