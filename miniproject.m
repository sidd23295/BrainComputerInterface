%% EEG data loading
clc;
clear all;
A = load('Run1.mat');
B = load('Run2.mat');
C = load('Run3.mat');
D = load('Run4.mat');
A = squeeze(A.y);
B = squeeze(B.y);
C = squeeze(C.y);
D = squeeze(D.y);
run1 = load('classrun1.mat');
run2 = load('classrun2.mat');   %good data
run3 = load('classrun3.mat');   %good data
run4 = load('classrun4.mat');

%% filtering of the data

time = A(1,:);
chan_r1 = eegfilt(A(2:17,:),256,8,20);
trigger_r1 = find(diff(A(18,:))==1);
trigger_r1_beg = [3075,trigger_r1];
% figure
% plot(time,chan_r1)
% hold on
% plot(time,chan_r1(7,:))

time = B(1,:);
chan_r2 = eegfilt(B(2:17,:),256,8,20);
trigger_r2 = find(diff(B(18,:))==1);
trigger_r2_beg = [3075,trigger_r2];
% figure
% plot(time,chan_r2)
% hold on
% plot(time,chan_r2(7,:))

time = C(1,:);
chan_r3 = eegfilt(C(2:17,:),256,8,20);
trigger_r3 = find(diff(C(18,:))==1);
trigger_r3_beg = [3075,trigger_r3];
% figure
% plot(time,chan_r3)
% hold on
% plot(time,chan_r3(7,:))

time = D(1,:);
chan_r4 = eegfilt(D(2:17,:),256,8,20);
trigger_r4 = find(diff(D(18,:))==1);
trigger_r4_beg = [3075,trigger_r4];
% figure
% plot(time,chan_r4)
% hold on
% plot(time,chan_r4(7,:))

%% finding the average across all trials for r1

left_eegdatar1 = [];
right_eegdatar1 = [];
for i = 1:40
    if(run1.z1(2,i) == 1)
        left_eegdatar1 = [left_eegdatar1 chan_r1(:,trigger_r1_beg(i)+128:trigger_r1_beg(i)+1028)];
    elseif(run1.z1(1,i) == 1)
        right_eegdatar1 = [right_eegdatar1 chan_r1(:,trigger_r1_beg(i)+128:trigger_r1_beg(i)+1028)];
    else
        continue
    end
end

% left_eegdatar1(8,:) = left_eegdatar1(6,:);
% right_eegdatar1(8,:) = right_eegdatar1(6,:);

%% finding the average across all trials for r2

left_eegdatar2 = [];
right_eegdatar2 = [];
for i = 1:40
    if(run2.z2(2,i) == 1)
        left_eegdatar2 = [left_eegdatar2 chan_r2(:,trigger_r2_beg(i)+128:trigger_r2_beg(i)+1028)];
    elseif(run2.z2(1,i) == 1)
        right_eegdatar2 = [right_eegdatar2 chan_r2(:,trigger_r2_beg(i)+128:trigger_r2_beg(i)+1028)];
    else
        continue
    end
end

% left_eegdatar2(8,:) = left_eegdatar2(6,:);
% right_eegdatar2(8,:) = right_eegdatar2(6,:);

%% finding the average across all trials for r3

left_eegdatar3 = [];
right_eegdatar3 = [];
for i = 1:40
    if(run3.z3(2,i) == 1)
        left_eegdatar3 = [left_eegdatar3 chan_r3(:,trigger_r3_beg(i)+128:trigger_r3_beg(i)+1028)];
    elseif(run3.z3(1,i) == 1)
        right_eegdatar3 = [right_eegdatar3 chan_r3(:,trigger_r3_beg(i)+128:trigger_r3_beg(i)+1028)];
    else
        continue
    end
end

% left_eegdatar3(8,:) = left_eegdatar3(6,:);
% right_eegdatar3(8,:) = right_eegdatar3(6,:);
%% finding the average across all trials for r4

left_eegdatar4 = [];
right_eegdatar4 = [];
for i = 1:40
    if(run2.z2(2,i) == 1)
        left_eegdatar4 = [left_eegdatar4 chan_r4(:,trigger_r4_beg(i)+128:trigger_r4_beg(i)+1028)];
    elseif(run2.z2(1,i) == 1)
        right_eegdatar4 = [right_eegdatar4 chan_r4(:,trigger_r4_beg(i)+128:trigger_r4_beg(i)+1028)];
    else
        continue
    end
end

% left_eegdatar4(8,:) = left_eegdatar4(6,:);
% right_eegdatar4(8,:) = right_eegdatar4(6,:);


%% all of the trials

left_eegdatar1 = left_eegdatar1(:,1:18000);
left_eegdatar2 = left_eegdatar2(:,1:18000);
left_eegdatar3 = left_eegdatar3(:,1:18000);
left_eegdatar4 = left_eegdatar4(:,1:18000);

right_eegdatar1 = right_eegdatar1(:,1:18000);
right_eegdatar2 = right_eegdatar2(:,1:18000);
right_eegdatar3 = right_eegdatar3(:,1:18000);
right_eegdatar4 = right_eegdatar4(:,1:18000);

right_tot = [right_eegdatar1 right_eegdatar2 right_eegdatar3 right_eegdatar4];
left_tot = [left_eegdatar1 left_eegdatar2 left_eegdatar3 left_eegdatar4];

%% avg over left
left_avg = left_tot(:,1:900);
l_trial = zeros(16,900,80);
l_trial(:,:,1) = left_tot(:,1:900);
for i = 1:79
    left_trial = left_tot(:,(i*900)+1:(i+1)*900);
    l_trial(:,:,i+1) = left_trial;
    left_avg = left_avg + left_trial;
end

left_avg = left_avg/80;

%% avg over right
right_avg = right_tot(:,1:900);
r_trial = zeros(16,900,80);
r_trial(:,:,1) = right_tot(:,1:900);

for i = 1:79
    right_trial = right_tot(:,(i*900)+1:(i+1)*900);
    r_trial(:,:,i+1) = right_trial;
    right_avg = right_avg + right_trial;
end

right_avg = right_avg/80;

%% mean about time
ltrialavg = mean(l_trial,2);
rtrialavg = mean(r_trial,2);

figure(1)
scatter(ltrialavg(1,:,:),ltrialavg(2,:,:),'r');
hold on
scatter(rtrialavg(1,:,:),rtrialavg(2,:,:),'b');

%% covariance using each trial

s1 = [];
s2 = [];

for i = 1:80
    s1(:,:,i) = cov(l_trial(:,:,i).');
    s2(:,:,i) = cov(r_trial(:,:,i).');
    [V(:,:,i), G(:,:,i)] = eig(s1(:,:,i),s2(:,:,i));
    W(:,:,i) = V(:,11:16,i);
end

left_CSP_trial = [];
right_CSP_trial = [];

for i = 1:80
    left_CSP_trial(:,:,i) = W(:,:,i)'*l_trial(:,:,i);
end
for i = 1:80
    right_CSP_trial(:,:,i) = W(:,:,i)'*r_trial(:,:,i);
end

%% Covariance using avg trials

S1 = [];
S2 = [];
l_std = std(l_trial,0,2);
r_std = std(r_trial,0,2);

l_avg = mean(l_trial,3);
r_avg = mean(r_trial,3);
S1 = cov(l_avg.');
S2 = cov(r_avg.');
[V, G] = eig(S1,S2);
W = V;

leftCSP = W'*squeeze(l_std);
rightCSP = W'*squeeze(r_std);


%% CSP avg trials - separability

% avg_r_CSP = mean(right_CSP_trial,2);
% avg_l_CSP = mean(left_CSP_trial,2);
% n = 1;
% for i = 1:6
%     for j = 1:6
%         if i ~= j
%             figure(n)
% %             subplot(1,2,1)
%             scatter(avg_r_CSP(i,:,:),avg_r_CSP(j,:,:),'r');
%             hold on
%             scatter(avg_l_CSP(i,:,:),avg_l_CSP(j,:,:),'b');
%             title(['i = ',num2str(i),' j = ',num2str(j)]);
% %             subplot(1,2,2)
%             n = n+1;
%         end
%     end
% end

 %% Channel separability
% n = 1;
% for i = 1:6
%     for j = 1:6
%         if i ~= j
%             figure(n)
%             scatter(left_avg(i,:),left_avg(j,:),'b');
%             hold on
%             scatter(right_avg(i,:),right_avg(j,:),'r');
%             title(['i = ',num2str(i),' j = ',num2str(j)]);
%             n = n+1;
%         end
%     end
% end
% 

%% Separability plots


% avg_r_CSP = std(right_CSP_trial,0,2);
% avg_l_CSP = std(left_CSP_trial,0,2);
n = 1;
for i = 1:6
    for j = 1:6
        if i ~= j
            figure(n)
            subplot(1,2,1)
            scatter(rightCSP(i,:),rightCSP(j,:),'r');
            hold on
            scatter(leftCSP(i,:),leftCSP(j,:),'b');
            title('CSP Space');
            subplot(1,2,2)
            scatter(ltrialavg(i,:,:),ltrialavg(j,:,:),'b');
            hold on
            scatter(rtrialavg(i,:,:),rtrialavg(j,:,:),'r');
            xlabel('Channel 1');
            ylabel('Channel 2');
            legend('Imagined Right', 'Imagined Left');
            title('Channel Space');
            n = n+1;
        end
    end
end

%% Standard deviation data

% std_l_trials = [];
% std_r_trials = [];
% 
% for i = 1:80
%     std_l_trials(:,:,i) = std(left_CSP_trial(:,:,i)');
% end
% 
% for i = 1:80
%     std_r_trials(:,:,i) = std(right_CSP_trial(:,:,i)');
% end
% 
% 
% std_l_trials = squeeze(std_l_trials);
% std_r_trials = squeeze(std_r_trials);
% 
% avg_std_l = mean(std_l_trials,2);
% avg_std_r = mean(std_r_trials,2);

% avg_std_l = std(left_CSP_trial,0,2);
% avg_std_r = std(right_CSP_trial,0,2);
% 
% savg_std_l = squeeze(avg_std_l);
% savg_std_r = squeeze(avg_std_r);
% 
% n = 1;
% for i = 1:6
%     for j = 1:6
%         if i ~= j
%             figure(n)
%             subplot(1,2,1)
%             scatter(savg_std_r(i,:),savg_std_r(j,:),'r');
%             hold on
%             scatter(savg_std_l(i,:),savg_std_l(j,:),'b');
%             title(['i = ',num2str(i),' j = ',num2str(j)]);
%             subplot(1,2,2)
%             scatter(left_avg(i,:),left_avg(j,:),'b');
%             hold on
%             scatter(right_avg(i,:),right_avg(j,:),'r');
%             n = n+1;
%         end
%     end
% end
% 

% scatter(t,avg_l_CSP(1:6),'r');
% hold on
% scatter(t,avg_r_CSP(1:6),'b');

%% topo plots of the Eigen vectors
% W dimension is 16 * 6 * 80 trials
% 
% W_avg = mean(W,3);
% eeglab
% EEG.chanlocs = readlocs('CSP.locs');
% for i = 1:6
%     subplot(3,2,i)
%     topoplot(W_avg(:,i),EEG.chanlocs);
%     title(['Plot of W(',num2str(i),')']);
% end

%% Part 2




%% covariance of avg trials
%
% S1 = cov(left_avg.');
% S2 = cov(right_avg.');
%
% [V, D] = eig(S1,S2);
%
% %% Data using the weights
%
% W = V(:,11:16);
%
%
% left_CSP = transpose(W)*left_avg;
% right_CSP = transpose(W)*right_avg;
%
% %% plot CSP
%
% for i = 1:6
%     figure(i)
%     subplot(1,2,1);
%     plot(left_CSP(i,:))
%     xlim([20 1028])
%     subplot(1,2,2);
%     plot(right_CSP(i,:))
%     xlim([20 1028])
%
% end
%
%
%
%


%%

avg_std_l = std(left_CSP_trial,0,2);

s1 = [];
s2 = [];

for i = 1:80
    s1(:,:,i) = cov(l_trial(:,:,i).');
    s2(:,:,i) = cov(r_trial(:,:,i).');
    [V(:,:,i), G(:,:,i)] = eig(s2(:,:,i),s1(:,:,i));
    W(:,:,i) = V(:,11:16,i);
end

right_CSP_trial = [];

for i = 1:80
    right_CSP_trial(:,:,i) = W(:,:,i)'*r_trial(:,:,i);
end



avg_std_r = std(right_CSP_trial,0,2);

savg_std_l = squeeze(avg_std_l);
savg_std_r = squeeze(avg_std_r);

n = 1;

for i = 1:6
    for j = 1:6
        if i ~= j
            figure(n)
            subplot(1,2,1)
            scatter(savg_std_r(i,:),savg_std_r(j,:),'r');
            hold on
            scatter(savg_std_l(i,:),savg_std_l(j,:),'b');
            xlabel('Feature 1');
            ylabel('Feature 2');
            legend('Imagined Right', 'Imagined Left');
            
%             title(['i = ',num2str(i),' j = ',num2str(j)]);
            subplot(1,2,2)
            scatter(left_avg(i,:),left_avg(j,:),'b');
            hold on
            scatter(right_avg(i,:),right_avg(j,:),'r');
            xlabel('Channel 1');
            ylabel('Channel 2');
            legend('Imagined Right', 'Imagined Left');
            n = n+1;
        end
    end
end

%% Part 2


