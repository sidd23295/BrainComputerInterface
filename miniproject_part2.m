%% EEG data loading
clc;
clear all;
A = load('Run1.mat');
B = load('Run2.mat');
C = load('Run3.mat');
D = load('Run4.mat');
E = load('Run1_1.mat');
F = load('Run2_1.mat');
G = load('Run3_1.mat');
A = squeeze(A.y);
B = squeeze(B.y);
C = squeeze(C.y);
D = squeeze(D.y);
E = squeeze(E.y);
F = squeeze(F.y);
G = squeeze(G.y);

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

time = B(1,:);
chan_r2 = eegfilt(B(2:17,:),256,8,20);
trigger_r2 = find(diff(B(18,:))==1);
trigger_r2_beg = [3075,trigger_r2];
% figure
% plot(time,chan_r2)

time = C(1,:);
chan_r3 = eegfilt(C(2:17,:),256,8,20);
trigger_r3 = find(diff(C(18,:))==1);
trigger_r3_beg = [3075,trigger_r3];
% figure
% plot(time,chan_r3)

time = D(1,:);
chan_r4 = eegfilt(D(2:17,:),256,8,20);
trigger_r4 = find(diff(D(18,:))==1);
trigger_r4_beg = [3075,trigger_r4];
% figure
% plot(time,chan_r4)

time = E(1,:);
chan_r5 = eegfilt(E(2:17,:),256,8,20);
trigger_r5 = find(diff(E(18,:))==1);
trigger_r5_beg = [3075,trigger_r5];
% figure
% plot(time,chan_r5)

time = F(1,:);
chan_r6 = eegfilt(F(2:17,:),256,8,20);
trigger_r6 = find(diff(F(18,:))==1);
trigger_r6_beg = [3075,trigger_r6];
% figure
% plot(time,chan_r6)

time = G(1,:);
chan_r7 = eegfilt(G(2:17,:),256,8,20);
trigger_r7 = find(diff(G(18,:))==1);
trigger_r7_beg = [3075,trigger_r7];
% figure
% plot(time,chan_r7)


%% total trials

tot_trial = [];
for i = 1:40
    tot_trial = [tot_trial chan_r1(:,trigger_r1_beg(i)+128 ...
        : trigger_r1_beg(i)+1028)];
end

for i = 1:40
    tot_trial = [tot_trial chan_r2(:,trigger_r2_beg(i)+128 ...
        : trigger_r2_beg(i)+1028)];
end

for i = 1:40
    tot_trial = [tot_trial chan_r3(:,trigger_r3_beg(i)+128 ...
        : trigger_r3_beg(i)+1028)];
end

for i = 1:40
    tot_trial = [tot_trial chan_r4(:,trigger_r4_beg(i)+128 ...
        : trigger_r4_beg(i)+1028)];
end

for i = 1:40
    tot_trial = [tot_trial chan_r5(:,trigger_r5_beg(i)+128 ...
        :trigger_r5_beg(i)+1028)];
end

for i = 1:40
    tot_trial = [tot_trial chan_r6(:,trigger_r6_beg(i)+128 ...
        :trigger_r6_beg(i)+1028)];
end

for i = 1:40
    tot_trial = [tot_trial chan_r7(:,trigger_r7_beg(i)+128 ...
        :trigger_r7_beg(i)+1028)];
    
end

%% grouping data
group(:,:,1) = tot_trial(:,1:900);
for i = 1:279
    group(:,:,i+1) = tot_trial(:,(i*900)+1:(i+1)*900);
end

%% Group factors
% key data where right is 1 and left is 0
keydata = [run1.z1(1,:) run2.z2(1,:) run3.z3(1,:) ...
    run4.z4(1,:) run1.z1(1,:) run2.z2(1,:) run3.z3(1,:)]';

%% Trial separation

train_data = zeros(252,16);
test_data = zeros(36,16);

train_data = group(:,:,1:252);
test_data = group(:,:,253:280);

train_std = squeeze(std(train_data,0,2))';
test_std = squeeze(std(test_data,0,2))';

[classifier, error] = classify(test_std,train_std,keydata(1:252));



