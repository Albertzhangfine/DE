clc,clear,close all
warning off
format longG
% feature jit off
tic
% ��ֽ����㷨DE Algorithm
% DE ����
maxiter = 200;  % ��������
sizepop = 20;   % ��Ⱥ����
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
F0 = 0.5;          % �Ǳ�����
CR = 0.9;          % �ӽ�����
% ��ʼ����Ⱥ
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
%     fitness(i) = fun([x1,x2]);
    fitness(i) = fun(pop(i,:));
end
% ��¼һ������ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % ȫ�����
fitnesszbest=bestfitness; % ȫ�������Ӧ��ֵ

% ����Ѱ��
for i=1:maxiter
   
    F = F0.^2.*exp(1-maxiter/(maxiter+1-i));
%     F = F0;
%     F = F0*i/maxiter;
    for j=1:sizepop
        
        % ���ѡ��������Ⱥ
        index = randperm(sizepop);
        r1 = index(1);
        r2 = index(2);
        r3 = index(3);
        if j==r1
            r1 = index(4);
        elseif j==r2
            r2 = index(4);
        elseif j==r3
            r3 = index(4);
        end
        
%         newpop = pop(r1,:) + F*(pop(r2,:) - pop(r3,:));
        newpop = zbest + F*(pop(r2,:) - pop(r3,:));
        
        if rand<=CR || isequal(j,r1)
            pop(j,:) = newpop;
        end
        
        % x1  Խ������
        if pop(j,1)>popmax1 || pop(j,1)<popmin1
            pop(j,1) = popmin1 + (popmax1-popmin1)*rand;
        end
        % x2  Խ������
        if pop(j,2)>popmax2 || pop(j,2)<popmin2
            pop(j,2) = popmin2 + (popmax2-popmin2)*rand;
        end
               
        % ��Ӧ�ȸ���
        fitness(j) = fun(pop(j,:));
        
        % �Ƚ�  �����Ƚ�
        if fitness(j)<bestfitness
            bestfitness = fitness(j);
            zbest =  pop(j,:);
        end
        
    end
    fitness_iter(i) = bestfitness;
end
disp('���Ž�')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight

toc