clc;
clear all;
NSGAII()
% define a function
function NSGAII()
generations=100;
popnum=100;
DNAsize=10;
M=2;
V=10;
% define the number of the objective function 
K=M+V;
% v is the number of decision variables
% ��Ⱥ�Ĵ�СΪ100*10
minvalue=repmat(zeros(1,DNAsize),popnum,1);%������Сֵ---B = repmat(A, m, n) %������A����m*n�飬��B��m*n��Aƽ�̶���
maxvalue=repmat(ones(1,DNAsize),popnum,1); %�������ֵ
% to initialize the pop
% ȥ����һ����Ⱥ
% �õ���һ����������ֵ�����
pop=rand(popnum,DNAsize).*rand(popnum,DNAsize).*(maxvalue-minvalue)+minvalue;    %�����µĳ�ʼ��Ⱥ
% ��ʼȥ��������
for i=1:generations
    % cross over
    % save for child
    newpop=zeros(popnum,DNAsize)
% ��ȡ���ཻ��DNA�ķ�ʽ��
    for j=1:popnum/2
        % define which position
        
        %������һ�ε�����ò�Ʋ���
        % ����������� ��������������У� 1-popnum 1-100��һ�������
        k=randperm(popnum);
        %  beta= randn(1,popnum) ����1-100���� ����̫�����
        
        %       beta=(randn(1,popnum)) *4
        %      % ������һ�д��벢û��ʲô�ô�
        %      % rand ��������0-1 ����������п���Ϊ�����п����Ǹ�����
        %      % ����0 ��1 �����֡�������ֱ������������ܸ���
        %      m=round(rand(1,popnum)); ����Ϊ1.481
        beta=(-1).^round(rand(1,DNAsize)).*abs(randn(1,DNAsize))*1.481;
        newpop(j*2-1,:)=(pop(k(1),:)+pop(k(2),:))/2+beta.*(pop(k(1),:)-pop(k(2),:))./2;
        %������һ���Ӵ�
        % ������ǽ���
        newpop(j*2,:)=(pop(k(1),:)+pop(k(2),:))/2-beta.*(pop(k(1),:)-pop(k(2),:))./2;    %�����ڶ����Ӵ�
    end
    % mutation
    % size (A) ���ص���ֵ��m*n
    k=rand(size(newpop));
    miu=rand(size(newpop));% ���ö���ʽ����
    % ������һ������֮��ıȽϡ���
    temp=k<1/DNAsize & miu< 0.5 % Ҫ����Ļ���λ�� ������һ���߼�����
    % �����ʼ���ص���һ���߼����������ն���ʽȥ����
    newpop(temp)=newpop(temp)+(maxvalue(temp)-minvalue(temp)).*((2.*miu(temp)+(1-2.*miu(temp)).*(1-(newpop(temp)-minvalue(temp))./(maxvalue(temp)-minvalue(temp))).^21).^(1/21)-1);        %�������һ
    newpop(temp)=newpop(temp)+(maxvalue(temp)-minvalue(temp)).*(1-(2.*(1-miu(temp))+2.*(miu(temp)-0.5).*(1-(maxvalue(temp)-newpop(temp))./(maxvalue(temp)-minvalue(temp))).^21).^(1/21));  %���������
end
% ���ȥ����Խ������
newpop(newpop>maxvalue)=maxvalue(newpop>maxvalue);
newpop(newpop<minvalue)=minvalue(newpop<minvalue);
%��ʼȥ����
%��ʼȥ�ϲ����������, 
newpop=[pop,newpop];% 200*10
%�ϲ�������Ⱥ  �����200*30 ��һ��������
%--����Ŀ�꺯��ֵ size(newpopulation,1) ���ص��ǵ�һ�е�ά��
%����ZDT1��ģ������
functionValue=zeros(size(newpop,1),2);

% �����һά��Ŀ�꺯�� �������һ��Ŀ����ֵ
functionValue(:,1)=newpop(:,1);
% ȥ����ڶ���Ŀ�꺯����λ��
% sum(A,2) the column vector contains the sum of each row
g=1+9*sum(newpop(:,2:DNAsize),2)./(DNAsize-1);

%������ϵڶ�ά�ȵ�Ŀ������
functionValue(:,2)=g.*((1-newpop(:,1)./g).^0.5);
% ��֧������

% ��ʼȥ�ϲ��������
% there is 22 ���еĽ��
x=[newpop;functionValue];
%%[functionValue_sorted,newsite]=sortrows(functionValue);

% ȥ������ԭ�ĵĴ���
% ȥ����������Ⱥ 200 �����壬�൱��200 ��solution
% each soultion has two enties. ȥ����Sq..
% ����һ���ṹ��ı���
% initialize the front number to 1
% matlab�е�struct�Ĵ�����
front=1;
% �����������൱�ڽṹ����
F(front).f=[]; % used to store the number of the front 
%individual=[];
% ���������൱�ڽṹ������
% we need to find the q and p
% the number of the variable.......
N=20;


for i=1:size(x,1)
    % the number of solutions which dominate the solution p
    individual(i).n=0;
    individual(i).p=[];
   %a set of solutions that the solution p dominates
    
    % begain the algorithm
    for j=1:size(x,1)
        dom_less=0;
        dom_more=0;
        dom_equal=0;
        % there is two objectives
        for k=1:2
            % judge the objective
            if(x(i,N+k))<x(j,N+k)
                dom_less=dom_less+1;
            elseif x(i,N+k)==x(j,N+k)
                dom_equal=dom_equal+1;
            else
                dom_more=dom_more+1;
            end
        end
        %ȥͳ�� dominated
        if dom_less==0&& dom_equal~=2
            individual(i).n=individual(i).n+1;
            % ����һ�־��ǵ���0�����������������
            % ���ֹ�ϵ��֧�����֧��
        elseif dom_more==0&& dom_equal~=2
            % ȥ��¼index    % ���������
            individual(i).p=[individual(i).p j]
        end
    end
    % p belongs to the first front  p û�б�������dominate,��¼����
    if individual(i).n==0
        x(i,N+2+1)=1;
        % p dominated
        F(front).f=[F(front).f i];
    end
end

% find the subsequent fronts
while ~isempty(F(front).f)
    Q=[];% which is a array to save the next
    % iterating through all of those individuals
    for i=1:length(F(front).f)
        %��ʼ����ÿ��front�� dominated set
        if ~(individual(F(front).f(i)).p)
            % how to do  �����������
            % ��ʾ����j.. ��������������
            for j=1:length(F(front).f(i).p)
                % ����һ��ȥ��ȥһ
                % ��ʾ���ǵ�i��pΪ����
                individual(indiviudal(F(front).f(i)).p(j)).n=individual(indiviudal(F(front).f(i)).p(j)).n-1;
                if  individual(indiviudal(F(front).f(i)).p(j)).n==0
                    % ��ӵ�Q֮��ȥ
                    % ȥ��¼������ֵ����к���
                    % ���ҽ�������������� rank+1.
                    x(individual(F(front).f(i)).p(j),23)=front+1;
                    Q=[Q,individual(F(front).f(i)).p(j)];
                    
                end
                
                
            end
            
        end
  end
   front=front+1;
   F(front).f=Q;
end
% ��һ�鵽����ȥ����һ���򵥵�����
% ��С�������в�����λ��

[temp,index_of_fronts]=sort(x(:,23));
% �õ�����ĺ��������������� rankȥ����
for i=1:length(index_of_fronts)
  sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end
%
current_index=0;
% finding the crowding distance for each individual in each front
for front=1:(length(F)-1)
% all distances is zero
distance=0;
y=[];
previous_index=current_index+1;
for i=1:length(F(front).f)
    % ���е����Ϊ�����
   y(i,:)=sorted_based_on_front(current_index+i);
end
current_index=current_index+i;
%sort each individual based on the objective 
sort_based_on_objective=[];

%�о��ϳɵı����ǱȽϴ������
for i=1:2
 % ���ȥ���õ������ˡ�����������������������������
 [sort_based_on_objective,index_of_objectives]=sort(y(:,20+i));
   % ���½������Ϊ��
 sorted_based_objective=[];
 for j=1:length(index_of_objectives)
  sorted_based_objective(j,:)=y(index_of_objectives(j),:);
 end 
% how to do it ........................
%Ѱ�ҵ�������ֵ��
f_max=sorted_based_on_objective(length(index_of_objectives),20+i);
%ȡ�õ�һ�У����һ�е�Ԫ��
f_min=sorted_based_on_objective(1,20+i);

% ��߽�Ϊ infinite ��ʾΪ������
y(index_of_objectives(length(index_of_objectives)),20+2+1+i)=Inf;
y(index_of_objectives(1),M + V + 1 + i) = Inf;
for j=2:length(index_of_objectives)-1
 
end
 
 
 
end
    
    
    
end

end










