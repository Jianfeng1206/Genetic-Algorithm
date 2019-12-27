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
% 种群的大小为100*10
minvalue=repmat(zeros(1,DNAsize),popnum,1);%个体最小值---B = repmat(A, m, n) %将矩阵A复制m*n块，即B由m*n块A平铺而成
maxvalue=repmat(ones(1,DNAsize),popnum,1); %个体最大值
% to initialize the pop
% 去产生一个种群
% 得到了一个随机的数字的序列
pop=rand(popnum,DNAsize).*rand(popnum,DNAsize).*(maxvalue-minvalue)+minvalue;    %产生新的初始种群
% 开始去迭代计算
for i=1:generations
    % cross over
    % save for child
    newpop=zeros(popnum,DNAsize)
% 采取互相交换DNA的方式，
    for j=1:popnum/2
        % define which position
        
        %怀疑这一段的内容貌似不对
        % 这个函数就是 随机打乱数字序列， 1-popnum 1-100的一个随机数
        k=randperm(popnum);
        %  beta= randn(1,popnum) 产生1-100个的 的正太随机数
        
        %       beta=(randn(1,popnum)) *4
        %      % 觉得这一行代码并没有什么用处
        %      % rand 产生的是0-1 的随机数，有可能为正，有可能是负数的
        %      % 会变成0 和1 的数字。。。。直接用随机数可能更好
        %      m=round(rand(1,popnum)); 方差为1.481
        beta=(-1).^round(rand(1,DNAsize)).*abs(randn(1,DNAsize))*1.481;
        newpop(j*2-1,:)=(pop(k(1),:)+pop(k(2),:))/2+beta.*(pop(k(1),:)-pop(k(2),:))./2;
        %产生第一个子代
        % 这里的是交叉
        newpop(j*2,:)=(pop(k(1),:)+pop(k(2),:))/2-beta.*(pop(k(1),:)-pop(k(2),:))./2;    %产生第二个子代
    end
    % mutation
    % size (A) 返回的数值是m*n
    k=rand(size(newpop));
    miu=rand(size(newpop));% 采用多项式变异
    % 矩阵与一个常数之间的比较。。
    temp=k<1/DNAsize & miu< 0.5 % 要变异的基因位置 生成了一个逻辑矩阵
    % 这个开始返回的是一个逻辑的数，按照多项式去变异
    newpop(temp)=newpop(temp)+(maxvalue(temp)-minvalue(temp)).*((2.*miu(temp)+(1-2.*miu(temp)).*(1-(newpop(temp)-minvalue(temp))./(maxvalue(temp)-minvalue(temp))).^21).^(1/21)-1);        %变异情况一
    newpop(temp)=newpop(temp)+(maxvalue(temp)-minvalue(temp)).*(1-(2.*(1-miu(temp))+2.*(miu(temp)-0.5).*(1-(maxvalue(temp)-newpop(temp))./(maxvalue(temp)-minvalue(temp))).^21).^(1/21));  %变异情况二
end
% 如何去处理越界问题
newpop(newpop>maxvalue)=maxvalue(newpop>maxvalue);
newpop(newpop<minvalue)=minvalue(newpop<minvalue);
%开始去变异
%开始去合并父类和子类, 
newpop=[pop,newpop];% 200*10
%合并父子种群  变成了200*30 的一个矩阵了
%--计算目标函数值 size(newpopulation,1) 返回的是第一列的维度
%计算ZDT1的模型问题
functionValue=zeros(size(newpop,1),2);

% 计算第一维的目标函数 计算出第一个目标数值
functionValue(:,1)=newpop(:,1);
% 去计算第二个目标函数的位置
% sum(A,2) the column vector contains the sum of each row
g=1+9*sum(newpop(:,2:DNAsize),2)./(DNAsize-1);

%计算完毕第二维度的目标数了
functionValue(:,2)=g.*((1-newpop(:,1)./g).^0.5);
% 非支配排序

% 开始去合并两个结果
% there is 22 个列的结果
x=[newpop;functionValue];
%%[functionValue_sorted,newsite]=sortrows(functionValue);

% 去看各个原文的代码
% 去遍历各个种群 200 个个体，相当于200 个solution
% each soultion has two enties. 去计算Sq..
% 这是一个结构体的变量
% initialize the front number to 1
% matlab中的struct的创建，
front=1;
% 这里的这个是相当于结构数组
F(front).f=[]; % used to store the number of the front 
%individual=[];
% 这里的这个相当于结构体数组
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
        %去统计 dominated
        if dom_less==0&& dom_equal~=2
            individual(i).n=individual(i).n+1;
            % 还有一种就是等于0的情况。。。。。。
            % 两种关系，支配与非支配
        elseif dom_more==0&& dom_equal~=2
            % 去记录index    % 数组的自增
            individual(i).p=[individual(i).p j]
        end
    end
    % p belongs to the first front  p 没有被其他的dominate,记录数据
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
        %开始遍历每个front的 dominated set
        if ~(individual(F(front).f(i)).p)
            % how to do  遍历这个集合
            % 表示的是j.. 而不是其他的了
            for j=1:length(F(front).f(i).p)
                % 遍历一次去减去一
                % 表示的是第i个p为集合
                individual(indiviudal(F(front).f(i)).p(j)).n=individual(indiviudal(F(front).f(i)).p(j)).n-1;
                if  individual(indiviudal(F(front).f(i)).p(j)).n==0
                    % 添加到Q之中去
                    % 去记录这个数字的序列号码
                    % 并且将所有满足的条件 rank+1.
                    x(individual(F(front).f(i)).p(j),23)=front+1;
                    Q=[Q,individual(F(front).f(i)).p(j)];
                    
                end
                
                
            end
            
        end
  end
   front=front+1;
   F(front).f=Q;
end
% 这一块到这里去进行一个简单的排序
% 从小到大排列并计算位置

[temp,index_of_fronts]=sort(x(:,23));
% 得到排序的函数。。。。根据 rank去排序
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
    % 其中的这个为输出了
   y(i,:)=sorted_based_on_front(current_index+i);
end
current_index=current_index+i;
%sort each individual based on the objective 
sort_based_on_objective=[];

%感觉合成的变量是比较错误的了
for i=1:2
 % 如何去更好的排序了。。。。。。。。。。。。。。。
 [sort_based_on_objective,index_of_objectives]=sort(y(:,20+i));
   % 重新将这个置为空
 sorted_based_objective=[];
 for j=1:length(index_of_objectives)
  sorted_based_objective(j,:)=y(index_of_objectives(j),:);
 end 
% how to do it ........................
%寻找到最大的数值。
f_max=sorted_based_on_objective(length(index_of_objectives),20+i);
%取得第一行，最后一列的元素
f_min=sorted_based_on_objective(1,20+i);

% 令边界为 infinite 表示为最大的了
y(index_of_objectives(length(index_of_objectives)),20+2+1+i)=Inf;
y(index_of_objectives(1),M + V + 1 + i) = Inf;
for j=2:length(index_of_objectives)-1
 
end
 
 
 
end
    
    
    
end

end










