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
newpop=[pop;newpop];% 200*10
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
x=[newpop,functionValue];
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
N=10;
% 与其他的进行相比较
for i=1:size(x,1)
    
    individual(i).n=0;
    individual(i).p=[];
    
    for j=1:size(x,1)
        com_equal=0;
        com_more=0;
        com_less=0;
        for k=1:2
            if(x(i,N+k))==((x(j,N+k)))
                com_equal=com_equal+1;
            elseif (x(i,N+k))<((x(j,N+k)))
                com_less=com_less+1;
            else
                com_more=com_more+1;
            end
        end
        % to juge the relationship of two individual
        if com_more==0 &&  ~com_equal==2
            individual(i).p=[individual(i).p j];
        elseif com_less==0 && ~ com_equal==2
            % the number of soultions which dominats it.
            individual(i).n=individual(i).n+1;
        end
        
    end
    if individual(i).n==0
        F(front).f=[F(front).f i];
        x(i,N+2+1)=front+1;
    end
    
end
% 开始去分层
while ~isempty(F(front).f)
    Q=[];
    for i=1:length(F(front).f)
        % 去除掉前一层的数据
        
        if ~isempty(individual(F(front).f(i)).p)
            % to get the length of it
            for j=1:length(individual(F(front.f(i))).p)
                individual(F(front).f(i)).p(j).n=individual(F(front).f(i)).p(j).n-1;
                
                if (individual(F(front.f(i)).p(j)).n==0)
                    Q=[Q  individual(F(front).f(i)).p(j)];
                    % 并且进行相对应, 加入rank数值
                    x(individual(F(front).f(i)).p(j),N+2+1)=front+1;
                    
                end
                
            end
            
        end
        
    end
    % to get next front
    front=front+1;
    F(front).f=Q;
end
%crowding-distance-assignment
% a non-dominated set I
% get the length of

% 要知道rank的具体的信息，从小到达开始排列，得到rank的地址。。。。

% selected ;
[temp,rank]=sort(x(:,N+2+1));

% 得到其具体的数值

% get the whole values
for i=1:length(rank)
    sorted_values(i,:)=x(rank(i),:);
end

% 根据rank 排序
% find the crowding distance for every individual

current_index=0;
for front=1:(length(F)-1)
    
    % to calculate every individual
    % 计算出当前front 的具体的函数数值
    
    % 去得到当前front 的排序
    previous_index=current_index;
    % to get the numbers in this front.....
    for i=1:length(F(front).f)
        
        % 这一层的population 的大小
        y(i,:)= sorted_values(current_index+i,:);
        
        if i==length(F(front).f)
            current_index=current_index+i;
        end
    end
    sorted_based_on_objective=[];
    for k=1:2
        
        [ sorted_based_on_objective, index_of_objective]=sort(y(:,N+k));
        
        % 记录排列后的索引
        
        distance(index_of_objective(1))=Inf;
        distance(length(index_of_objective))=Inf;
        % normalized function value.
        fmin=sorted_based_on_objective(1,N+k);
        fmax=sorted_based_on_objective(length(index_of_objective),N+k);
        % for each objective funtion, the boundary solutions(solutions with the smallest values )
        % to get the other distances
        % 去记录当前索引的距离信息
        y(1,N+2+1+k)=Inf;
        y(length(index_of_objective),N+2+k)=Inf;
        
        for j=2:(length(index_of_objective)-1)
            y(j,N+2+1+k)=(sorted_based_on_objective(j+1,N+k)-sorted_based_on_objective(j-1,N+k))/(fmax-fmin);
        end
        % to calculate the distance
        distance=[];
        % initilize the distance array
        distance(:,1)=zeros(length(F(front).f),1);
        
        % to add the distance to each individual
        for i=1:2
            distance(:,1)=distance(:,1)+ y(:,N+2+1+i);
        end
        % the sum of distance;
        y(:,N+2+2)=distance;
        y=y(:,1:N+2+2);
        
        % to get this
        z(previous_index:current_index,:)=y;
    end
    
    
    
end
% sorting the population according to each objective function value.




end




