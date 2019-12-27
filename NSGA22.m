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
newpop=[pop;newpop];% 200*10
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
x=[newpop,functionValue];
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
N=10;
% �������Ľ�����Ƚ�
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
% ��ʼȥ�ֲ�
while ~isempty(F(front).f)
    Q=[];
    for i=1:length(F(front).f)
        % ȥ����ǰһ�������
        
        if ~isempty(individual(F(front).f(i)).p)
            % to get the length of it
            for j=1:length(individual(F(front.f(i))).p)
                individual(F(front).f(i)).p(j).n=individual(F(front).f(i)).p(j).n-1;
                
                if (individual(F(front.f(i)).p(j)).n==0)
                    Q=[Q  individual(F(front).f(i)).p(j)];
                    % ���ҽ������Ӧ, ����rank��ֵ
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

% Ҫ֪��rank�ľ������Ϣ����С���￪ʼ���У��õ�rank�ĵ�ַ��������

% selected ;
[temp,rank]=sort(x(:,N+2+1));

% �õ���������ֵ

% get the whole values
for i=1:length(rank)
    sorted_values(i,:)=x(rank(i),:);
end

% ����rank ����
% find the crowding distance for every individual

current_index=0;
for front=1:(length(F)-1)
    
    % to calculate every individual
    % �������ǰfront �ľ���ĺ�����ֵ
    
    % ȥ�õ���ǰfront ������
    previous_index=current_index;
    % to get the numbers in this front.....
    for i=1:length(F(front).f)
        
        % ��һ���population �Ĵ�С
        y(i,:)= sorted_values(current_index+i,:);
        
        if i==length(F(front).f)
            current_index=current_index+i;
        end
    end
    sorted_based_on_objective=[];
    for k=1:2
        
        [ sorted_based_on_objective, index_of_objective]=sort(y(:,N+k));
        
        % ��¼���к������
        
        distance(index_of_objective(1))=Inf;
        distance(length(index_of_objective))=Inf;
        % normalized function value.
        fmin=sorted_based_on_objective(1,N+k);
        fmax=sorted_based_on_objective(length(index_of_objective),N+k);
        % for each objective funtion, the boundary solutions(solutions with the smallest values )
        % to get the other distances
        % ȥ��¼��ǰ�����ľ�����Ϣ
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




