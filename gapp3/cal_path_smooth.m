% 计算路径平滑度函数
function [path_smooth] = cal_path_smooth(pop, x)
[n, ~] = size(pop);
path_smooth = zeros(1, n);
for i = 1 : n
    single_pop = pop{i, 1};
    [~, m] = size(single_pop);
    for j = 1 : m - 2
        % 点i所在列（从左到右编号1.2.3...）
        x_now = mod(single_pop(1, j), x) + 1; 
        % 点i所在行（从上到下编号行1.2.3...）
        y_now = fix(single_pop(1, j) / x) + 1;
        % 点i+1所在列、行
        x_next1 = mod(single_pop(1, j + 1), x) + 1;
        y_next1 = fix(single_pop(1, j + 1) / x) + 1;
        % 点i+2所在列、行
        x_next2 = mod(single_pop(1, j + 2), x) + 1;
        y_next2 = fix(single_pop(1, j + 2) / x) + 1;
        %path_smooth(1, i) = path_smooth(1, i) + abs(atan(abs(x_now - x_next1)/abs(y_now - y_next1))-atan(abs(x_next2 - x_next1)/abs(y_next2 - y_next1)));
        %a2 = (x_now - x_next1)^2 + (y_now - y_next1)^2;
        %b2 = (x_next2 - x_next1)^2 + (y_next2 - y_next1)^2;
        c2 = (x_now - x_next2)^2 + (y_now - y_next2)^2;
        %angle = (a2 + c2 - b2) / (2 * sqrt(a2) *  sqrt(c2));
        if c2 < 8 && c2 > 4
            path_smooth(1, i) = path_smooth(1, i) + 5;
        elseif c2 <= 4 && c2 > 1
            path_smooth(1, i) = path_smooth(1, i) + 30;
        elseif    c2 <= 1
            path_smooth(1, i) = path_smooth(1, i) + 5000;
        end
    end
end
