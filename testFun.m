function [a] = testFun()


for i = 1:20
    a = i;
    if a>11
        return
    end
end