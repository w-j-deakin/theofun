function [output] = Rotate(input, angle)

n = length(input(:,1));

for i = 1:n
    output(i,1) = (input(i,1) * cos(angle)) - (input(i,2) * sin(angle));
    output(i,2) = (input(i,1) * sin(angle)) + (input(i,2) * cos(angle));
end

end

