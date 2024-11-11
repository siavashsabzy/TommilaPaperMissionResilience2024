function ParforProgressbar3300(~)
persistent nbitr3300
if isempty(nbitr3300)
    nbitr3300 = 0;
end
nbitr3300 = nbitr3300 + 1;
clc;
fprintf('progress3300: %0.2f \n', nbitr3300 * 100 / 3300);
end