function y = expNonlinear(x, alph)
y = 2 ./ (1 + exp(-2*(x + alph * exp(x)))) -1 - alph;
end