function input = symb2input(input_symb)

% convert input matrix from symbolic to function/double

if isnumeric(input_symb)
    input = input_symb;
else
    input = matlabFunction(input_symb);
    if nargin(input) == 0
        input = input();
    end
end