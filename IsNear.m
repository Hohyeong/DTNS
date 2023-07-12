function tf=IsNear(a,b)
%IsNear a와 b가 거의 같은지 확인

tf = false;

if ~isnumeric(a) || isempty(a) || ~isnumeric(b) || isempty(b)
    warning('Inputs Must be Numeric.');
    tf = false;
    return
else
end

% Tolerance 
tol=sqrt( eps(max(abs(a),abs(b))) );

if any(size(a)~=size(b)) && numel(a)>1 && numel(b)>1
   warning('A and B Must be the Same Size or Either can be a Scalar.')
   tf = false;
   return
end

% main line
tf = abs((a-b)) <= abs(tol);

end

