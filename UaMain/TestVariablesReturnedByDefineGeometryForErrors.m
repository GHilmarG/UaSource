function TestVariablesReturnedByDefineGeometryForErrors(MUA,s,b,S,B)

if any(isnan(s)) ; save TestSave ; error(' s returned by DefineGeometry contains NaN') ; end
if any(isnan(b)) ; save TestSave ; error(' b returned by DefineGeometry contains NaN') ; end
if any(isnan(S)) ; save TestSave ; error(' S returned by DefineGeometry contains NaN') ; end
if any(isnan(B)) ; save TestSave ; error(' B returned by DefineGeometry contains NaN') ; end

if ~isa(s,'double') ; save TestSave ; error(' s must be a double floating-point variable') ; end
if ~isa(b,'double') ; save TestSave ; error(' b must be a double floating-point variable') ; end
if ~isa(S,'double') ; save TestSave ; error(' S must be a double floating-point variable') ; end
if ~isa(B,'double') ; save TestSave ; error(' B must be a double floating-point variable') ; end




if numel(s)~=MUA.Nnodes ; save TestSave ; error(' Number of elements in s (%-i) differs from number of nodes (%-i) \n',numel(s),MUA.Nnodes) ; end
if numel(b)~=MUA.Nnodes ; save TestSave ; error(' Number of elements in b (%-i) differs from number of nodes (%-i) \n',numel(b),MUA.Nnodes) ; end
if numel(S)~=MUA.Nnodes ; save TestSave ; error(' Number of elements in S (%-i) differs from number of nodes (%-i) \n',numel(S),MUA.Nnodes) ; end
if numel(B)~=MUA.Nnodes ; save TestSave ; error(' Number of elements in B (%-i) differs from number of nodes (%-i) \n',numel(B),MUA.Nnodes) ; end

end