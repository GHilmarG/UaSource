classdef (Sealed) decomposition < matlab.mixin.CustomDisplay & matlab.mixin.internal.Scalar
    %DECOMPOSITION Matrix decomposition
    %   DA = DECOMPOSITION(A) returns a decomposition of matrix A, which
    %   can be used to solve the linear system A*X = B efficiently. The
    %   call X = DA\B returns the same vector as A\B, but is typically
    %   faster.
    %
    %   DA = DECOMPOSITION(A, TYPE) specifies the TYPE of decomposition to
    %   use. Possible values of TYPE depend on the matrix A:
    %
    %   For any matrix A:
    %                   'auto' - (default) Type is based on input matrix.
    %                     'qr' - QR decomposition
    %                    'cod' - COD decomposition
    %
    %   For square matrix A:
    %                     'lu' - LU decomposition
    %                    'ldl' - LDL decomposition (A must be symmetric)
    %                   'chol' - Cholesky decomposition (A must be SPD)
    %             'triangular' - Triangular matrix (A must be triangular)
    %     'permutedTriangular' - Permuted triangular matrix (A must be
    %                               a permutation of a triangular matrix)
    %                 'banded' - Banded LU decomposition
    %             'hessenberg' - LU decomposition of dense Hessenberg matrix
    %               'diagonal' - Diagonal matrix (A must be diagonal)
    %
    %   DECOMPOSITION(A, TYPE, 'upper') and DECOMPOSITION(A, TYPE, 'lower')
    %   specify that only the upper or lower triangular part of A is to be
    %   used. With this syntax, TYPE must be 'ldl', 'chol' or 'triangular'.
    %   For 'ldl' and 'chol', the matrix is assumed to be symmetric; for
    %   'triangular', the matrix is assumed to be triangular.
    %
    %   DA = DECOMPOSITION(..., 'CheckCondition', TF) determines if the
    %   call DA\B throws a warning when A is badly conditioned or of low
    %   rank. Default is true.
    %
    %   DA = DECOMPOSITION(..., 'RankTolerance', TOL) specifies the
    %   tolerance used to determine the rank of matrix A. This option is
    %   used only when TYPE is 'qr' or 'cod', or when TYPE is 'auto' and A
    %   is rectangular; otherwise, this option is ignored.
    %
    %   DA = DECOMPOSITION(..., Name, Value) specifies additional options
    %   for the sparse decomposition functions used, for detailed control
    %   of the decomposition.
    %
    %   decomposition properties:
    %              MatrixSize - Size of the matrix.
    %                    Type - Type of decomposition.
    %   IsConjugateTransposed - Was the object DA (complex conjugate) transposed?
    %             ScaleFactor - Stores scalar multipliers
    %                  IsReal - Is the matrix real?    
    %                IsSparse - Is the matrix sparse?
    %                Datatype - Datatype of the matrix
    %          CheckCondition - Do DA\b and b/DA warn for badly conditioned matrices?
    %
    %   decomposition methods:
    %        isIllConditioned - Determine whether matrix is ill-conditioned.
    %         mldivide (DA\B) - Solve linear system A*X = B.
    %         mrdivide (B/DA) - Solve linear system X*A = B.
    %                    rank - Numerical rank used in mldivide.
    %                   rcond - Reverse condition number of input matrix A.
    %
    %   Decompositions can also be transposed (DA'), negated (-DA), multiplied by
    %   a scalar (c*DA) or divided by a scalar (DA/c).
    
    %   Copyright 2017-2019 The MathWorks, Inc.
    
    properties (GetAccess = public, SetAccess = private)
        %TYPE - Type of the decomposition.
        %   Type is a string describing the type of the decomposition. Type can
        %   be one of:
        %     'lu', 'qr', 'chol', 'ldl', 'triangular', 'permutedTriangular',
        %     'cod', 'banded', 'diagonal', 'hessenberg'
        %
        %   The type is chosen on construction and cannot be modified. On
        %   construction, some of the types can be chosen explicitly.
        %
        %   See also DECOMPOSITION
        Type
        
        %SCALEFACTOR - Scaling factor by which matrix is multiplied
        %   ScaleFactor represents a numeric scalar by which the decomposition
        %   is multiplied. On construction, ScaleFactor is 1. To change
        %   dA.ScaleFactor, use 3*dA or dA/2, for example.
        %
        %   See also DECOMPOSITION
        ScaleFactor = 1
        
        %ISSPARSE - Whether the matrix is sparse
        %   IsSparse indicates whether the original matrix is sparse.
        %
        %   See also DECOMPOSITION
        IsSparse
        
        %ISREAL - Whether the matrix is real
        %   IsReal indicates whether the original matrix is real.
        %
        %   See also DECOMPOSITION
        IsReal
        
        %MATRIXSIZE - Size of the matrix.
        %   MatrixSize is a vector giving the size of the decomposed matrix.
        %   Note that size(DA) does not represent the size of the decomposed
        %   matrix.
        %
        %   See also DECOMPOSITION
        MatrixSize
        
        %ISCONJUGATETRANSPOSED - Transposition mode of the matrix
        %   IsConjugateTransposed indicates whether the decomposition object
        %   has been transposed.
        %
        %   To change dA.IsConjugateTransposed, use dA'.
        %
        %   See also DECOMPOSITION
        IsConjugateTransposed = false
    end
    
    properties (Dependent, GetAccess = public, SetAccess = private)
        %DATATYPE - datatype of the matrix.
        %   Datatype is either 'single' or 'double'.
        Datatype
    end
    
    properties (Dependent, GetAccess = public, SetAccess = public)
        %CHECKCONDITION - Whether a warning is given in \ and /
        %   CheckCondition indicates whether a warning is given in \ and / if
        %   the matrix is ill conditioned. CheckCondition can be set on
        %   construction, or the property can be set in the object.
        %
        %   See also DECOMPOSITION, ISILLCONDITIONED
        CheckCondition
    end
    
    properties (Access = private)
        Underlying
        IsSingle
        useRcond % true unless type is 'qr' or 'cod'
        CheckConditionP = false;
        rcondSaved = []
    end
    
    methods
        function f = decomposition(A, varargin)
            
            if nargin == 0
                A = [];
            elseif ~isfloat(A) || ~ismatrix(A)
                error(message('MATLAB:decomposition:InvalidA'));
            end
            
            [type, uplo, ...
                checkcond, banddensity, mindegree, luPivot, ...
                luSymPivot, ldlPivot, ranktol] = parseInputs(varargin);
            
            switch type
                case 'auto'
                    [f.Underlying, type, f.ScaleFactor] = chooseDecomposition(A, ...
                        banddensity, mindegree, luPivot, luSymPivot, ldlPivot, ranktol);
                case 'qr'
                    if issparse(A)
                        f.Underlying = matlab.internal.decomposition.SparseQR(A, ranktol, mindegree);
                    else
                        f.Underlying = matlab.internal.decomposition.DenseQR(A, ranktol);
                    end
                case 'lu'
                    if size(A, 1) ~= size(A, 2)
                        error(message('MATLAB:decomposition:InvalidAForLU'));
                    end
                    if issparse(A)
                        f.Underlying = matlab.internal.decomposition.SparseLU(A, luPivot, luSymPivot, mindegree);
                    else
                        f.Underlying = matlab.internal.decomposition.DenseLU(A);
                    end
                case 'chol'
                    if size(A, 1) ~= size(A, 2)
                        error(message('MATLAB:decomposition:InvalidAForCholSquare'));
                    end
                    if strcmp(uplo, 'none')
                        if ~ishermitian(A)
                            error(message('MATLAB:decomposition:InvalidAForChol'));
                        end
                    else
                        if strcmp(uplo, 'lower')
                            A = A';
                        end
                    end
                    if issparse(A)
                        f.Underlying = matlab.internal.decomposition.SparseCholesky(A, mindegree);
                    else
                        f.Underlying = matlab.internal.decomposition.DenseCholesky(A);
                    end
                    if ~success(f.Underlying)
                        error(message('MATLAB:decomposition:InvalidAForCholSPD'));
                    end
                case 'ldl'
                    if size(A, 1) ~= size(A, 2)
                        error(message('MATLAB:decomposition:InvalidAForLDLSquare'));
                    end
                    if strcmp(uplo, 'none')
                        if ~ishermitian(A)
                            error(message('MATLAB:decomposition:InvalidAForLDL'));
                        end
                    else
                        if ~isreal(A) && ~isreal(diag(A))
                            error(message('MATLAB:decomposition:InvalidAForLDLDiag'));
                        end
                        if strcmp(uplo, 'upper')
                            A = A';
                        end
                    end
                    if issparse(A)
                        if ~isreal(A)
                            error(message('MATLAB:ldl:complexSparseLDLnotAvailable'));
                        end
                        f.Underlying = matlab.internal.decomposition.SparseLDL(A, ldlPivot, mindegree);
                    else
                        f.Underlying = matlab.internal.decomposition.DenseLDL(A);
                    end
                case 'triangular'
                    if size(A, 1) ~= size(A, 2)
                        error(message('MATLAB:decomposition:InvalidAForTriangSquare'));
                    end
                    
                    if strcmp(uplo, 'none')
                        if istriu(A)
                            uplo = 'upper';
                        elseif istril(A)
                            uplo = 'lower';
                        else
                            error(message('MATLAB:decomposition:InvalidAForTriang'));
                        end
                    end
                    
                    if issparse(A)
                        f.Underlying = matlab.internal.decomposition.SparseTriangular(A, uplo);
                    else
                        f.Underlying = matlab.internal.decomposition.DenseTriangular(A, uplo);
                    end
                    
                case 'permutedTriangular'
                    if size(A, 1) ~= size(A, 2)
                        error(message('MATLAB:decomposition:InvalidAForPermTriang'));
                    end
                    
                    [ist, p, q] = matlab.internal.decomposition.builtin.isPermutedTriangle(A);
                    if ist
                        A = A(p, q);
                        uplo = 'lower';
                    else
                        error(message('MATLAB:decomposition:InvalidAForPermTriang'));
                    end
                    
                    if issparse(A)
                        f.Underlying = matlab.internal.decomposition.SparseTriangular(A, uplo, p, q);
                    else
                        f.Underlying = matlab.internal.decomposition.DenseTriangular(A, uplo, p, q);
                    end
                    
                case 'diagonal'
                    if size(A, 1) ~= size(A, 2) || ~isdiag(A)
                        error(message('MATLAB:decomposition:InvalidAForDiagonal'));
                    end
                    f.Underlying = matlab.internal.decomposition.Diagonal(A);
                case 'hessenberg'
                    if bandwidth(A) > 1 || issparse(A) || size(A, 1) ~= size(A, 2)
                        error(message('MATLAB:decomposition:InvalidAForHessenberg'));
                    end
                    f.Underlying = matlab.internal.decomposition.DenseHessenberg(A);
                case 'banded'
                    
                    if size(A, 1) ~= size(A, 2) || ~isa(A, 'double')
                        error(message('MATLAB:decomposition:InvalidAForBanded'));
                    end
                    
                    % Special case for tridiagonal matrices:
                    [kl, ku] = bandwidth(A);
                    nnzFullBand = size(A,1)*(kl+ku+1)-kl*(kl+1)/2-ku*(ku+1)/2;
                    
                    % Tridiagonal A
                    if kl==1 && ku==1 && isreal(A) && nnz(A) == nnzFullBand
                        f.Underlying = matlab.internal.decomposition.Tridiagonal(sparse(A));
                        if ~success(f.Underlying)
                            f.Underlying = matlab.internal.decomposition.Banded(A, kl, ku);
                        end
                    else
                        f.Underlying = matlab.internal.decomposition.Banded(A, kl, ku);
                    end
                case 'cod'
                    if issparse(A)
                        f.Underlying = matlab.internal.decomposition.SparseCOD(A, ranktol, mindegree);
                    else
                        f.Underlying = matlab.internal.decomposition.DenseCOD(A, ranktol);
                    end
            end
            
            f.Type = type;
            f.MatrixSize = size(A);
            f.IsSparse = issparse(A);
            f.IsReal = isreal(A);
            f.IsSingle = isa(A, 'single');
            
            f.useRcond = ~strcmp(type, 'qr') && ~strcmp(f.Type, 'cod');
            
            if checkcond
                f.CheckCondition = true;
            end
        end
        
        function dt = get.Datatype(f)
            if f.IsSingle
                dt = 'single';
            else
                dt = 'double';
            end
        end
        
        
        %% Addition, as suggested by https://merzba.ch/dw/blg:matlab_decomposition_parfor
        function f = set.Datatype(f, dt)
            assert(any(strcmpi(dt, {'single', 'double'})));
            f.IsSingle = strcmpi(dt, 'single');
        end
        
        
        function doCheck = get.CheckCondition(f)
            doCheck = f.CheckConditionP;
        end
        
        function f = set.CheckCondition(f, doCheck)
            f.CheckConditionP = doCheck;
            if doCheck
                if f.useRcond && isempty(f.rcondSaved)
                    f.rcondSaved = rcond(f.Underlying);
                end
                % Nothing to do in the ~useRcond case; the rank is always
                % computed, no matter what CheckCondition was initially.
            end
        end
        
        function x = mldivide(f,b)
            %\ Matrix left divide for decomposition
            %   X = DA \ B, where DA is a decomposition object, solves the linear
            %   system A*X = B, where DA = DECOMPOSITION(A).
            %
            %   DA2 = S \ DA, where DA is a decomposition object and S is a numeric
            %   scalar, divides DA by the scalar S.
            %
            %   See also: DECOMPOSITION, MRDIVIDE, CTRANSPOSE
            %
            
            if isa(b, 'decomposition')
                if isscalar(f) && isfloat(f)
                    if b.IsConjugateTransposed
                        f = conj(f);
                    end
                    b.ScaleFactor = b.ScaleFactor / double(f);
                    x = b;
                else
                    error(message('MATLAB:decomposition:InvalidDivisor'));
                end
            else % f is a decomposition
                
                if ~isfloat(b) || ~ismatrix(b)
                    error(message('MATLAB:decomposition:InvalidB'));
                end
                
                transp = f.IsConjugateTransposed;
                if transp && strcmp(f.Type, 'qr')
                    error(message('MATLAB:decomposition:QRmldivideTransp'));
                end
                
                castBack = false;
                if f.IsSingle ~= isa(b, 'single')
                    if issparse(b) || f.IsSparse
                        error(message('MATLAB:mldivide:sparseSingleNotSupported'));
                    end
                    if isa(b, 'single')
                        % Note: A\b would cast matrix A to single - but here the decomposition
                        % is already done in double, so return x = single( dA\double(b) ) instead.
                        b = double(b);
                        castBack = true;
                    else
                        b = single(b);
                    end
                end
                
                if issparse(b) && ~f.IsSparse
                    b = full(b);
                end
                
                if (~transp && size(b, 1) ~= f.MatrixSize(1)) || ...
                        (transp && size(b, 1) ~= f.MatrixSize(2))
                    error(message('MATLAB:decomposition:mldivide'));
                end
                
                x = solve(f.Underlying, b, transp);
                
                if f.ScaleFactor ~= 1
                    if transp
                        x = x / conj(f.ScaleFactor);
                    else
                        x = x / f.ScaleFactor;
                    end
                end
                
                if castBack
                    x = single(x);
                end
                
                if f.CheckConditionP
                    warnIfIllConditioned(f);
                end
            end
        end
        
        function x = mrdivide(b,f)
            %/ Matrix right divide for decomposition
            %   X = B / DA, where DA is a decomposition object, solves the linear
            %   system X*A = B, where DA = DECOMPOSITION(A).
            %
            %   DA2 = DA / S, where DA is a decomposition object and S is a numeric
            %   scalar, divides DA by the scalar S.
            %
            %   See also: DECOMPOSITION, MLDIVIDE, CTRANSPOSE
            %
            
            if isa(b, 'decomposition')
                if isscalar(f) && isfloat(f)
                    if b.IsConjugateTransposed
                        f = conj(f);
                    end
                    b.ScaleFactor = b.ScaleFactor / double(f);
                    x = b;
                else
                    error(message('MATLAB:decomposition:InvalidDivisor'));
                end
            else % f is a decomposition
                
                if ~isfloat(b) || ~ismatrix(b)
                    error(message('MATLAB:decomposition:InvalidB'));
                end
                
                transp = f.IsConjugateTransposed;
                if ~transp && strcmp(f.Type, 'qr')
                    error(message('MATLAB:decomposition:QRmrdivideTransp'));
                end
                
                castBack = false;
                if f.IsSingle ~= isa(b, 'single')
                    if f.IsSparse || issparse(b)
                        error(message('MATLAB:mrdivide:sparseSingleNotSupported'));
                    end
                    if isa(b, 'single')
                        % Note: A\b would cast matrix A to single - but here the decomposition
                        % is already done in double, so return x = single( dA\double(b) ) instead.
                        b = double(b);
                        castBack = true;
                    else
                        b = single(b);
                    end
                end
                
                if issparse(b) && ~f.IsSparse
                    b = full(b);
                end
                
                if (~transp && size(b, 2) ~= f.MatrixSize(2)) || ...
                        (transp && size(b, 2) ~= f.MatrixSize(1))
                    error(message('MATLAB:decomposition:mrdivide'));
                end
                
                x = solve(f.Underlying, b', ~transp)';
                
                if f.ScaleFactor ~= 1
                    if transp
                        x = x / conj(f.ScaleFactor);
                    else
                        x = x / f.ScaleFactor;
                    end
                end
                
                if castBack
                    x = single(x);
                end
                
                if f.CheckConditionP
                    warnIfIllConditioned(f);
                end
            end
        end
        
        function f = mtimes(b, f)
            %* Scalar multiply for decomposition
            %    B = DA * S or B = S * DA, where DA is a decomposition and S is a
            %    numeric scalar, returns the scaled decomposition B.
            %
            %    See also TIMES
            %
            if isobject(b)
                tmp = b;
                b = f;
                f = tmp;
            end
            if ~isscalar(b) || ~isfloat(b)
                error(message('MATLAB:decomposition:InvalidMult'));
            end
            if f.IsConjugateTransposed
                b = conj(b);
            end
            f.ScaleFactor = double(b)*f.ScaleFactor;
        end
        
        function f = times(b, f)
            %.* Scalar multiply for decomposition
            %    B = DA .* S or B = S .* DA, where DA is a decomposition and S is a
            %    numeric scalar, returns the scaled decomposition B.
            %
            %    See also MTIMES
            %
            if isobject(b)
                tmp = b;
                b = f;
                f = tmp;
            end
            if ~isscalar(b) || ~isfloat(b)
                error(message('MATLAB:decomposition:InvalidMult'));
            end
            if f.IsConjugateTransposed
                b = conj(b);
            end
            f.ScaleFactor = double(b)*f.ScaleFactor;
        end
        
        function f = rdivide(f, b)
            %./ Scalar right divide for decomposition
            %   DA2 = DA ./ S, where DA is a decomposition object and S is a
            %   numeric scalar, divides DA by the scalar S.
            %
            %   See also: DECOMPOSITION, MLDIVIDE, CTRANSPOSE
            %
            if ~isscalar(b) || ~isfloat(b)
                error(message('MATLAB:decomposition:InvalidDivisor'));
            end
            if f.IsConjugateTransposed
                b = conj(b);
            end
            f.ScaleFactor = f.ScaleFactor/double(b);
        end
        
        function f = ldivide(b, f)
            %.\ Scalar right divide for decomposition
            %   DA2 = S .\ DA, where DA is a decomposition object and S is a
            %   numeric scalar, divides DA by the scalar S.
            %
            %   See also: DECOMPOSITION, MLDIVIDE, CTRANSPOSE
            %
            if ~isscalar(b) || ~isfloat(b)
                error(message('MATLAB:decomposition:InvalidDivisor'));
            end
            if f.IsConjugateTransposed
                b = conj(b);
            end
            f.ScaleFactor = f.ScaleFactor/double(b);
        end
        
        function f = ctranspose(f)
            %' Complex conjugate transpose
            %    DA' is the conjugate transpose of decomposition DA. CTRANSPOSE(DA)
            %    is equivalent to DA'.
            %
            %    Example:
            %    X = DA'\B solves the linear system A'*X = B.
            %
            %    See also DECOMPOSITION, MLDIVIDE
            %
            f.IsConjugateTransposed = ~f.IsConjugateTransposed;
        end
        
        function f = uminus(f)
            %- Uminus Unary minus of decomposition
            %    -DA is a decomposition scaled by -1.
            %
            %    X = -DA \ B is equal to X = -(DA \ B).
            %
            %    See also DECOMPOSITION, MLDIVIDE
            %
            f.ScaleFactor = -f.ScaleFactor;
        end
        
        function isIllCond = isIllConditioned(f)
            %ISILLCONDITIONED Checks whether matrix is ill-conditioned
            %    TF = ISILLCONDITIONED(DA) returns true if the original matrix A is
            %    ill-conditioned.
            %
            %    If ISILLCONDITIONED(DA) is true, then DA\B and B/DA display a
            %    warning. This warning can be turned off with the 'CheckCondition'
            %    property.
            %
            %    A matrix is diagnosed as well-conditioned if
            %      - Type is 'qr' or 'cod', and rank(DA) is equal to min(size(A)).
            %      - Type is not 'qr' or 'cod', and rcond(DA) is greater than eps.
            %
            %    See also RCOND, RANK
            %
            if f.useRcond
                rc = rcond(f);
                isIllCond = rc < eps(class(rc)) || isnan(rc);
            else
                isIllCond = f.Underlying.rank_ < min(f.MatrixSize);
            end
        end
        
        function rc = rcond(f)
            %RCOND Estimated reverse condition number
            %    RC = RCOND(DA) returns an estimate of the reverse condition number
            %    of the original matrix A. This is the same estimate that A\B uses
            %    to warn if matrix A is ill-conditioned.
            %
            %    If DA.Type is 'qr' or 'cod', then RCOND(DA) is not supported.
            %
            %    See also ISILLCONDITIONED, RANK
            %
            if ~f.useRcond
                error(message('MATLAB:decomposition:RcondNotSupported'));
            elseif ~isempty(f.rcondSaved)
                rc = f.rcondSaved;
            else
                rc = rcond(f.Underlying);
            end
        end
        
        function rk = rank(f)
            %RANK Estimated rank of decomposition object
            %    RC = RANK(DA) returns an estimate of the rank of the original
            %    matrix A.
            %
            %    If DA.Type is 'qr' or 'cod', then this is the same estimate that
            %    A\B uses to warn if matrix A is ill-conditioned. This estimate is
            %    based on a QR decomposition of A, and therefore not equivalent to
            %    RANK(A), which uses the SVD of A.
            %
            %    If DA.Type is not 'qr' or 'cod', then RANK is not supported.
            %
            %    See also ISILLCONDITIONED, RCOND.
            %
            if f.useRcond
                error(message('MATLAB:decomposition:RankNotSupported'));
            else
                rk = f.Underlying.rank_;
            end
        end
        
    end
    
    
    methods (Hidden)
        function f = uplus(f)
            % Do nothing
        end
        
        function f = transpose(~) %#ok<STOUT>
            %.' Transpose
            %    Not supported for decomposition.
            error(message('MATLAB:decomposition:TransposeNotSupported'));
        end
        
        %% Changes as suggested by https://merzba.ch/dw/blg:matlab_decomposition_parfor
        %function s = saveobj(~)
        %   s = struct;
        %   warning(message('MATLAB:decomposition:SaveNotSupported'));
        % end
         function s = saveobj(obj,~)
   
            mc = metaclass(obj);
            props = {mc.PropertyList.Name}';
            s = struct;
            for ii = 1 : numel(props)
                s.(props{ii}) = obj.(props{ii});
            end
        end
        %%
    end
    
    methods (Hidden, Static)
        %% Changes as suggested by https://merzba.ch/dw/blg:matlab_decomposition_parfor
        % function obj = loadobj(~)
        %   warning(message('MATLAB:decomposition:LoadNotSupported'));
        %   obj = decomposition;
        %  end
        function obj = loadobj(s)
            fns = fieldnames(s);
            obj = decomposition([]);
            for ii = 1 : numel(fns)
                obj.(fns{ii}) = s.(fns{ii});
            end
        end
        
    end
    
    methods (Access = protected)
        
        function str = getPropertyGroups(obj)
            
            propertyList = {'MatrixSize', 'Type'};
            
            if obj.IsConjugateTransposed
                propertyList{end+1} = 'IsConjugateTransposed';
            end
            if obj.ScaleFactor ~= 1
                propertyList{end+1} = 'ScaleFactor';
            end
            
            str = matlab.mixin.util.PropertyGroup(propertyList);
        end
        
        function str = getAllPropertyGroups(~)
            propertyList = {'MatrixSize', 'Type', 'IsConjugateTransposed', ...
                'ScaleFactor', 'IsReal', 'IsSparse', ...
                'Datatype', 'CheckCondition'};
            
            str = matlab.mixin.util.PropertyGroup(propertyList);
        end
        
        function footer = getFooter(~)
            name = inputname(1);
            linktext = getString(message('MATLAB:graphicsDisplayText:AllPropertiesText'));
            msg = getString(message('MATLAB:graphicsDisplayText:FooterLinkFailureMissingVariable',name));
            cmd = sprintf('if exist(''%s'',''var''),displayAllProperties(%s),else,disp(''%s'');end',name,name,msg);
            footer = sprintf(['  Show <a href="matlab:%s">' ...
                '%s</a>\n'],cmd,linktext);
        end
        
    end
    
    methods (Access = private)
        function warnIfIllConditioned(f)
            if f.useRcond
                if f.rcondSaved == 0.0
                    warning(message('MATLAB:singularMatrix'));
                elseif f.rcondSaved < eps(char(f.Datatype))
                    warning(message('MATLAB:nearlySingularMatrix', sprintf('%13.6e',f.rcondSaved)));
                elseif isnan(f.rcondSaved)
                    warning(message('MATLAB:illConditionedMatrix'));
                end
            else
                fRank = f.Underlying.rank_;
                if fRank < min(f.MatrixSize)
                    fTol = f.Underlying.ranktol_;
                    warning(message('MATLAB:rankDeficientMatrix',sprintf('%d',fRank),sprintf('%13.6e',fTol)));
                end
            end
        end
    end
    
    methods (Hidden)
        function displayAllProperties(obj)
            propgroups = getAllPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroups);
        end
    end
    
end

function [type, uplo, checkcond, banddensity, mindegree, ...
    luPivot, luSymPivot, ldlPivot, ranktol] = parseInputs(args)

type = 'auto';
uplo = 'none'; % values are 'none' (not assigned), 'upper' or 'lower'

% Note: not storing these in a struct to reduce overhead for small matrices
checkcond = true;
banddensity = 0.5;
mindegree = true;
luPivot = 0.1;
luSymPivot = 0.001;
ldlPivot = 0.01;
ranktol = []; % Default behavior: determine tolerance from the matrix.

if isempty(args)
    return;
end

typeList = {'auto', 'qr', 'lu', 'chol', 'ldl', 'triangular', 'diagonal', 'permutedTriangular', 'hessenberg', 'banded', 'cod'};
NameList = {'CheckCondition', 'BandDensity', 'LUPivotTolerance', 'LDLPivotTolerance', 'RankTolerance'};
uploList = {'upper', 'lower'};
typeAndNameList = {'auto', 'qr', 'lu', 'chol', 'ldl', 'triangular', 'diagonal', 'permutedTriangular', 'hessenberg', 'banded', 'cod', ...
    'CheckCondition', 'BandDensity', 'LUPivotTolerance', 'LDLPivotTolerance', 'RankTolerance'};
uploAndNameList = {'upper', 'lower', 'CheckCondition', 'BandDensity', 'LUPivotTolerance', 'LDLPivotTolerance', 'RankTolerance'};
argOffset = 0;

indType = strncmpiWithInputCheck(args{1}, typeList);
indName = strncmpiWithInputCheck(args{1}, NameList);

if nnz(indType) + nnz(indName) ~= 1
    % Let validatestring give the error message
    validatestring(args{1}, typeAndNameList);
end

if nnz(indType) == 1
    
    type = typeList{indType};
    argOffset = 1;
    
    if argOffset == length(args)
        return;
    end
    
    indName = strncmpiWithInputCheck(args{2}, NameList);
    
    if ismember(type, {'chol', 'ldl', 'triangular'})
        % Check for 'upper' / 'lower'
        indUpLo = strncmpiWithInputCheck(args{2}, uploList);
        
        if nnz(indName) + nnz(indUpLo) ~= 1
            % Let validatestring give the error message
            validatestring(args{2}, uploAndNameList);
        end
        
        if any(indUpLo)
            uplo = uploList(indUpLo);
            argOffset = 2;
            
            if length(args) >= 3
                indName = strncmpiWithInputCheck(args{3}, NameList);
            end
        end
    end
end

if argOffset == length(args)
    return;
end

for ii=1:2:(length(args)-argOffset)
    
    if ii > 1
        % In the first iteration, indName was already computed
        indName = strncmpiWithInputCheck(args{argOffset+ii}, NameList);
    end
    
    if nnz(indName) ~= 1
        % Let validatestring give the error message
        validatestring(args{argOffset+ii}, NameList);
    else
        indName = find(indName);
    end
    
    if argOffset+ii+1 > length(args)
        error(message('MATLAB:decomposition:NameWithoutValue'));
    end
    
    value = args{argOffset+ii+1};
    
    switch indName
        case 1 %'CheckCondition'
            validateattributes(value, {'logical'}, {'scalar'});
            checkcond = value;
            
        case 2 %'BandDensity'
            validateattributes(value, {'numeric'}, {'scalar', 'real', '>=', 0, '<=', 1});
            banddensity = double(full(value));
            
        case 3 %'LUPivotTolerance'
            validateattributes(value, {'numeric'}, {'vector', 'real', '>=', 0, '<=', 1});
            value = double(value);
            if isscalar(value)
                luPivot = value;
            elseif length(value) == 2
                luPivot = value(1);
                luSymPivot = value(2);
            else
                error(message('MATLAB:decomposition:InvalidLUPivot'));
            end
            
        case 4 %'LDLPivotTolerance'
            validateattributes(value, {'double'}, {'scalar', 'real', '>=', 0, '<=', 0.5});
            ldlPivot = value;
            
        case 5 %'RankTolerance'
            validateattributes(value, {'double'}, {'scalar', 'real', 'nonnegative', 'nonnan'});
            ranktol = value;
            
    end
end
end


function [decA, type, scal] = chooseDecomposition(A, ...
    banddensity, mindegree, luPivot, luSymPivot, ldlPivot, ranktol)
scal = 1;
if issparse(A)
    
    % Rectangular A -> QR decomposition
    if size(A, 1) ~= size(A, 2)
        decA = matlab.internal.decomposition.SparseQR(A, ranktol, mindegree);
        type = 'qr';
        return
    end
    
    % Special cases for banded matrices:
    [kl, ku] = bandwidth(A);
    nnzFullBand = size(A,1)*(kl+ku+1)-kl*(kl+1)/2-ku*(ku+1)/2;
    
    % Diagonal A
    if kl == 0 && ku == 0
        decA = matlab.internal.decomposition.Diagonal(A);
        type = 'diagonal';
        return
    end
    
    % Tridiagonal A (fast version of banded LU if no pivoting is needed)
    if kl==1 && ku==1 && isreal(A) && nnz(A) == nnzFullBand
        decA = matlab.internal.decomposition.Tridiagonal(A);
        type = 'banded';
        if success(decA)
            return
        end
    end
    
    % Banded A -> use banded LU decomposition
    if kl > 0 && ku > 0 && nnz(A) > nnzFullBand*banddensity
        decA = matlab.internal.decomposition.Banded(A, kl, ku);
        type = 'banded';
        return
    end
    
    % Upper triangular A -> triangular solver
    if ku == 0 && all(diag(A) ~= 0)
        decA = matlab.internal.decomposition.SparseTriangular(A, 'lower');
        type = 'triangular';
        return
    end
    
    % Lower triangular A -> triangular solver
    if kl == 0 && all(diag(A) ~= 0)
        decA = matlab.internal.decomposition.SparseTriangular(A, 'upper');
        type = 'triangular';
        return
    end
    
    % Permutation of triangular A -> permuted triangular solver
    [ist, p, q] = matlab.internal.decomposition.builtin.isPermutedTriangle(A);
    if ist
        decA = matlab.internal.decomposition.SparseTriangular(A(p, q), 'lower', p, q);
        type = 'permutedTriangular';
        return
    end
    
    % Symmetric A
    if kl == ku && ishermitian(A)
        
        % Symmetric positive definite A -> Cholesky decomposition
        d = diag(A);
        if all(d > 0)
            try
                decA = matlab.internal.decomposition.SparseCholesky(A, mindegree);
                type = 'chol';
                if success(decA)
                    return
                end
            catch
                % Continue and use LDL decomposition
            end
        elseif all(d < 0)
            try
                decA = matlab.internal.decomposition.SparseCholesky(-A, mindegree);
                type = 'chol';
                if success(decA)
                    scal = -1;
                    return
                end
            catch
                % Continue and use LDL decomposition
            end
        end
        
        % Real symmetric A -> LDL decomposition
        if isreal(A) && matlab.internal.decomposition.builtin.sparseLDLsupported()
            try
                decA = matlab.internal.decomposition.SparseLDL(A, ldlPivot, mindegree);
                type = 'ldl';
                return
            catch
                % Continue and use LU decomposition
            end
        end
    end
    
    % General square A -> LU decomposition
    decA = matlab.internal.decomposition.SparseLU(A, luPivot, luSymPivot, mindegree);
    type = 'lu';
    
else % A is dense
    
    % Rectangular A -> QR decomposition
    if size(A, 1) ~= size(A, 2)
        decA = matlab.internal.decomposition.DenseQR(A, ranktol);
        type = 'qr';
        return
    end
    
    [kl, ku] = bandwidth(A);
    
    % Lower triangular A -> triangular solve
    if ku == 0
        decA = matlab.internal.decomposition.DenseTriangular(A, 'lower');
        type = 'triangular';
        return
    end
    
    % Upper triangular A -> triangular solve
    if kl == 0
        decA = matlab.internal.decomposition.DenseTriangular(A, 'upper');
        type = 'triangular';
        return
    end
    
    % Permutation of triangular A -> permuted triangular solver
    [ist, p, q] = matlab.internal.decomposition.builtin.isPermutedTriangle(A);
    if ist
        decA = matlab.internal.decomposition.DenseTriangular(A(p, q), 'lower', p, q);
        type = 'permutedTriangular';
        return
    end
    
    % Symmetric A
    if ishermitian(A)
        
        % Symmetric positive definite A -> Cholesky decomposition
        d = diag(A);
        if all(d > 0)
            decA = matlab.internal.decomposition.DenseCholesky(A);
            type = 'chol';
            if success(decA)
                return
            end
        elseif all(d < 0)
            decA = matlab.internal.decomposition.DenseCholesky(-A);
            type = 'chol';
            if success(decA)
                scal = -1;
                return
            end
        end
        
        % Symmetric indefinite A -> LDL decomposition
        decA = matlab.internal.decomposition.DenseLDL(A);
        type = 'ldl';
        return
    end
    
    % Hessenberg A -> LU decomposition adapted for Hessenberg matrices
    if kl == 1
        decA = matlab.internal.decomposition.DenseHessenberg(A);
        type = 'hessenberg';
        return
    end
    
    % General square A -> LU decomposition
    decA = matlab.internal.decomposition.DenseLU(A);
    type = 'lu';
end
end

function index = strncmpiWithInputCheck(arg, list)
% Check that the input is either a row char vector or a scalar string. Use
% that valid input to compare it against a provided list and return the
%index. Return [] if validation failed.
index = [];
if (ischar(arg) && isrow(arg)) || (isstring(arg) && isscalar(arg))
    index = strncmpi(arg, list, max(strlength(arg),1));
end
end