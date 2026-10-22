function h=ArrayChecksum(X)

%% Exact, platform-independent checksum of the contents of a numerical array
%
%   h=ArrayChecksum(X)
%
% returns h=[h1 h2 n], where h1 and h2 are two independent checksums (integers between 0 and 2^31) of the bit patterns of the
% elements of X, and n is the number of 32-bit words checksummed.
%
% (17 Oct 2026) Used by MeshKey.m to detect changes of the mesh. The bit patterns of the elements (eg the two 32-bit words of a double)
% are multiplied by weights that depend on their position, and summed modulo two primes below 2^31. All intermediate values stay below
% 2^53, so the arithmetic is exact in double precision, and h is the same on any computer and with any MATLAB version. Any change of X,
% eg in the last bit of one value, or a permutation of values, changes h, except with a probability of about 1 in 4e18. This is a
% checksum for detecting accidental changes, not a cryptographic hash.
%
% The weights are computed for blocks of 2^20 words and kept between calls (16 MB); each block is multiplied by a further,
% block-dependent factor, so that positions in different blocks are distinguished.
%
% double, single, int32, uint32, int64 and uint64 arrays are checksummed via their bit patterns; logical, char, int8, uint8, int16 and
% uint16 arrays via their values converted to double; complex arrays via [real(X(:));imag(X(:))]; sparse arrays via full(X).
%%

persistent W

B=2^20 ; Q=2^21 ; P=[2147483647 2147483629] ;
K=[1073741789 805306457 12345 ; 939524087 671088637 54321] ;

if isempty(W)
    i=(1:B)' ;
    W=zeros(B,2) ;
    for k=1:2
        W(:,k)=1+mod(mod(i.*i,Q)*K(k,1)+mod(i,Q)*K(k,2)+K(k,3),Q) ;   % weights in [1,2^21]; all terms below 2^52
    end
end

x=X(:) ;
if issparse(x) ; x=full(x) ; end
if ~isreal(x) ; x=[real(x);imag(x)] ; end
if islogical(x) || ischar(x) || any(class(x)==["int8","uint8","int16","uint16"])
    x=double(x) ;
end
if isempty(x)
    h=[0 0 0] ;
    return
end

w=double(typecast(x,'uint32')) ;          % exact 32-bit words, below 2^32
n=numel(w) ; nf=floor(n/B) ; r=n-nf*B ;          % full blocks, and the length of the last, partial block
if nf>0 ; wf=reshape(w(1:nf*B),B,nf) ; end
b=1:(nf+(r>0)) ;
h=[0 0 n] ;
for k=1:2
    s=zeros(1,numel(b)) ;
    if nf>0
        s(1:nf)=mod(sum(mod(wf.*W(:,k),P(k)),1),P(k)) ;  % products below 2^53, block sums below 2^51
    end
    if r>0
        s(end)=mod(sum(mod(w(nf*B+1:n).*W(1:r,k),P(k))),P(k)) ;
    end
    R=1+mod(mod(b.*b,Q)*K(k,2)+b*K(k,1)+K(k,3),Q) ;        % block factors in [1,2^21]
    h(k)=mod(sum(mod(s.*R,P(k))),P(k)) ;
end

end
