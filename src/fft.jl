# Radix-2 Cooley–Tukey FFT; no FFTW or binary artifact dependency.
function fft!(x::Vector{ComplexF64};inverse=false)
    n=length(x)
    n>0 && ispow2(n) || throw(ArgumentError("FFT length must be a positive power of two"))
    j=0
    for i in 0:n-1
        i<j && ((x[i+1],x[j+1])=(x[j+1],x[i+1]))
        bit=n>>1
        while bit>0 && (j & bit)!=0
            j ⊻= bit; bit >>= 1
        end
        j ⊻= bit
    end
    len=2
    while len<=n
        root=cis((inverse ? 2 : -2)*pi/len)
        for start in 1:len:n
            w=1.0+0.0im
            for k in 0:(len>>1)-1
                u=x[start+k]; v=w*x[start+k+(len>>1)]
                x[start+k]=u+v; x[start+k+(len>>1)]=u-v
                w*=root
            end
        end
        len*=2
    end
    inverse && (x ./= n)
    x
end

# Grantham volume and polarity, same biochemical encoding used in MAFFT's
# original FFT scheme; normalized across the 20 standard amino acids.
const PHYSICOCHEMICAL=let
    chars=collect("SRLPTAVGIFYCHQNKDEMW")
    volumes=[32,124,111,32.5,61,31,84,3,111,132,136,55,96,85,56,119,54,83,105,170]
    polarities=[9.2,10.5,4.9,8,8.6,8.1,5.9,9,5.2,5.2,6.2,5.5,10.4,10.5,11.6,11.3,13,12.3,5.7,5.4]
    normalize(v)=((v .- sum(v)/length(v)) ./ sqrt(sum((v .- sum(v)/length(v)).^2)/length(v)))
    vn,pn=normalize(volumes),normalize(polarities)
    Dict(c=>(vn[i],pn[i]) for (i,c) in enumerate(chars))
end

function encoding(a::Alignment,alphabet)
    channels=alphabet==:protein ? 2 : 4
    values=zeros(channels,alignment_width(a))
    for (i,r) in enumerate(a.records),(j,c) in enumerate(r.sequence)
        if alphabet==:protein
            v=get(PHYSICOCHEMICAL,c,(0.0,0.0))
            values[1,j]+=a.weights[i]*v[1]; values[2,j]+=a.weights[i]*v[2]
        else
            k=findfirst(==(c),"ACGT")
            isnothing(k) || (values[k,j]+=a.weights[i])
        end
    end
    values
end

"""FFT correlations of profile volume/polarity (protein) or base indicators (DNA).
Return `(lag, correlation, overlap)` peaks. Lag is second-column minus first-column.
This identifies candidate diagonals, not an alignment or a homology probability.
"""
function fft_lags(a::Alignment,b::Alignment;alphabet=:protein,max_peaks=20,min_overlap=nothing)
    alphabet in (:protein,:dna) || throw(ArgumentError("Invalid alphabet"))
    max_peaks>=0 || throw(ArgumentError("max_peaks must be nonnegative"))
    n,m=alignment_width(a),alignment_width(b)
    (n==0 || m==0 || max_peaks==0) && return NamedTuple{(:lag,:correlation,:overlap),Tuple{Int,Float64,Int}}[]
    required=isnothing(min_overlap) ? max(1,min(8,min(n,m))) : min_overlap
    required>0 || throw(ArgumentError("min_overlap must be positive"))
    av,bv=encoding(a,alphabet),encoding(b,alphabet)
    sizefft=nextpow(2,n+m-1)
    correlations=zeros(n+m-1)
    x=zeros(ComplexF64,sizefft); y=similar(x)
    for channel in axes(av,1)
        fill!(x,0); fill!(y,0)
        x[1:n]=reverse(av[channel,:]); y[1:m]=bv[channel,:]
        fft!(x); fft!(y); x .*= y; fft!(x;inverse=true)
        @inbounds for i in eachindex(correlations)
            correlations[i]+=real(x[i])
        end
    end
    peaks=NamedTuple{(:lag,:correlation,:overlap),Tuple{Int,Float64,Int}}[]
    for lag in -(n-1):m-1
        overlap=max(0,min(n,m-lag)-max(1,1-lag)+1)
        overlap>=required && push!(peaks,(lag=lag,correlation=correlations[n+lag],overlap=overlap))
    end
    sort!(peaks;by=p->(-p.correlation,abs(p.lag),p.lag))
    peaks[1:min(max_peaks,length(peaks))]
end

function consensus(a::Alignment,model)
    p=profiles(a,model)
    [sum(col)==0 ? '-' : model.symbols[argmax(col)] for col in eachcol(p)]
end

"""Find a monotone chain of high-identity window centres on FFT candidate diagonals.
Anchors force column matches and are approximate; full DP is the default engine.
"""
function fft_anchors(a::Alignment,b::Alignment;model=ScoringModel(),window=12,
                     identity_threshold=0.85,max_peaks=20)
    window>=2 || throw(ArgumentError("window must be >=2"))
    0<identity_threshold<=1 || throw(ArgumentError("Invalid anchor identity threshold"))
    ac,bc=consensus(a,model),consensus(b,model)
    candidates=Tuple{Int,Int}[]
    for p in fft_lags(a,b;alphabet=model.alphabet,max_peaks,min_overlap=window)
        lo=max(1,1-p.lag); hi=min(length(ac),length(bc)-p.lag)
        for start in lo:window:hi-window+1
            matches=count(k->ac[k]!='-' && ac[k]==bc[k+p.lag],start:start+window-1)
            centre=start+(window-1)÷2
            matches/window>=identity_threshold && ac[centre]!='-' && ac[centre]==bc[centre+p.lag] &&
                push!(candidates,(centre,centre+p.lag))
        end
    end
    sort!(unique!(candidates))
    isempty(candidates) && return candidates
    lengthbest=ones(Int,length(candidates)); predecessor=zeros(Int,length(candidates))
    for i in eachindex(candidates),j in 1:i-1
        if candidates[j][1]<candidates[i][1] && candidates[j][2]<candidates[i][2] && lengthbest[j]+1>lengthbest[i]
            lengthbest[i]=lengthbest[j]+1; predecessor[i]=j
        end
    end
    chain=Tuple{Int,Int}[]; i=argmax(lengthbest)
    while i>0
        push!(chain,candidates[i]); i=predecessor[i]
    end
    reverse!(chain)
end
