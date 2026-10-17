using Test, Random, ProfileAlignments, BioSequences
const PA=ProfileAlignments
ungapped(a)=replace.(getfield.(a.records,:sequence),"-"=>"")

# Independent oracle: enumerate every monotone alignment path for tiny inputs.
function exhaustive_score(a,b,mode;match=2.0,mismatch=-3.0,open=-4.0,extend=-1.0,terminal=-0.5)
    best=-Inf
    function score(x,y)
        total=0.0; prev=0
        nongapx=findall(!=('-'),x); nongapy=findall(!=('-'),y)
        spanx=isempty(nongapx) ? (length(x)+1,0) : extrema(nongapx)
        spany=isempty(nongapy) ? (length(y)+1,0) : extrema(nongapy)
        for k in eachindex(x)
            if x[k]!='-' && y[k]!='-'
                total+=x[k]==y[k] ? match : mismatch
                prev=0
            else
                state=x[k]=='-' ? 1 : 2
                free=mode==:overlap ? !(spanx[1]<=k<=spanx[2] && spany[1]<=k<=spany[2]) :
                    mode==:query && state==1 && !(spanx[1]<=k<=spanx[2])
                total+=free ? terminal : prev==state ? extend : open
                prev=state
            end
        end
        total
    end
    function visit(i,j,x,y)
        if i>length(a) && j>length(b)
            best=max(best,score(x,y)); return
        end
        i<=length(a) && j<=length(b) && visit(i+1,j+1,vcat(x,a[i]),vcat(y,b[j]))
        i<=length(a) && visit(i+1,j,vcat(x,a[i]),vcat(y,'-'))
        j<=length(b) && visit(i,j+1,vcat(x,'-'),vcat(y,b[j]))
    end
    visit(1,1,Char[],Char[])
    best
end

@testset "Affine optimality against exhaustive paths" begin
    words=[""]
    for n in 1:3
        append!(words,[String(collect(chars)) for chars in Iterators.product(fill(['A','C'],n)...)])
    end
    for mode in (:global,:overlap,:query)
        model=ScoringModel(;alphabet=:dna,gap_open=-4,gap_extend=-1,terminal_gap=-0.5,mode)
        for a in words,b in words
            result=pairalign(a,b;model)
            @test result.score≈exhaustive_score(a,b,mode)
            @test ungapped(result.alignment)==[a,b]
            @test PA.affine_path(Alignment([Sequence("a",a)]),Alignment([Sequence("b",b)]),model)[2]≈result.score
        end
    end
end

@testset "Scoring and validation" begin
    model=ScoringModel()
    @test model.matrix[PA.residue_index(model,'W'),PA.residue_index(model,'W')]==11
    @test model.matrix[PA.residue_index(model,'D'),PA.residue_index(model,'E')]==2
    @test model.matrix[PA.residue_index(model,'I'),PA.residue_index(model,'L')]==2
    @test pairalign("NNN","NNN";model=ScoringModel(;alphabet=:dna,mode=:global)).score==0
    @test_throws ArgumentError ScoringModel(;gap_open=1)
    @test_throws ArgumentError ScoringModel(;gap_extend=NaN)
    @test_throws ArgumentError ScoringModel(;mode=:wrong)
    @test_throws ArgumentError ScoringModel(;matrix=ones(2,2))
    @test_throws ArgumentError Alignment([Sequence("a","AC"),Sequence("a","AC")])
    @test_throws ArgumentError Alignment([Sequence("a","A"),Sequence("b","AC")])
    @test_throws ArgumentError Alignment([Sequence("a","A")];weights=[0])
    @test_throws ArgumentError pairalign("AC!","AC")
    @test_throws ArgumentError pairalign("AC","AC";max_cells=1)
    globalmodel=ScoringModel(;alphabet=:dna,mode=:global,gap_open=-4,gap_extend=-1)
    long=pairalign("AAAAACCCCCGGGGG","AAAAAGGGGG";model=globalmodel)
    @test long.score==10*2-4-4*1
    @test count(==('-'),long.alignment[2].sequence)==5
    fragment=pairalign("CDEFG","ACDEFGHIK";model=ScoringModel(;mode=:query))
    @test fragment.alignment[1].sequence=="-CDEFG---"
    @test fragment.score==sum(model.matrix[PA.residue_index(model,c),PA.residue_index(model,c)] for c in "CDEFG")
    bio_pair=pairalign(aa"ACDEFG",aa"CDEFG";model=ScoringModel(;mode=:overlap))
    @test ungapped(bio_pair.alignment)==["ACDEFG","CDEFG"]
    @test Sequence("dna",dna"ACGN-").sequence=="ACGN-"
    @test ungapped(align_sequences([aa"ACDEFG",aa"CDEFG"]))==["ACDEFG","CDEFG"]
end

function column_sequences(a)
    [Tuple(col) for col in eachcol(PA.alignment_matrix(a)) if any(!=('-'),col)]
end
@testset "Profile column preservation and weights" begin
    a=Alignment([Sequence("a","ACD-EFGH"),Sequence("b","ACDKEFGH")];weights=[1,3])
    b=Alignment([Sequence("c","CDEFG")])
    merged=merge_profiles(a,b)
    @test ungapped(merged)==["ACDEFGH","ACDKEFGH","CDEFG"]
    @test column_sequences(Alignment(merged.records[1:2]))==column_sequences(a)
    @test merged.weights[2]/merged.weights[1]==3
    @test all(isfinite,Alignment(a.records;weights=[1e308,1e308]).weights)
    @test sum(merged.weights)≈1
    @test merge_profiles(Alignment(Sequence[]),a)===a
    doublegap=Alignment([Sequence("a","AC--DE"),Sequence("b","AC--DE")])
    clean=Alignment([Sequence("a","ACDE"),Sequence("b","ACDE")])
    @test alignment_score(doublegap)≈alignment_score(clean)
    @test PA.trim_columns(doublegap).records[1].sequence=="ACDE"
    @test eltype(PA.alignment_matrix(a))==UInt8
end

@testset "FFT against direct transforms and correlations" begin
    rng=MersenneTwister(107)
    for n in (1,2,4,8,16,32)
        x=complex.(randn(rng,n),randn(rng,n))
        direct=[sum(x[j+1]*cis(-2*pi*j*k/n) for j in 0:n-1) for k in 0:n-1]
        transformed=PA.fft!(copy(x))
        @test transformed≈direct atol=1e-10
        @test PA.fft!(transformed;inverse=true)≈x atol=1e-10
    end
    @test_throws ArgumentError PA.fft!(zeros(ComplexF64,3))
    a=Alignment([Sequence("a","ACDEFGHIKLMNPQRSTVWY")])
    b=Alignment([Sequence("b","DEFGHIKLMNPQ")])
    peaks=fft_lags(a,b;max_peaks=100,min_overlap=1)
    av,bv=PA.encoding(a,:protein),PA.encoding(b,:protein)
    for p in peaks
        direct=sum(av[c,i]*bv[c,i+p.lag] for c in axes(av,1),i in max(1,1-p.lag):min(size(av,2),size(bv,2)-p.lag))
        @test p.correlation≈direct atol=1e-9
    end
    @test first(fft_lags(a,b)).lag==-2
    anchors=fft_anchors(a,b;window=4)
    @test !isempty(anchors)
    @test all(j->anchors[j][1]<anchors[j+1][1] && anchors[j][2]<anchors[j+1][2],1:length(anchors)-1)
    approximate=merge_profiles(a,b;anchors=:fft,anchor_window=4)
    @test ungapped(approximate)==ungapped(a) ∪ ungapped(b)
    @test approximate[2].sequence=="--DEFGHIKLMNPQ------"
    dna_a=Alignment([Sequence("a","AACGTACG")]); dna_b=Alignment([Sequence("b","CGTAC")])
    @test first(fft_lags(dna_a,dna_b;alphabet=:dna)).lag==-2
end

@testset "Guide trees, duplicates, refinement and MSA" begin
    @test guide_tree([0.0 0.1 1.0;0.1 0.0 1.0;1.0 1.0 0.0]).left.leaves==[3]
    @test_throws ArgumentError guide_tree([0.0 1;2 0])
    @test_throws ArgumentError guide_tree(zeros(0,0))
    @test_throws ArgumentError kmer_distances([Sequence("a","A")];k=0)
    legacy_kmers=function(records;k=6,alphabet=:protein)
        raw=[replace(r.sequence,"-"=>"") for r in records]
        effective=isempty(raw) ? k : min(k,max(1,minimum(length.(raw))))
        counts=[PA.tuple_kmer_counts(s,effective,alphabet) for s in raw]
        PA.distances_from_kmers(raw,counts,true,1)
    end
    kmer_fixture=[Sequence("k$i",s) for (i,s) in enumerate((
        "ACDEFGHIKLMNPQ","ACDEYGHIKLMNPQ","NNNNNNNNNNNNNN","ACD---HIKLMNPQ"))]
    for k in (1,6)
        @test kmer_distances(kmer_fixture;k,threads=1)==legacy_kmers(kmer_fixture;k)
    end
    long_kmers=[Sequence("long_$i",repeat(s,3)) for (i,s) in enumerate((
        "ACDEFGHIKLMNPQ","ACDEYGHIKLMNPQ","NNNNNNNNNNNNNN"))]
    for k in (22,23)
        @test kmer_distances(long_kmers;k,threads=1)==legacy_kmers(long_kmers;k)
    end
    records=[Sequence("full","ACDEFGHIKLMNPQ"),Sequence("fragment","DEFGHIK"),
        Sequence("insert","ACDEFGWHIKLMNPQ"),Sequence("duplicate","ACDEFGHIKLMNPQ")]
    distance=kmer_distances(records)
    @test distance[1,4]==0
    @test distance==transpose(distance)
    result=align_sequences(records;return_trace=true)
    @test getfield.(result.alignment.records,:id)==getfield.(records,:id)
    @test ungapped(result.alignment)==getfield.(records,:sequence)
    @test result.alignment[1].sequence==result.alignment[4].sequence
    @test all(diff(result.scores).>=0)
    @test last(result.scores)≈alignment_score(result.alignment)
    @test length(result.representatives)==3
    @test result.alignment.weights[1]==result.alignment.weights[4]
    @test result.strategy==:accurate
    @test align_sequences(records).records[1].sequence==result.alignment[1].sequence
    balanced=align_sequences(records;strategy=BalancedStrategy(),return_trace=true)
    large=align_sequences(records;strategy=:large,return_trace=true)
    @test balanced.strategy==:balanced
    @test large.strategy==:large
    @test large.cycles==0
    @test ungapped(balanced.alignment)==getfield.(records,:sequence)
    @test ungapped(large.alignment)==getfield.(records,:sequence)
    automatic=align_sequences(records;strategy=:auto,strategy_thresholds=(1,2),return_trace=true)
    @test automatic.strategy==:large
    @test align_sequences(records;strategy=:large,anchors=:none,retree=2,
        max_iterations=1,return_trace=true).cycles>=1
    @test_throws ArgumentError align_sequences(records;strategy=:unknown)
    @test_throws ArgumentError align_sequences(records;strategy=:auto,strategy_thresholds=(10,2))
    @test isempty(align_sequences(Sequence[]))
    @test align_sequences(["AC"]).records[1].sequence=="AC"
    @test align_sequences(["","AC"]).records[1].sequence=="--"
    # A deliberately shifted internal block has a strictly improvable score.
    bad=Alignment([Sequence("a","ACDEFGHIK--"),Sequence("b","ACD--EFGHIK"),Sequence("c","ACDEFGHIK--")])
    refined=refine_alignment(bad;model=ScoringModel(;mode=:global),return_trace=true,max_iterations=5)
    @test refined.accepted>0
    @test all(diff(refined.scores).>0)
    @test ungapped(refined.alignment)==ungapped(bad)
    @test all(r->r.sequence=="ACDEFGHIK",refined.alignment.records)
    # Stress varied substitutions, indels and terminal fragments reproducibly.
    rng=MersenneTwister(41)
    for trial in 1:12
        reference=rand(rng,collect("ACDEFGHIKLMNPQRSTVWY"),35)
        seqs=[String(reference)]
        for j in 1:5
            mutated=copy(reference)
            mutated[rand(rng,1:35)]=rand(rng,collect("ACDEFGHIKLMNPQRSTVWY"))
            deleteat!(mutated,rand(rng,10:20))
            push!(seqs,String(mutated))
        end
        out=align_sequences(seqs;return_trace=true)
        @test ungapped(out.alignment)==seqs
        @test all(diff(out.scores).>=0)
        @test length(unique(length.(getfield.(out.alignment.records,:sequence))))==1
    end
    dna=align_sequences(["ACGTACGT","CGTAC","ACGTTACGT"];model=ScoringModel(;alphabet=:dna))
    @test ungapped(dna)==["ACGTACGT","CGTAC","ACGTTACGT"]
end

@testset "Deterministic threaded execution" begin
    rng=MersenneTwister(812)
    reference=rand(rng,collect("ACDEFGHIKLMNPQRSTVWY"),70)
    records=Sequence[Sequence("reference",String(reference))]
    for i in 1:11
        sequence=copy(reference)
        sequence[rand(rng,1:length(sequence))]=rand(rng,collect("ACDEFGHIKLMNPQRSTVWY"))
        deleteat!(sequence,rand(rng,15:55))
        push!(records,Sequence("sample_$i",String(sequence)))
    end
    serial=align_sequences(records;threads=1,return_trace=true)
    threaded=align_sequences(records;threads=:auto,return_trace=true)
    @test getfield.(threaded.alignment.records,:sequence)==
        getfield.(serial.alignment.records,:sequence)
    @test threaded.alignment.weights==serial.alignment.weights
    @test threaded.scores==serial.scores
    @test kmer_distances(records;threads=:auto)==kmer_distances(records;threads=1)
    @test alignment_score(threaded.alignment;threads=:auto)==
        alignment_score(threaded.alignment;threads=1)
    @test_throws ArgumentError align_sequences(records;threads=0)
    @test_throws ArgumentError kmer_distances(records;threads=:invalid)
    groups=[records[1:4],records[5:8],records[9:12]]
    batched=align_sequence_sets(groups;threads=:auto,strategy=:large)
    serial_groups=[align_sequences(group;threads=1,strategy=:large) for group in groups]
    @test [getfield.(a.records,:sequence) for a in batched]==
        [getfield.(a.records,:sequence) for a in serial_groups]
    @test isempty(align_sequence_sets(Vector{Vector{Sequence}}();threads=:auto))
end

@testset "FASTA standalone workflow" begin
    mktempdir() do dir
        path=joinpath(dir,"in.fasta"); output=joinpath(dir,"out.fasta")
        records=[Sequence("one complete ID","ACDEFGHIK"),Sequence("two","CDEFG")]
        write_fasta(records,path;width=4)
        @test getfield.(read_fasta(path),:id)==getfield.(records,:id)
        result=align_fasta(path,output)
        @test getfield.(read_fasta(output),:sequence)==getfield.(result.records,:sequence)
        write(path,"ACGT\n")
        @test_throws ArgumentError read_fasta(path)
    end
end
