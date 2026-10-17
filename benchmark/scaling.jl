using ProfileAlignments, Random, Printf

function homologs(n=300,len=300;seed=2026,duplicate_fraction=0.15)
    rng=MersenneTwister(seed)
    alphabet=collect("ACDEFGHIKLMNPQRSTVWY")
    reference=rand(rng,alphabet,len)
    distinct=max(2,round(Int,n*(1-duplicate_fraction)))
    pool=Sequence[]
    for i in 1:distinct
        sequence=copy(reference)
        for _ in 1:max(1,len÷80)
            sequence[rand(rng,eachindex(sequence))]=rand(rng,alphabet)
        end
        for position in sort(rand(rng,20:len-20,3);rev=true)
            deleteat!(sequence,position)
        end
        push!(pool,Sequence("pool_$i",String(sequence)))
    end
    [Sequence("sequence_$i",pool[mod1(i,distinct)].sequence) for i in 1:n]
end

function measure(records,strategy;threads=1)
    GC.gc()
    timed=@timed align_sequences(records;strategy,threads)
    (seconds=timed.time,bytes=timed.bytes,width=length(first(timed.value.records).sequence))
end

sizes=isempty(ARGS) ? [100,300,600] : parse.(Int,ARGS)
warmup=homologs(8,80)
for strategy in (:accurate,:balanced,:large)
    align_sequences(warmup;strategy,threads=1)
end
println("Julia threads: ",Threads.nthreads())
println("n\tunique\tstrategy\tseconds\tallocated_bytes\twidth")
for n in sizes
    records=homologs(n)
    unique_count=length(unique(r->r.sequence,records))
    for strategy in (:accurate,:balanced,:large)
        result=measure(records,strategy;threads=:auto)
        @printf "%d\t%d\t%s\t%.6f\t%d\t%d\n" n unique_count String(strategy) result.seconds result.bytes result.width
    end
end
