using ProfileAlignments, Random, Printf

function homologs(n=96,len=300;seed=2026)
    rng=MersenneTwister(seed)
    alphabet=collect("ACDEFGHIKLMNPQRSTVWY")
    reference=rand(rng,alphabet,len)
    records=Sequence[]
    for i in 1:n
        sequence=copy(reference)
        for _ in 1:max(1,len÷100)
            sequence[rand(rng,eachindex(sequence))]=rand(rng,alphabet)
        end
        for position in sort(rand(rng,20:len-20,3);rev=true)
            deleteat!(sequence,position)
        end
        push!(records,Sequence("sequence_$i",String(sequence)))
    end
    records
end

n=length(ARGS)>=1 ? parse(Int,ARGS[1]) : 96
len=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 300
records=homologs(n,len)
warmup=records[1:8]
align_sequences(warmup;threads=1,retree=1,max_iterations=0)
align_sequences(warmup;threads=:auto,retree=1,max_iterations=0)

function timed(f)
    GC.gc()
    @elapsed f()
end

serial=timed() do
    align_sequences(records;threads=1,retree=1,max_iterations=0)
end
threaded=timed() do
    align_sequences(records;threads=:auto,retree=1,max_iterations=0)
end

@printf "Julia threads: %d\n" Threads.nthreads()
@printf "Input: %d sequences, ~%d residues\n" length(records) length(records[1].sequence)
@printf "Serial:   %.3f s\n" serial
@printf "Threaded: %.3f s\n" threaded
@printf "Speedup:  %.2fx\n" serial/threaded
