using ProfileAlignments

records=[Sequence("reference","ACDEFGHIKLMNPQRSTVWY"),
         Sequence("fragment","DEFGHIKLMNPQ"),
         Sequence("insertion","ACDEFGWHIKLMNPQRSTVWY"),
         Sequence("duplicate","ACDEFGHIKLMNPQRSTVWY")]
result=align_sequences(records;return_trace=true)
for record in result.alignment
    println(record.id, "\t",record.sequence)
end
println("Accepted scores: ",result.scores)
println("Refinement improvements: ",result.accepted)

a=Alignment([records[1]])
b=Alignment([records[2]])
println("FFT candidate diagonals: ",fft_lags(a,b;max_peaks=3))
