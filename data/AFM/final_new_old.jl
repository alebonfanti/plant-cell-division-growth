using CSV
using DataFrames
using Plots
using DelimitedFiles
using Statistics
using StatsPlots
using HypothesisTests
using LaTeXStrings




function main()

    p = plot()


    font(family= "san-serif")

    young_marchantia = readdlm("Young_marchantia.txt", '\t');
    Fst_marchantia = readdlm("1st_division_marchantia.txt", '\t');
    Snd_marchantia = readdlm("2nd_division_marchantia.txt", '\t');

    young_leaf = readdlm("Young_leaf.txt", '\t');
    mature_leaf = readdlm("Mature_leaf.txt", '\t');
print(young_marchantia)

    boxplot!(p, [1], vec(young_marchantia), fillalpha=0.1, linewidth=2, color = :grey, legend = false)
    #dotplot!(p, [1], vec(young_marchantia), marker=(:orange, stroke(0)), legend = false)

    boxplot!(p, [2], vec(Fst_marchantia), fillalpha=0.1, linewidth=2, color = :grey, legend = false)
    #dotplot!(p, [2], vec(Fst_marchantia), marker=(:orange, stroke(0)), legend = false)

    boxplot!(p, [3], vec(Snd_marchantia), fillalpha=0.1, linewidth=2, color = :grey, legend = false, grid = :off, framestyle = :box)
    #dotplot!(p, [3],  vec(Snd_marchantia), marker=(:orange, stroke(0)), legend = false)

    boxplot!(p, [4], vec(young_leaf), fillalpha=0.1, linewidth=2, color = :grey, legend = false)
    #dotplot!(p, [4],  vec(young_leaf), marker=(:orange, stroke(0)), legend = false)

    boxplot!(p, [5], vec(mature_leaf), fillalpha=0.1, linewidth=2, color = :grey, legend = false)
    #dotplot!(p, [5],  vec(mature_leaf), marker=(:orange, stroke(0)), legend = false)

    plot!(p, [0, 6], [1, 1], linewidth=2, color = :grey, linestyle = :dash)

    #dotplot!(p, [4], conv_area', marker=(:black, stroke(0)), legend = false)
    ylabel!("Young's modulus ratio"); 
    xticks!(1:5, ["3-6h", "16-24h", "Second division", "3-6h", "16-24h"], xrotation=50)
    ylims!(0.4, 1.6)
    plot!(size=(500,300))
    savefig(p, "old_new.pdf")


    print("Number of samples \n")
    print("Marchantia 1st young = ", size(young_marchantia,2), "\n")
    print("Marchantia 1st old = ", size(Fst_marchantia,2), "\n")
    print("Marchantia 2nd young = ", size(Snd_marchantia,2), "\n")
    print("Leaf young = ", size(young_leaf,2), "\n")
    print("Mature young = ", size(mature_leaf,2), "\n")

    Fst_marchantia = vec(Fst_marchantia)
    Snd_marchantia = vec(Snd_marchantia)
    young_marchantia = vec(young_marchantia)
    young_leaf = vec(young_leaf)
    mature_leaf = vec(mature_leaf)

    Fst_Snd = MannWhitneyUTest(Fst_marchantia, Snd_marchantia)
    Fst_young = MannWhitneyUTest(Fst_marchantia, young_marchantia)
    young_Fst_leaf = MannWhitneyUTest(young_leaf, mature_leaf)

    print(Fst_young, "\n")
    print(Fst_Snd, "\n")
    print(young_Fst_leaf, "\n")


end

main()