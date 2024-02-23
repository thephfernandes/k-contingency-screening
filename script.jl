include("cfb/cfb.jl")
using Main.CFB, BenchmarkTools

# Escolha do arquivo a ser utilizado
if length(ARGS) < 1
    println("Por favor, forneça um sistema como argumento da linha de comando.")
    exit(1)
end

sistema = ARGS[1]
arquivo = "$sistema.txt"

# Chamadas ao cálculo da CFB
println("------ $sistema -------")
g = read_edgelist(arquivo)
@time cfb_main(g, 1, sistema)
