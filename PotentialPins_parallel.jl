include("/home/trung/_qhe-julia/FQH_state_v2.jl")
include("/home/trung/_qhe-julia/Potentials.jl")

using .FQH_states
using .Potentials
using LinearAlgebra
using SparseArrays

#Ne = 10
#No = 2*Ne-1
BLAS.set_num_threads(1)  # Set BLAS to use the same threads as main program
println("Using $(Threads.nthreads()) threads")


function get_Lz(all_basis::Vector{BitVector}, all_coef_ortho::Matrix{T} where T<: Number,θ::Float64, ϕ::Float64)
    dim = length(all_basis)
    #mat = spzeros(Complex{Float64},(dim,dim))

    Ne = count(all_basis[1])
    No = length(all_basis[1])

    #mat = diracdelta_matrix(all_basis, [pos1,pos2,pos3,pos4])
    #@time for i in 1:npins
    #    mat += sphere_point_matrix(all_basis, θ_list[i], ϕ_list[i])
    #end

    @time mat = sphere_point_matrix(all_basis, θ,ϕ) + sphere_point_matrix(all_basis, θ,ϕ+π) + 5 * sparse(I, dim, dim)
    #@time mat = sphere_point_matrix(all_basis,θ_list, ϕ_list,1.0,5.0)

    @time mat_MR = conj.(all_coef_ortho) * mat * transpose(all_coef_ortho)


    @time ED = eigen(mat_MR)

    println()

    gap = abs(ED.values[2])-abs(ED.values[1])

    #println("Eigenvalues = ")
    #if length(ED.values) > 3
    #    display(ED.values[1:5])
    #    println("\t.\n\t.\n\t.")
    #else
    #    display(ED.values)
    #end


    vecs = ED.vectors[:,sortperm(abs.(ED.values))] # Sort by increasing abs(E)


    gs_coef = transpose(all_coef_ortho) * vecs[:,1]

    gs = FQH_state(all_basis, gs_coef)
    #printwf(gs;fname="ground_states/$(Ne)e$(No)_gs_0_0_$(θ)_$(ϕ)")

    LZ = get_Lz_sphere(gs)
    return gap, LZ

end

function main()
    println("Input root file name:")
    fname = readline()
    roots = readwf(fname)
    jack_list = Vector{FQH_state}()

    for root in roots.basis
        rootstring = prod(string.((Int.(root))))
        print("\r$(rootstring)\t")
        jack = sphere_normalize(readwf("jacks/J_$(rootstring)"))
        push!(jack_list, jack)
    end

    Ne = count(jack_list[1].basis[1])
    No = length(jack_list[1].basis[1])
    println("\n\n$(Ne) electrons and $(No) orbitals\n")


    all_basis, all_coef = collate_many_vectors(jack_list; separate_out=true, collumn_vector=true)

    println("Orthonormalizing basis using QR decomposition")
    @time all_coef_ortho = collect(transpose(Matrix(qr(all_coef).Q)))
    println("------")

    gap_list = Float64[]
    LZ_list  = Float64[]


    θ_list   = π * (0:0.05:0.5)
    LZ_list  = zero(θ_list)
    gap_list = zero(θ_list)
    npts = length(θ_list)

    Threads.@threads for i in 1:npts
        println(θ_list[i])
        gap, LZ = get_Lz(all_basis, all_coef_ortho, θ_list[i], 0.0)
        gap_list[i]= gap
        LZ_list[i] = LZ
        println("-*-*-*-*-*-")
    end

    npts = length(gap_list)
    open("output_Lz/$(fname)_Lz.dat", "w") do f
        for i in 1:npts
            write(f, "$(θ_list[i])\t$(LZ_list[i])\n")
        end
    end

    open("output_gap/$(fname)_gap.dat", "w") do f
        for i in 1:npts
            write(f, "$(θ_list[i])\t$(gap_list[i])\n")
        end
    end
end

@time main()