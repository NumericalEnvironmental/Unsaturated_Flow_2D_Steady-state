########################################################################
#
# Unsat.jl
#
# by Walt McNab (May 2018)
#
# A numerical model for steady-state, unsaturated flow
# in porous media, assuming a single fluid phase and passive gas phase
#
# Model domain is a 2-D vertical cross-section
#
########################################################################


### defined types ###


type Material
    name::AbstractString
    Kh::Float64                         # (saturated) hydraulic conductivity components
    Kz::Float64
    alpha::Float64                      # Van Genuchten unsaturated hydraulic parameters
    N::Float64
    m::Float64
    sr::Float64                         # residual saturation
end


type Cell
    x::Float64                          # grid cell center coordinates
    z::Float64
    P::Float64                          # pressure head (psi)
    Q::Float64                          # source/sink flux
    matNum::Int64                       # associated material index number
    fixed::Bool                         # fixed-head cell, or not (boolean)
    xConnect::Array{Int64, 1}           # list of connecting cells
    zConnect::Array{Int64, 1}
end


type Domain
    xLength::Float64                    # model discretization
    zLength::Float64
    nx::Int64
    nz::Int64
    dx::Float64
    dz::Float64
    Ax::Float64
    Az::Float64
    defaultMat::AbstractString          # default material assignment (heterogeneities handle in blocks.txt file)
end


type Solver
    psi0::Float64                       # initial pressure across entire model domain (to start iterations)
    f::Float64                          # pressure correction weighting factor, per iteration
    maxIter::Int64                      # number of iterations
end


### support functions (physics) ###


function kr(psi, material)
    # relative permeability (Van Genuchten/Mualem model)
    if psi < 0.
        y1 = 1. + abs(material.alpha * psi)^material.N
        y2 = 1. - abs(material.alpha * psi)^(material.N - 1.) * y1^(-material.m)
        k = y2^2 / y1^(material.m/2.)
    else
        k = 1.
    end
    return k
end


function Se(psi, material)
    # effective saturation (Van Genuchten/Mualem model)
    if psi < 0.
        satEff = (1. + abs(material.alpha * psi)^material.N)^(-material.m)
    else
        satEff = 1.
    end
    return satEff
end


function S(psi, material)
    # saturation vs relative saturation
    return Se(psi, material) * (1. - material.sr) + material.sr
end


### support functions (utility) ###


function HMean2(a, b)
    # harmonic mean of two numbers
    return 2./(1./a + 1./b)
end


function GetIndex(x, z, domain)
    # find the index number of nearest cell object associated with location (x, y)
    row = maximum([round(Int64, z/domain.dz + 0.5), 1])
    col = maximum([round(Int64, x/domain.dx + 0.5), 1])
    return (row-1)*domain.nx + col
end


function GetMatNum(material, name)
    # return the index number of name in the materials list
    index = 0
    for (i, mat) in enumerate(material)
        if mat.name == name
            index = i
            break
        end
    end
    return index
end


function CreateCells(domain, matNum, psi0)
    # create cell objects with default material assignment and initial pressure head
    cell = Cell[]
    for j = 1:domain.nz, i = 1:domain.nx
        x = (i-0.5) * domain.dx
        z = (j-0.5) * domain.dz
        xConnect = Int64[]
        zConnect = Int64[]
        push!(cell, Cell(x, z, psi0, 0.0, matNum, false, xConnect, zConnect))
    end
    # populate cell connection arrays
    for j = 0:domain.nz-1, i = 1:domain.nx-1
        push!(cell[i+j*domain.nx].xConnect, i+j*domain.nx+1)
        push!(cell[i+j*domain.nx+1].xConnect, i+j*domain.nx)
    end
    for j = 0:domain.nz-2, i = 1:domain.nx
        push!(cell[i+j*domain.nx].zConnect, i+(j+1)*domain.nx)
        push!(cell[i+(j+1)*domain.nx].zConnect, i+j*domain.nx)
    end
    return cell
end


function UpdatePressures(cell, newP, solverParams)
    for (i, ce) in enumerate(cell)
        # weighting between old modeled head and new modeled head, per iteration
        ce.P = solverParams.f*newP[i] + (1.-solverParams.f)*ce.P        
    end
    return cell
end


function SumFluxes(cell, material, domain)
    # compute fluxes (for model output and final flow balance estimates)
    flowSum = Float64[]
    Qx = Float64[]
    Qz = Float64[]
    errorF = Float64[]
    for ce in cell
        push!(flowSum, ce.Q)                # initialize summations
        push!(Qx, 0.)
        push!(Qz, 0.)
    end
    for (i, ce) in enumerate(cell)
        matI = material[ce.matNum]        
        for xcon in ce.xConnect         # for each horizontal connection (up to 2 columns)
            matJ = material[cell[xcon].matNum]
            conduct = HMean2(kr(ce.P, matI)*matI.Kh, kr(cell[xcon].P, matJ)*matJ.Kh) * domain.Ax / domain.dx
            flux = conduct * (cell[xcon].P - ce.P)
            flowSum[i] += flux
            Qx[i] += sign(ce.x - cell[xcon].x) * flux
        end       
        Qx[i] = Qx[i]/length(ce.xConnect)
        for zcon in ce.zConnect         # for each vertical connection (up to two columns)
            matJ = material[cell[zcon].matNum]
            conduct = HMean2(kr(ce.P, matI)*matI.Kz, kr(cell[zcon].P, matJ)*matJ.Kz) * domain.Az / domain.dz
            flux = conduct * ((cell[zcon].P + cell[zcon].z) - (ce.P + ce.z))
            flowSum[i] += flux
            Qz[i] += sign(ce.z - cell[zcon].z) * flux            
        end  
        Qz[i] = Qz[i]/length(ce.zConnect)
        flowThru = abs(Qx[i]) + abs(Qz[i])
        push!(errorF, flowSum[i]/flowThru)
    end
    return Qx, Qz, errorF
end


### matrix operations functions


function LHS_matrix(cell, material, domain)
    # fill out the LHS of the equation matrix
    row_index = Int64[]                     # indexing system for sparse matrix
    col_index = Int64[]
    data = Float64[]
    for (i, ce) in enumerate(cell)          # for each row
        diag = 0.0
        if ce.fixed == false
            matI = material[ce.matNum]        
            for xcon in ce.xConnect         # for each horizontal connection (up to 2 columns)
                matJ = material[cell[xcon].matNum]
                conduct = HMean2(kr(ce.P, matI)*matI.Kh, kr(cell[xcon].P, matJ)*matJ.Kh) * domain.Ax / domain.dx
                push!(row_index, i)
                push!(col_index, xcon)
                push!(data, conduct)
                diag -= conduct
            end       
            for zcon in ce.zConnect         # for each vertical connection (up to two columns)
                matJ = material[cell[zcon].matNum]           
                conduct = HMean2(kr(ce.P, matI)*matI.Kz, kr(cell[zcon].P, matJ)*matJ.Kz) * domain.Az / domain.dz 
                push!(row_index, i)
                push!(col_index, zcon)
                push!(data, conduct)                
                diag -= conduct
            end   
        else
            diag = 1.0             # fixed pressure head cell
        end    
        # set diagonal term
        push!(row_index, i)
        push!(col_index, i)
        push!(data, diag)
    end
    return data, row_index, col_index
end


function RHS_vector(cell, material, domain)
    # construct explicit matrix
    b = Float64[]
    for (i, ce) in enumerate(cell)          # for each row    
        if ce.fixed == false
            bSum = -ce.Q
            matI = material[ce.matNum]        
            for zcon in ce.zConnect         # for each vertical connection (up to two columns)
                matJ = material[cell[zcon].matNum]           
                conduct = HMean2(kr(ce.P, matI)*matI.Kz, kr(cell[zcon].P, matJ)*matJ.Kz) * domain.Az / domain.dz 
                bSum -= conduct * (cell[zcon].z - ce.z)
            end  
            push!(b, bSum)
        else
            push!(b, ce.P)                  # fixed-head cell
        end
    end
    return b
end


### input and output functions


function ReadMaterials()
    # read various model parameters from file
    material = Material[]
    data = readdlm("materials.txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        name = data[1][i, 1]
        Kh = Float64(data[1][i, 2])
        Kz = Float64(data[1][i, 3])
        alpha = Float64(data[1][i, 4])
        N = Float64(data[1][i, 5])
        m = 1-1/N
        sr = Float64(data[1][i, 6])
        push!(material, Material(name, Kh, Kz, alpha, N, m, sr))
    end
    println("Read material properties.")
    return material
end


function ReadSolverParams()
    # read numerical model "knobs" from file
    data = readdlm("solver.txt", '\t', header=false)
    psi0 = Float64(data[1, 2])
    f = Float64(data[2, 2])
    maxIter = Int64(data[3, 2])
    solverParams = Solver(psi0, f, maxIter)
    println("Read solver parameters.")
    return solverParams
end


function ReadDomain()
    # model geometry
    data = readdlm("domain.txt", '\t', header=false)
    xLength = Float64(data[1, 2])
    zLength = Float64(data[2, 2])
    nx = Int64(data[3, 2])
    nz = Int64(data[4, 2])
    dx = xLength/nx
    dz = zLength/nz
    Ax = dz             # placeholders for later expansion to 3-D
    Az = dx
    defaultMat = data[5, 2]
    domain = Domain(xLength, zLength, nx, nz, dx, dz, Ax, Az, defaultMat)
    println("Read model domain properties.")
    return domain
end


function ReadBlocks(cell, material)
    # read and process rectangular material property heterogeneities
    data = readdlm("blocks.txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        x0 = Float64(data[1][i, 1])
        z0 = Float64(data[1][i, 2])
        xf = Float64(data[1][i, 3])
        zf = Float64(data[1][i, 4])
        matBlock = data[1][i, 5]
        matNum = GetMatNum(material, matBlock)
        for ce in cell
            if ce.x >= x0 && ce.x <= xf && ce.z >= z0 && ce.z <= zf
                ce.matNum = matNum
            end
        end
    end
    println("Read material blocks.")
    return cell
end


function ReadSources(cell, domain)
    # read and process line sources/boundary conditions
    data = readdlm("sources.txt", '\t', header=true)
    for i = 1:size(data[1], 1)
        BC = data[1][i, 1]
        direction = data[1][i, 2]
        pos = Float64(data[1][i, 3])
        start = Float64(data[1][i, 4])
        fin = Float64(data[1][i, 5])
        value = Float64(data[1][i, 6])
        if direction == "x"
            cellStart = GetIndex(start, pos, domain)
            cellEnd = GetIndex(fin, pos, domain)
            stepSize = 1
            A = domain.Az
        else
            cellStart = GetIndex(pos, start, domain)
            cellEnd = GetIndex(pos, fin, domain)
            stepSize = domain.nx
            A = domain.Ax
        end
        if BC == "flux"
            for i = cellStart:stepSize:cellEnd
                cell[i].Q = value * A
            end
        else
            for i = cellStart:stepSize:cellEnd
                cell[i].fixed = true
                cell[i].P = value
            end
        end
    end
    println("Read sources/boundary condition constraints.")
    return cell
end


function Spatial(domain, material, solverParams)
    # set up the model spatial domain
    matNum = GetMatNum(material, domain.defaultMat)                 # default material index number
    cell = CreateCells(domain, matNum, solverParams.psi0)           # create grid cells
    cell = ReadBlocks(cell, material)                               # read blocks file and update cell material assignments, as warranted
    cell = ReadSources(cell, domain)                                # read and assign sources/boundary conditions
    println("Set up model cells and distributed properties.")
    return cell
end


function WriteOutput(cell, material, domain, fileName)
    # write property distribution to file
    csvfile = open(fileName,"w")
    flowImbalance = Float32[]
    Qx, Qz, errorF = SumFluxes(cell, material, domain)              # update flow balance estimates
    line_out = "cell" * "," * "x" * "," * "z" * "," * "psi" * "," * "S" * "," * "Qx" * "," * "Qz" * "," * "error"
    println(csvfile, line_out)
    for (i, ce) in enumerate(cell)
        sat = S(ce.P, material[ce.matNum])
        line_out = string(i) * "," * string(ce.x) * "," * string(ce.z) * "," * string(ce.P) * "," * string(sat) * "," * string(Qx[i]) * "," * string(Qz[i]) * "," * string(errorF[i])
        println(csvfile,line_out)
    end
    close(csvfile)
    println("Wrote " * fileName * ".")
    return mean(abs(errorF)), maximum(abs(errorF))
end


### main ###


function Unsat()

    # problem setup
    domain = ReadDomain()
    material = ReadMaterials()
    solverParams = ReadSolverParams()
    cell = Spatial(domain, material, solverParams)
    
    # write initial conditions to file
    fileName = "starting_estimates.csv"
    WriteOutput(cell, material, domain, fileName)

    println("Solving ...")
    
    # iterative solution of coupled non-linear flow balance PDEs
    numCells = length(cell)    
    iter = 0
    while iter < solverParams.maxIter
        data, row_index, col_index = LHS_matrix(cell, material, domain)         # assemble and solve matrix of flux balance equations
        A = sparse(row_index, col_index, data, numCells, numCells) 
        b = RHS_vector(cell, material, domain) 
        newP = \(A, b)
        cell = UpdatePressures(cell, newP, solverParams)
        iter += 1
    end
    
    # write model results to file
    fileName = "final_estimates.csv"
    errorMean, errorMax = WriteOutput(cell, material, domain, fileName)
    println("Mean absolute flow imbalance = ", errorMean)
    println("Max. absolute flow imbalance = ", errorMax)    
    println("Finished.")

end


Unsat()             ### run model ###
