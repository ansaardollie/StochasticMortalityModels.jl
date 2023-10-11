export expected_lifetime, lexpectancies

function expected_lifetime(
    mx, ages;
    sex::Sex,
    at_age=[0],
    mode::CalculationChoice=CC_JULIA,
    radix=1_000_000,
    debug=false)

    Δx = ages[2:end] - ages[1:end-1]
    single_ages = all(Δx .== 1)

    function a0(m0, sex::Sex, mode::CalculationChoice)
        config = a0_config(mode, sex)

        (lb, slope, intercept) = config

        idx = mode == CC_JULIA ? 3 : 2
        for i in eachindex(lb)
            if m0 >= lb[i]
                idx = i
                break
            end
        end

        return intercept[idx] + slope[idx] * m0
    end

    n = length(mx)

    ax = single_ages ? Vector{Any}(fill(0.5, n)) : [0, (0.5 * Δx)...]

    if ages[1] == 0
        m0 = mx[1]
        a0 = a0(m0, sex, mode)
        ax[1] = a0
    end


    if ages[n] == 110 || !single_ages
        aend = 1 / mx[n]
        ax[n] = aend
    end

    if single_ages
        qx = mx ./ (1 .+ (1 .- ax) .* mx)
    else
        mx′ = mx[1:end-1]
        ax′ = ax[1:end-1]
        qx = (Δx .* mx′) ./ (1 .+ (Δx .- ax′) .* mx′)
        qx = [qx..., 1]
    end

    if ages[n] == 110 || !single_ages
        qx[n] = 1
    end

    px = 1 .- qx

    lx = cumprod(cat([1], px[1:end-1], dims=1)) .* radix
    dx = lx .* qx

    if single_ages
        Lx = lx .- (1 .- ax) .* dx
        Tx = reverse(cumsum(reverse(Lx)))
    else
        lx′ = lx[1:end-1]
        dx′ = dx[1:end-1]
        ax′ = ax[1:end-1]
        Lx = Δx .* (lx′ .- dx′) .+ ax′ .* dx′
        Lx = [Lx..., lx[end] / mx[end]]
        Tx = reverse(cumsum(reverse(Lx)))
    end

    ex = Tx ./ lx

    if debug
        df = DataFrame(:Age => ages, :mx => mx, :ax => ax, :qx => qx, :ex => ex)

        println(df)
    end
    return ex[indexin(at_age, ages)]

end


function lexpectancies(
    mxt::Matrix{Float64},
    xages,
    tyears;
    gender,
    at_age=[0],
    mode::CalculationChoice=CC_JULIA,
    radix=1_000_000,
    debug=false
)

    output = Matrix{Float64}(undef, length(at_age), length(tyears))

    for i in eachindex(tyears)
        year = tyears[i]
        mx = mxt[:, i]

        les = expected_lifetime(mx, xages; sex=gender, mode=mode, at_age=at_age, radix=radix, debug=debug)
        output[:, i] = les
    end

    return output

end