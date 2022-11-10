
"""
sample_theta()

Description of sample_theta
------------------------------
Returns θ in radians Uniformly distributed. θ is the polar angle. See Distributions.Uniform

"""
function sample_theta()
    return acos(1 - 2 * rand(Uniform())) ## polar angle in radians
end

"""
    sample_phi()

Description of sample_phi
------------------------------
Returns ϕ in radians Uniformly distributed. phi is the azimuthal angle. See Distributions.Uniform

Sampling equation is:
    + First θ1 is sampled uniformly from 2π uniform pdf.
    + Second angle θdif (the angle between the two electrons) is sampled from pdf = 0.5 - 0.5x. This pdf represents the 1 - cos(θ) angular distribution. 
        It is sampled through inverse CDF method. Integrating pdf and finding the inverse function we get ``` x = 1 - 2\\sqrt{1 - \\gamma}```
    + θ2 is determined as θ2 = θdif + θ1.

"""
function sample_phi()
    return 2*π*rand(Uniform()) ## azimuthal angle in radians
end


"""
    sample_theta_dif()

Description of sample_theta_dif
------------------------------
Returns θdif in radians distrtibuted according to pdf = 0.5 - 0.5*k*x (where x = cos(θdif)). θdif is the angle between the two electrons.

Sampling equation is:
    + θdif (the angle between the two electrons) is sampled from pdf = 0.5 - 0.5x. This pdf represents the 1 - cos(θ) angular distribution. 
        It is sampled through inverse CDF method. Integrating pdf and finding the inverse function we get ``` x = 1 - 2\\sqrt{1 - \\gamma}```
    + returned angle is the acos(x)

If no k is specified, it is assumed to be k = -1.0. 
"""
function sample_theta_dif()
    a = _k 
    b = 2.
    c = 2-_k -4*rand(Uniform())
    cosθdif = solve_quadratic(a,b,c)[1]
    return acos(cosθdif)
end

function sample_theta_dif(_k::Real)
    
    if (_k == 0.0)
        return 0.5
    end
    
    a = _k/4 
    b = 1/2
    c = 1/2 - _k/4 - rand(Uniform())
    θ = solve_quadratic(a,b,c)
    
    if( -1.0 <= θ[1] <= 1.0 )
        return acos(θ[1])
    elseif ( -1.0 <= θ[2] <= 1.0 )
        return acos(θ[2])
    else
        @show "soulution is outside of range u ∈ (-1.0, 1.0)"
    end
    return -100.0
end


function sample_theta_dif_uniform()
    return rand(Uniform(0, π))
end



"""
#### function ```sample_energies(df::DataFrame, thickness = 0.001)```
<br>

    Description of ```sample_energies```
    ------------------------------
Returns a tuple of energies (E1, E2) given by the distribution specified in dataframe ```df```. 

```df``` must have fields: ```:E1, :minE, :maxE, :a, :b, :cdf (given by get_cdf!(df))```

Sampling is done in the following way. 
E2 is sampled using inverse CDF method - where cdf is obtained from ```get_cdf!(df)```. 
E1 is dependent on E2 as it can only be sampled in the given row where E2 was chosen.
E1 is sampled as a random uniform number between (row of E1, next row).

<br>
It takes roughly 100s to sample 1e5 energies. 

"""
function sample_energies(df::DataFrame, thickness = 0.001)
    gamma = rand(Uniform(0.0, maximum(df.cdf))) # random uniform number from 0 to maximum of cdf column (it is not exactly 1.0 due to discrete numbers used in PC)
    
    r_id = findfirst(x -> x .>= gamma,df.cdf)
        
    r_id > 1 ? Pᵢ₋₁ = df[r_id-1, 8] : Pᵢ₋₁ = 0 
    
    df_Pᵢ= df[r_id, :]    
    
    a = 0.5* df_Pᵢ.a * thickness
    b = df_Pᵢ.b * thickness
    c = Pᵢ₋₁ - gamma - a*df_Pᵢ.minE^2 - b*df_Pᵢ.minE
    
    x1, x2 = solve_quadratic(a, b, c)
    
    if (df_Pᵢ.minE < x1) && (df_Pᵢ.maxE > x1)
        return rand(Uniform(df_Pᵢ.E1, df_Pᵢ.E1+thickness)), x1  #E1 is sampled as uniform number between steps, 
                                                                #E2 is sampled by CDF
    elseif (df_Pᵢ.minE < x2) && (df_Pᵢ.maxE > x2)
        return rand(Uniform(df_Pᵢ.E1, df_Pᵢ.E1+thickness)), x2 
    else
        sample_energies(df)
    end
    
end

function sample_energies_uniform(_max)
    T  = (rand(Uniform(0,_max)), rand(Uniform(0,_max))) 
    while( sum(T) >=  _max)
        T  = (rand(Uniform(0,_max)), rand(Uniform(0,_max))) 
    end
    return T
end


function solve_quadratic(a, b, c)
    d  = sqrt(b^2 - 4*a*c)
    return (-b - d) / (2*a), (-b + d) / (2*a)
end

