
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
Returns θdif in radians distrtibuted according to pdf = 0.5 - 0.5x (where x = cos(θdif)). θdif is the angle between the two electrons.

Sampling equation is:
    + θdif (the angle between the two electrons) is sampled from pdf = 0.5 - 0.5x. This pdf represents the 1 - cos(θ) angular distribution. 
        It is sampled through inverse CDF method. Integrating pdf and finding the inverse function we get ``` x = 1 - 2\\sqrt{1 - \\gamma}```
    + returned angle is the acos(x)

"""
function sample_theta_dif()
    a = 1. 
    b = -2.
    c = -3 + 4*rand(Uniform())
    cosθdif = solve_quadratic(a,b,c)[1]
    return acos(cosθdif)
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
    
    # while gamma > maximum(df.cdf)  # sanity check that gamma isnt some crazy number
    #     @show "incorrectly sampled gamma"
    #     @show gamma
    #     @show minimum(df.cdf), maximum(df.cdf)

    #     gamma = rand(Uniform(0.0, maximum(df.cdf)))
    #     # sample_energies(df)
    # end

    
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

function solve_quadratic(a, b, c)
    d  = sqrt(b^2 - 4*a*c)
    return (-b - d) / (2*a), (-b + d) / (2*a)
end

