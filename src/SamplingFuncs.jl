
"""
    sample_phi()

Description of sample_phi
------------------------------
Returns (ϕ₁, ϕ₂) tuple in radians Uniformly distributed. Phi is the polar angle. See Distributions.Uniform

"""
function sample_phi()
    return acos.(1 .- 2 .* rand(Uniform(),2))
end

"""
    sample_theta()

Description of sample_theta
------------------------------
Returns (θ₁, θ₂) tuple in radians Uniformly distributed. See Distributions.Uniform

Sampling equation is:
``` acos(1 - 2*U) ; where U is uniform random number; U ∈ (0,1). ```

"""
function sample_theta()
    γ = rand(Uniform())
    
    θ1 = 2*π*rand(Uniform())

    a = 1. #fix correct values here and check what x1,x2 give as values! 
    b = 2.
    c = 3 - 4*γ
    @show x1,x2 = sample_quadratic(a,b,c)

    

    return acos.(1 .-2*rand(Distributions.Uniform(),2))
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
    gamma = rand(Uniform()) # random uniform number from 0 to 1
    
    if gamma > maximum(df.cdf) || gamma < minimum(df.cdf) # sanity check that gamma isnt some crazy number
        @show "sanity check"
        sample_energies(df)
        
    end
    
    r_id = findfirst(x -> x .>= gamma,df.cdf)
        
    r_id > 1 ? Pᵢ₋₁ = df[r_id-1, 8] : Pᵢ₋₁ = 0 
    
    df_Pᵢ= df[r_id, :]    
    
    a = 0.5*df_Pᵢ.a * thickness
    b = df_Pᵢ.b * thickness
    c = Pᵢ₋₁ - gamma - a*df_Pᵢ.minE^2 - b*df_Pᵢ.minE
    
    x1, x2 = solve_quadratic(a, b, c)
    
    if (df_Pᵢ.minE < x1) && (df_Pᵢ.maxE) > x1
        return rand(Uniform(df_Pᵢ.E1, df_Pᵢ.E1+thickness)), x1  #E1 is sampled as uniform number between steps, 
                                                                #E2 is sampled by CDF
    else
        return rand(Uniform(df_Pᵢ.E1, df_Pᵢ.E1+thickness)), x2 
    end
    
end

function solve_quadratic(a, b, c)
    d  = sqrt(b^2 - 4*a*c)
    return (-b - d) / (2*a), (-b + d) / (2*a)
end
