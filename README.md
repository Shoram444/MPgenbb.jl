# MPGenbb

[![Build Status](https://travis-ci.com/Shoram444/MPGenbb.jl.svg?branch=main)](https://travis-ci.com/Shoram444/MPGenbb.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/Shoram444/MPGenbb.jl?svg=true)](https://ci.appveyor.com/project/Shoram444/MPGenbb-jl)
[![Coverage](https://codecov.io/gh/Shoram444/MPGenbb.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Shoram444/MPGenbb.jl)
[![Coverage](https://coveralls.io/repos/github/Shoram444/MPGenbb.jl/badge.svg?branch=main)](https://coveralls.io/github/Shoram444/MPGenbb.jl?branch=main)


### MPGenbb generator is a module used for generating custom .genbb files, which serve as an input for custom generated events into Falaise simulation module of SuperNEMO collaboration. 
<br>

#### The following example shows how to use the module to simply generate a file with momenta of each particle generated using a $2\nu\beta\beta$ decay spectrum provided by R. Dvornicky rebinned using MPRebinSpectra module (<a href="https://github.com/Shoram444/MPRebinSpectra.jl">MPRebinSpectra.jl</a>)


##### The structure of .genbb files is as follows: 
+ Each event is separated by ``#``
+ Each event consists of:
    1. Header line which has fields: ``ID Time nP`` (where ID is event ID from 0 to N; Time is starting time; nP is the number of particles within the event)
    2. Particle line with fields: ``Type px py pz dt`` (where particle ``Type`` is 1- gamma, 2- positron, 3-electron; ``px py pz`` are relativistic momentum components; ``dt`` is time difference)

An example event - 2 positrons, starting at ``Time`` = 0, with ``dt`` = 0 and momenta ``px1 = py1 = pz1 = 1, px2 = py2 = pz2 = 1`` - is generated as:
<br>
``
0 0 2 
2 1 1 1 0
2 2 2 2 0
#
``


```julia
using MPGenbb, DataFrames, CSV
Base.displaysize() = (5, 100) # this line sets display of rows 5 and cols 100
```

### The input spectrum must be rebinned with MPRebinSpectra module and saved as a ``DataFrame`` type.


```julia
inFile = string("spectrumG0_Rebinned_prec0001.csv")
```




    "spectrumG0_Rebinned_prec0001.csv"




```julia
df = CSV.File(inFile) |> DataFrame
```




<div class="data-frame"><p>608,481 rows × 8 columns</p><table class="data-frame"><thead><tr><th></th><th>E1</th><th>minE</th><th>maxE</th><th>minG</th><th>maxG</th><th>a</th><th>b</th><th>cdf</th></tr><tr><th></th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>0.00093</td><td>0.00093</td><td>0.01493</td><td>0.402172</td><td>0.4122</td><td>0.71628</td><td>0.401506</td><td>5.70061e-6</td></tr><tr><th>2</th><td>0.00093</td><td>0.01493</td><td>0.02693</td><td>0.4122</td><td>0.42152</td><td>0.776611</td><td>0.400606</td><td>1.07029e-5</td></tr><tr><th>3</th><td>0.00093</td><td>0.02693</td><td>0.04093</td><td>0.42152</td><td>0.433781</td><td>0.87581</td><td>0.397934</td><td>1.669e-5</td></tr><tr><th>4</th><td>0.00093</td><td>0.04093</td><td>0.05893</td><td>0.433781</td><td>0.45062</td><td>0.935488</td><td>0.395491</td><td>2.46496e-5</td></tr><tr><th>5</th><td>0.00093</td><td>0.05893</td><td>0.10993</td><td>0.45062</td><td>0.498454</td><td>0.937931</td><td>0.395348</td><td>4.8851e-5</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>



### Initial conditions need to be set. 


```julia
const MASS        = 0.511 # particle mass in MeV/c²
const PART_TYPE   = 3     # particle type: 1- gamma, 2-positron, 3-electron
const FILENAME    = "input_module.genbb" # name of the output file
const PROCESS     = "82Se 0νββ - Energy generated via SpectrumG0 - from R. Dvornicky";
```


```julia
NEVENTS          = 100  # how many events should be created
```




    100




```julia
open(FILENAME, "w") do file
    for id in 0:NEVENTS
        T::Tuple{Float64, Float64} = sample_energies(df)
        ϕ::Tuple{Float64, Float64} = sample_phi()
        θ::Tuple{Float64, Float64} = sample_theta()
        
        if id%10_000 == 0 && id >1
            @show "generated $id events!"
        end

        write(file, get_event_string(id, PART_TYPE , T, MASS, ϕ, θ))
    end
end
```
