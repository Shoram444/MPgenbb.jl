### MPGenbb generator is a module used for generating custom .genbb files, which serve as an input for custom generated events into Falaise simulation module of SuperNEMO collaboration. 
<br>

#### The following example shows how to use the module to simply generate a file with momenta of each particle generated using a $2\nu\beta\beta$ decay spectrum provided by R. Dvornicky rebinned using MPRebinSpectra module (<a href="https://github.com/Shoram444/MPRebinSpectra.jl">MPRebinSpectra.jl</a>)


##### The structure of .genbb files is as follows: 
+ Each event is separated by ``#``
+ Each event consists of:
    1. Header line which has fields: ``ID Time nP`` (where ID is event ID from 0 to N; Time is starting time; nP is the number of particles within the event)
    2. Particle line with fields: ``Type px py pz dt`` (where particle ``Type`` is 1- gamma, 2- positron, 3-electron; ``px py pz`` are momentum components; ``dt`` is time difference)

An example event - 2 positrons, starting at ``Time = 0``, with ``dt = 0`` and momenta ``px1 = py1 = pz1 = 1, px2 = py2 = pz2 = 1`` - is generated as:
<br>
``0 0 2 `` <br>
``2 1 1 1 0`` <br>
``2 2 2 2 0`` <br>
``#``

#### Docstrings for individual functions are found in ``/src``. Or by using ``?`` in julia REPL. 



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




<div class="data-frame"><p>614,474 rows × 8 columns</p><table class="data-frame"><thead><tr><th></th><th>E1</th><th>minE</th><th>maxE</th><th>minG</th><th>maxG</th><th>a</th><th>b</th><th>cdf</th></tr><tr><th></th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>0.00093</td><td>0.0</td><td>0.00093</td><td>0.0</td><td>0.402019</td><td>432.278</td><td>0.0</td><td>1.86939e-7</td></tr><tr><th>2</th><td>0.00093</td><td>0.00093</td><td>0.01493</td><td>0.402019</td><td>0.412043</td><td>0.716006</td><td>0.401353</td><td>5.88537e-6</td></tr><tr><th>3</th><td>0.00093</td><td>0.01493</td><td>0.02693</td><td>0.412043</td><td>0.421358</td><td>0.776314</td><td>0.400452</td><td>1.08858e-5</td></tr><tr><th>4</th><td>0.00093</td><td>0.02693</td><td>0.04093</td><td>0.421358</td><td>0.433615</td><td>0.875475</td><td>0.397782</td><td>1.68706e-5</td></tr><tr><th>5</th><td>0.00093</td><td>0.04093</td><td>0.05893</td><td>0.433615</td><td>0.450447</td><td>0.93513</td><td>0.39534</td><td>2.48272e-5</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>



### Setting up initial conditions .


```julia
const MASS        = 0.511 # particle mass in MeV/c²
const PART_TYPE   = 3     # particle type: 1- gamma, 2-positron, 3-electron
const OUTFILE     = "input_module.genbb" # name of the output file
const PROCESS     = "82Se 0νββ - Energy generated via SpectrumG0 - from R. Dvornicky";
```


```julia
NEVENTS          = 100  # how many events should be created
```




    100



#### To save the file use a simple for loop, at each step generate electron energies ```T``` and momentum vectors ``p1`` and ``p2`` (where direction of ``p1`` is random unifrom in sphere and direction of ``p2`` is correlated to ``p1``). 
#### Then generate ``get_event_string``. And save.


```julia
open(OUTFILE, "w") do file
    Threads.@threads for id in 0:NEVENTS-1
        T  = sample_energies(df)

        p1 = get_first_vector(T[1], MASS)
        p2 = get_second_vector(T[2], MASS, p1)

        if id%10_000 == 0 && id >1
            println("generated $id events!")
        end

        write(file, get_event_string(id, PART_TYPE , p1, p2))
    end
end
```

$$\begin{bmatrix} 1 & 2 & 1 \\ 3 & 0 & 1 \\ 0 & 2 & 4 \end{bmatrix}$$
