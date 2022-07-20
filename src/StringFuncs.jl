"""
get_particle_string(PART_TYPE, p::Vector{<:Real})

Description of get_particle_string
--------------------------------
Returns one-line string for particle in form "TYPE pₓ p𝔂 p𝐳 0"

example:

"electron px = 1.0 py = 1.0 pz = 1.0 dt = 0"

"3 1.0 1.0 1.0 0"
"""
function get_particle_string(PART_TYPE, p::Vector{<:Real})

    return string(PART_TYPE, " ", p[2], " ", p[3], " ", p[4], " 0")
end

function get_particle_string(PART_TYPE, p::DataFrameRow)

    return string(PART_TYPE, " ", p.px, " ", p.py, " ", p.pz, " 0")
end

"""
    get_header_string(ID)

Description of get_header_string
------------------------------
Returns one-line string for header in form "ID 0 2" (2 stands for number of particles)

example:

"ID = 0 t = 0 particles = 2"

"0 0 2"
"""
function get_header_string(ID)
    return string(ID, " 0 2")
end


"""
```get_event_string(ID, PART_TYPE, p1::Vector{<:Real}, p2::Vector{<:Real})```
```get_event_string(ID, PART_TYPE, p1::DataFrameRow, p2::DataFrameRow)```

Description of get_event_string
------------------------------
Returns one event, separated at the end by #.
"""
function get_event_string(ID, PART_TYPE, p1::Vector{<:Real}, p2::Vector{<:Real})
    sHeader = get_header_string(ID)
    sPart1 = get_particle_string(PART_TYPE, p1)
    sPart2 = get_particle_string(PART_TYPE, p2)

    return sHeader * "\n" * sPart1 * "\n" * sPart2 * "\n#\n"
end

function get_event_string(ID, PART_TYPE, p1::DataFrameRow, p2::DataFrameRow)
    sHeader = get_header_string(ID)
    sPart1 = get_particle_string(PART_TYPE, p1)
    sPart2 = get_particle_string(PART_TYPE, p2)

    return sHeader * "\n" * sPart1 * "\n" * sPart2 * "\n#\n"
end
