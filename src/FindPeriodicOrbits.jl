function find_approx_periodic_orbits(mapping, period, search_grid, accuracy)
    idxs = findall(<(accuracy), abs.((mapping^period).(search_grid.elems) - search_grid.elems))
    search_grid.marks[idxs] .= true
    return search_grid
end