export getcommonedge

# returns two integers (i1,i2). For e.g. i1 the return value is
# +/-{1,2,3} depending whether the edge common to cell1 and cell2
# is local edge 1,2,3 in cell1 and whether it is along or against
# the internal orientation of cell1. Similar for i2.
function getcommonedge(cell1, cell2)

    cv = intersect(cell1,cell2)
    @assert length(cv) == 2

    I1 = findfirst(cell1, cv)
    I2 = findfirst(cell2, cv)

    mod1 = (i,n) -> (i-1)%n + 1

    if I1[2] == I1[1]+1
        i1 = mod1(I1[1]+2,3)
    else
        i1 = mod1(I1[1]+1,3)
    end

    if I2[2] == I2[1]+1
        i2 = mod1(I2[1]+2,3)
    else
        i2 = mod1(I2[1]+1,3)
    end

    return i1, i2
end