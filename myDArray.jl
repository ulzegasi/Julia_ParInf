## chk_indexes is an array of tuples defining the lengths of the chunks on each process
## ex:
## Out [46]: 3-element Array{Any,1}:
## (1:9,)
## (10:14,)
## (15:23,)
##

function myDArray(init, dims, procs, dist, chk_indexes)
    np      = prod(dist)
    procs   = procs[1:np]
    cuts    = [[chk_indexes[i][1][1] for i = 1:num_workers] for ii = 1:1]
    idxs    = chk_indexes
    chunks  = Array(RemoteRef, dist...)
    for i = 1:np
        chunks[i] = remotecall(procs[i], init, idxs[i])
    end
    p = 0
    for i = 1:length(procs)
        if procs[i] == myid()
            p = i
        end
    end
    p = max(1 , p)
    A = remotecall_fetch(procs[p], r->typeof(fetch(r)), chunks[p])
    DArray{eltype(A),length(dims),A}(dims, chunks, procs, idxs, cuts)
end