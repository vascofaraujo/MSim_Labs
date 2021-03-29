function [index_l] = findClosest(BPM_matrix, BPMtoFind, index_m, size)

oldDiference = abs(BPM_matrix(1,1) - BPMtoFind);
index_l = 1;
    for j = 1:size
        currentDiference = abs(BPM_matrix(index_m,j) - BPMtoFind);
        if currentDiference < oldDiference
            oldDiference = currentDiference;
            index_l = j;         
    end
    



end