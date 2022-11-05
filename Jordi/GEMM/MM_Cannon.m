function C = MM_Cannon(A, B, p)
    A_nrows = length(A);
    m = A_nrows / p;
    

    A_cell = cell(p);
    B_cell = cell(p);
    for i = 1:p
        for j = 1:p
            ii = (i-1)*m + (1:m);
            jj = (j-1)*m + (1:m);
            A_cell{i,j} = A(ii,jj);
            B_cell{i,j} = B(ii,jj);
        end
    end
    A = A_cell;
    B = B_cell;
    
    c = cell(p);
    a = cell(p);
    b = cell(p);
    A_comm = cell(p);
    B_comm = cell(p);
    
    for i = 1:p
        for j = 1:p
            k = mod(i+j, p)+1;
            a{i,j} = A{i,k};
            b{i,j} = B{k,j};
            c{i,j} = 0;
        end
    end
    
    for l = 1:p
        for i = 0:(p-1)
            for j = 0:(p-1)
                c{i+1,j+1} = c{i+1,j+1} + a{i+1,j+1} * b{i+1,j+1};
                
                A_comm{i+1                , mod(j + p - 1, p)+1} = a{i+1,j+1};
                B_comm{mod(i + p - 1, p)+1, j+1                } = b{i+1,j+1};
            end
        end
        
        for i = 1:p
            for j = 1:p
                a{i,j} = A_comm{i,j};
                b{i,j} = B_comm{i,j};
            end
        end
    end
    

%     c = cell(p);
%     for i = 1:p
%         for j = 1:p
%             c{i,j} = 0;
%             for k = 1:p
%                 c{i,j} = c{i,j} + A{i,k}*B{k,j};
%             end
%         end
%     end

    C = zeros(A_nrows);
    for i = 1:p
        for j = 1:p
            ii = (i-1)*m + (1:m);
            jj = (j-1)*m + (1:m);
            C(ii,jj) = c{i,j};
        end
    end
end