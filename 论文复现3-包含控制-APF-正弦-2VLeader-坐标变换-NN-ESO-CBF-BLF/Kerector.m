function C = Kerector( A, B )
%KERECTOR 此处显示有关此函数的摘要
%   计算克罗尼可积
col_A = length(A(1,:));
row_A = length(A(:,1));
col_B = length(B(1,:));
row_B = length(B(:,1));
k_col_C = 1;
k_row_C = 1;
for i=1:row_A
    for j=1:col_A
        C(k_row_C:k_row_C+row_B-1,k_col_C:k_col_C+col_B-1)=A(i,j)*B;
        k_col_C=k_col_C+col_B;
    end
    k_col_C = 1;
    k_row_C=k_row_C+row_B;
    
end

end

