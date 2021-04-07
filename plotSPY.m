
formatSpec = '%d,%d,%f';
sizeA = [3 Inf];
file1 = fopen('mathatchel_342_alfa0.txt');
dados1 = fscanf(file1, formatSpec, sizeA);
dados1 = dados1';

sizeI = max(dados1(:,1));
sizeJ = max(dados1(:,2));

A_sparse1 = zeros(sizeI,sizeJ);

ndados = size(dados1,1);
for cnt = 1:ndados
    indexI = dados1(cnt,1) + 1; 
    indexJ = dados1(cnt,2) + 1;
    valor = dados1(cnt,3);
    
    A_sparse1(indexI, indexJ) = valor;
end

gain = A_sparse1'*A_sparse1;

%%
%Ordernacao da matriz ganho
amd_P = amd(A_sparse1);
colpermP = colperm(A_sparse1);
symrcmP = symrcm(A_sparse1);
%dissectP = dissect(gain);


gain_amd = sparse(A_sparse1(amd_P, amd_P));
gain_colp = sparse(A_sparse1(colpermP, colpermP));
gain_sym = sparse(A_sparse1(symrcmP,symrcmP));
%gain_dis = sparse(gain(dissectP,dissectP));

[q, ramd] = qr(gain_amd);
[q2, rcolp] = qr(gain_colp);
[q3, rsym] = qr(gain_sym);

[q_no_order, r_no_order] = qr(A_sparse1);

%%
figure
set(gca, 'FontSize', 20);
subplot(1,2,1)
spy(sparse(r_no_order))
title('R Factor - No ordering', 'FontSize', 20);
subplot(1,2,2)
spy(sparse(ramd))
title('R Factor - AMD', 'FontSize', 20);

% %%
% figure
% set(gca, 'FontSize', 20);
% spy(gain_amd)
% title('AMD Ordered Matrix', 'FontSize', 20);
%%
% figure
% set(gca, 'FontSize', 20);
% spy(A_sparse1)
% title('Matrix without ordering', 'FontSize', 20);
