function spyKKT(name, ns, mpc)

    %name = 'mat-ipopt_000-02.iajaa';
    %ns = 9;
    %mpc = case118;

    addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

    A = readcsr(name, 0, 1);

    [P Pinv, npart] = KKTpermute(mpc, ns);
    AP = A(P,P'); 
    if (A - AP(Pinv, Pinv') ~= sparse(size(A,1),size(A,2)))
       error('Inverse permutation does not result in original matrix.') 
    end
    
    spy(A);
    figure; spy(AP);
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'ytick',[])
    %set(gca,'FontSize',20)
    
    %visualize magnitude using surf()
%   [X,Y] = meshgrid(1:size(A,1),1:size(A,2));
%   figure; surf(X,Y,A)
    
    %visualize magnitude using patch()
%   [m, n] = size(A);
%   nonzeroInd = find(A);
%   [x, y] = ind2sub([m n], nonzeroInd);
% 
%   figure();
%   hp = patch(x, y, A(nonzeroInd), ...
%            'Marker', 's', 'MarkerFaceColor', 'flat', 'MarkerSize', 4, ...
%            'EdgeColor', 'none', 'FaceColor', 'none');
%   set(gca, 'XLim', [0, n + 1], 'YLim', [0, m + 1], 'YDir', 'reverse', ...
%   'PlotBoxAspectRatio', [n + 1, m + 1, 1]);
% 
%   colorbar();

end