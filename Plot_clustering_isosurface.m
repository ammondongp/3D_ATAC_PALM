function [Fv, Vol, Vol_nucleus, Dens, Dens_nucleus] = Plot_clustering_isosurface(X, IDX)

    k=max(IDX);

    Fv=repmat({struct('vertices',[0 0 0], 'faces',[0 0 0])}, 1, k+1);
    Vol=zeros(1,k);
    Dens=zeros(1,k);
    
    for i=1:k
        Xi=X(IDX==i,:);

        [Fv{i+1}.faces, Vol(i)] = boundary(Xi);
        Fv{i+1}.vertices = Xi;
        
        Dens(i)=double(size(Xi,1)/Vol(i));

    end
    
        [Fv{1}.faces, Vol_nucleus] = convhull(X);
        Fv{1}.vertices = X;
        
        Dens_nucleus=double(size(X,1)/Vol_nucleus);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Colors=hsv(k);
    
    figure;
        p=patch(Fv{1});
        p.FaceColor = [0 0 0];
        p.EdgeColor = 'none';
        hold on;
    
    for j=1:k
        p=patch(Fv{j+1});
    %     isonormals(image_data,p);
        p.FaceColor = 'green';% Green color instead of  Colors(j,:);
        p.EdgeColor = 'none';
        
        hold on;
    end
    
    alpha(0.2);
    
        daspect([1 1 1]);
        view(-45,45); 
        axis tight;
        camlight;
        lighting gouraud;
        hold off;
    
end