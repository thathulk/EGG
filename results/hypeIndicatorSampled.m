% HYPEINDICATORSAMPLED calculates the HypE fitness
%   HYPEINDICATORSAMPLED( POINTS, BOUNDS, K, NROFSAMPLES  ) 
%   samples the HypE fitness values for objective vectors POINTS. 
%
%   POINTS:      objective vectors as rows, all to be minimized
%   BOUNDS:      reference point
%   K:           parameter of HypE
%   NROFSAMPLES: nr of samples to be used
%
%   Example: f =  hypeIndicatorSampled( [1 3; 3 1], [4 4], 1, 10000 )

function F = hypeIndicatorSampled( points, bounds, nrOfSamples)

    Atemp = points;
    %remove the dominated solution
    for i=1:size(Atemp,2)
        btemp=repmat(bounds,size(Atemp,1),1);
        Atemp = Atemp(Atemp(:,i)<btemp(:,i),:);
    end
    Atemp=unique(Atemp,'rows');
    points = Atemp;
    clear Atemp btemp
    if(size(unique(points,'rows'),1)>1)
    [nrP, dim] = size(points);
    F = zeros(1,nrP);
    k=nrP;
    alpha = zeros(1,nrP);
    for i = 1 : k
        j = 1:i-1;
        alpha(i) = prod( (k-j) ./ (nrP - j ) )./i;
    end  
    
    
    if( length(bounds) == 1 )
        bounds = repmat( bounds,1, dim );
    end
    
    BoxL = min(points);
    
    S = rand(nrOfSamples,dim)*diag( bounds - BoxL) ...
    + ones( nrOfSamples,dim)*diag(BoxL);

    dominated = zeros( nrOfSamples, 1 );
    for j = 1 : nrP
        B = S - repmat(points(j,:),nrOfSamples,1);
        ind = find( sum( B >= 0, 2) == dim);
        dominated(ind) = dominated(ind) + 1;
    end
    
    for j = 1 : nrP
        B = S - repmat(points(j,:),nrOfSamples,1);
        ind = find( sum( B >= 0, 2) == dim);
        x = dominated(ind);        
        F(j) = sum(  alpha(x)  );
    end    
    F = F'*prod( bounds - BoxL)/nrOfSamples;
    else
        F=0.0;
    end
    

    