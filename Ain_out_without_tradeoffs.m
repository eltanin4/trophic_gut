function [m2b, b2m] = Ain_out(b_real, i_all_filt, j_all_filt, v_all_filt)
% This part is devoted for calculating metabolite consumption matrix m2b and
% byproduct production matrix b2m. Here m2b is determined by the real bacterial 
% abundance from an individual whose metagenomic bacterial abundance is b_real.
% i_all_filt, j_all_filt, v_all_filt is a triplet of the network edges.
    %i1=find((v_all_filt==2)+(v_all_filt==5)+(v_all_filt==6));
    i1=find((v_all_filt==2)+(v_all_filt==5)); % find indexes of all metabolite uptake edges
    m2b=ones(size(i1));    % m2b is the metabolite uptake rates
                            % and it is treated as 1 for all metabolites.
    %m2b=rand(size(i1));   % randomizing metabolite uptakes
    %m2b=exp(6.*randn(size(i1)));  % another version of randomization
    
    m2b=sparse(i_all_filt(i1),j_all_filt(i1),m2b,2244,2244);
    in_degree=sum(m2b,1)';  % calculate in degrees of all bacterial species
    [i2,j2,v2]=find(m2b);

    %m2b=sparse(i2,j2,v2./(in_degree(j2)),2244,2244);  % trade-offs

    %i1=find((v_all_filt==3)+(v_all_filt==5)+(v_all_filt==8)); 
    i1=find((v_all_filt==3)+(v_all_filt==5));  % find indexes of all metabolite produced
    b2m=ones(size(i1));
    b2m=sparse(i_all_filt(i1),j_all_filt(i1),b2m,2244,2244);
    out_degree=full(sum(b2m, 1)'); % calculate out degrees of all bacterial species
    [i2,j2,v2]=find(b2m);
    b2m=sparse(i2,j2,v2./out_degree(j2),2244,2244);  % byproducts equally distributed

    [i3,j3,v3]=find(m2b);
    m2b=sparse(i3,j3,v3.*b_real(j3),2244,2244);  % metabolite uptake is weighted 
                                                 % by bacterial abundance
    kc=sum(m2b,2);         % normaliztion factor for split of each metabolite
    [i2,j2,v2]=find(m2b);
    m2b=sparse(i2,j2,v2./kc(i2),2244,2244);  % metabolite splitting is normalized 
end