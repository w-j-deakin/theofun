function POL = ParetoLandscape(lsin, optin)

    varin = zeros(length(lsin(1).zvec),length(lsin));

    for i = 1:length(lsin)
        varin(:,i) = lsin(i).zvec;
    end

    P = ParetoND(varin, optin);

    R = P.RankRatio();

    rGrid = lsin(1).grid(R);

    POL = landscape(lsin(1).xgrid, lsin(1).ygrid, rGrid);

end