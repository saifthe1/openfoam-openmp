    piso = new dictionary(mesh.solutionDict().subDict("PISO"));

    nCorr = new int(readInt(piso->lookup("nCorrectors")));

    nNonOrthCorr = new int(piso->lookupOrDefault<int>("nNonOrthogonalCorrectors", 0));

    momentumPredictor = new bool(piso->lookupOrDefault<Switch>("momentumPredictor", true));

    transonic = new bool(piso->lookupOrDefault<Switch>("transonic", false));

    nOuterCorr = new int(piso->lookupOrDefault<int>("nOuterCorrectors", 1));
