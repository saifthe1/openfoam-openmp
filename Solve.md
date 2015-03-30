# Solve ! #
src » finiteVolume » fvMatrices » fvMatrix
```
template<class Type>
Foam::lduMatrix::solverPerformance Foam::fvMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::solve(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    lduMatrix::solverPerformance solverPerfVec
    (
        "fvMatrix<Type>::solve",
        psi_.name()
    );

    scalarField saveDiag = diag();

    Field<Type> source = source_;

    // At this point include the boundary source from the coupled boundaries.
    // This is corrected for the implict part by updateMatrixInterfaces within
    // the component loop.
    addBoundarySource(source);

    typename Type::labelType validComponents
    (
        pow
        (
            psi_.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    for(direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        // copy field and source

        scalarField psiCmpt = psi_.internalField().component(cmpt);
        addBoundaryDiag(diag(), cmpt);

        scalarField sourceCmpt = source.component(cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        FieldField<Field, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(cmpt)
        );

        lduInterfaceFieldPtrsList interfaces =
            psi_.boundaryField().interfaces();

        // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
        // bouCoeffsCmpt for the explicit part of the coupled boundary
        // conditions
        initMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        updateMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        lduMatrix::solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        solverPerf.print();

        if
        (
            solverPerf.initialResidual() > solverPerfVec.initialResidual()
         && !solverPerf.singular()
        )
        {
            solverPerfVec = solverPerf;
        }

        psi_.internalField().replace(cmpt, psiCmpt);
        diag() = saveDiag;
    }

    psi_.correctBoundaryConditions();

    return solverPerfVec;
}
```