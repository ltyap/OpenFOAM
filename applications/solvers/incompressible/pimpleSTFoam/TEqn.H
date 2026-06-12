{
  alphat = turbulence->nut()/Prt;
  alphat.correctBoundaryConditions();
  volScalarField alphaEff("alphaEff", alpha + alphat);

  fvScalarMatrix TEqn	    
  (
      fvm::ddt(T) 
    + fvm::div(phi, T)
    - fvm::laplacian(alphaEff, T)
    == 
      fvModels.source(T)
  );

  TEqn.relax();
  fvConstraints.constrain(TEqn);
  TEqn.solve();
  T.correctBoundaryConditions();

  
}
