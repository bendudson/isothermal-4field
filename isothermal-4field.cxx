
#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>
//#include <invert_laplace.hxx>
#include <bout/invert/laplacexz.hxx>      // Laplacian inversion (BSTING only has XZ-PETSc)

class MergingFlux : public PhysicsModel {
protected:
  int init(bool restarting) {

    // Read options
    Options *opt = Options::getRoot()->getSection("model");
    // mesh->periodicZ=true;
    
    /////////////////////////////////////////////////////////////////////
    // Normalisations
    
    OPTION(opt, Nnorm, 1e19);  // Reference density [m^-3]
    OPTION(opt, Tnorm, 100);   // Reference temperature [eV]
    OPTION(opt, Bnorm, 1.0);   // Reference magnetic field [T]
    OPTION(opt, AA, 2.0); // Atomic mass number

    BoutReal mi = SI::Mp * AA; // Ion mass [kg]

    mi_me = mi / SI::Me; // Ion mass to electron mass ratio
    
    SAVE_ONCE4(Nnorm, Tnorm, Bnorm, mi_me);
    
    // Normalisation frequency
    Omega_ci = SI::qe*Bnorm / mi;
    Cs = sqrt(SI::qe * Tnorm / mi);
    rho_s = Cs / Omega_ci;
    beta = SI::mu0 * SI::qe * Nnorm * Tnorm / SQ(Bnorm);

    output.write("\tOmega_ci = %e, Cs = %e, rho_s = %e, beta = %e\n", Omega_ci, Cs, rho_s, beta);
    output.write("\tmi_me = %e\n", mi_me);
    
    SAVE_ONCE4(Omega_ci, Cs, rho_s, beta);
    
    /////////////////////////////////////////////////////////////////////
    // Read and normalise the temperatures
    
    OPTION(opt, Te, 100.);   // Electron temperature [eV]
    OPTION(opt, Ti, 100.);   // Ion temperature [eV]
    Te /= Tnorm;
    Ti /= Tnorm;
    SAVE_ONCE2(Te, Ti);

    /////////////////////////////////////////////////////////////////////
    // Collisional coefficients
    
    OPTION(opt, resistivity, 0.0);
    OPTION(opt, viscosity, 0.0);     // [m^2/s] . Normalised later
    OPTION(opt, viscosity_par, 0.0); // [m^2/s] . Normalised later
    
    BoutReal Coulomb = 6.6 - 0.5*log(Nnorm/1e20) + 1.5*log(Te*Tnorm); // Coulomb logarithm
    
    output.write("\tCoulomb logarithm = %e\n", Coulomb);
    // Electron collision time [seconds]
    BoutReal tau_e = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Te*Tnorm, -3./2));
    // ion-ion collision time [seconds]
    BoutReal tau_i = sqrt(AA) / (4.80e-8 * (Nnorm / 1e6) * Coulomb * pow(Te*Tnorm, -3./2));
    
    output.write("\tCollision times [sec]: tau_e = %e, tau_i = %e\n", tau_e, tau_i);
    
    // Parallel conductivity
    BoutReal sigma_par = 1.96*Nnorm*SQ(SI::qe)*tau_e/SI::Me;

    output.write("\tBraginskii resistivity: %e [Ohm m]\n", 1./sigma_par);
    output.write("\tUsing resistivity: %e\n", resistivity);
    
    // Perpendicular viscosity
    BoutReal eta_perp = 0.3 * Nnorm*SI::qe*Te*Tnorm/ ( SQ(Omega_ci)*tau_i );
    
    // Perpendicular gyro-viscosity
    BoutReal eta_gyro = 0.5*Nnorm*SI::qe*Te/Omega_ci;
    
    output.write("\tViscosity: %e [Pa s], Gyro-viscosity: %e [Pa s]\n", eta_perp, eta_gyro);
    output.write("\tKinematic viscosity: %e [m^2/s], Gyro-viscosity: %e [m^2/s]\n", eta_perp/(mi*Nnorm), eta_gyro/(mi*Nnorm));
    output.write("\tUsing viscosity: %e\n", viscosity);

    BoutReal eta_par = tau_i * SQ(Cs)*(Ti/Tnorm);
    output.write("\tParallel kinematic viscosity: %e [m^2/s]\n", eta_par);
    output.write("\tUsing parallel viscosity: %e\n", viscosity_par);
    
    // Normalise
    viscosity     *= Bnorm / Tnorm;
    viscosity_par *= Bnorm / Tnorm;
    resistivity *= SI::qe*Nnorm / Bnorm;
    
    output.write("\tNormalised viscosity = %e, parallel viscosity = %e, resistivity = %e\n", viscosity, viscosity_par, resistivity);

    OPTION(opt, vacuum_density, 2e-2);
    OPTION(opt, vacuum_trans, 5e-4);
    OPTION(opt, vacuum_mult, 1e6);
    OPTION(opt, vacuum_damp, 1.0);
    
    /////////////////////////////////////////////////////////////////////
    // Coordinates and mesh
    
    Coordinates *coord = mesh->coordinates();
    
    mesh->get(Bxyz, "B");
    SAVE_ONCE(Bxyz);

    mesh->communicate(Bxyz); // To get yup/ydown fields
    
    logB = log(Bxyz);
    SAVE_ONCE(logB);
    
    // Metric factor appearing in front of the bracket operator
    // For Clebsch coordinates this is 1
    bracket_factor = sqrt(coord->g_yy) / (coord->J * Bxyz);
    mesh->communicate(bracket_factor);
    // Normalise by leaving metric as-is (SI units), but dividing dx,dy and dz by rho_s
    coord->dx /= rho_s;
    coord->dy /= rho_s;
    coord->dz /= rho_s;
    
    R = sqrt(coord->g_yy) / rho_s;

    inv_dy = 1. / (sqrt(coord->g_yy) * coord->dy);
    
    /////////////////////////////////////////////////////////////////////
    // Equilibrium profiles
    
    // Read parallel current density
    mesh->get(J0, "Jpar");
    J0 /= SI::qe*Nnorm*Cs; // Normalise
    
    // Read equilibrium pressure
    Field3D P0;
    mesh->get(P0, "pressure"); // In Pascals
    P0 /= SI::qe*Tnorm*Nnorm;
    
    P0 = 0.0;
    
    // Calculate density from the pressure, assuming isothermal
    N0 = P0 / (Te + Ti);
    
    SAVE_ONCE2(J0, N0);
    
    /////////////////////////////////////////////////////////////////////
    //
    // Solve for electromagnetic potential, vorticity, density and parallel momentum
    
    OPTION(opt, electromagnetic, true);
    OPTION(opt, FiniteElMass, false);

    OPTION(opt, Diffusion, false);
    OPTION(opt, curvilinear, false);

    if(curvilinear){
      mesh->periodicZ=true;
    }
      
    
    omega.mergeYupYdown();
    n.mergeYupYdown();
    nvi.mergeYupYdown();
    Ajpar.mergeYupYdown();
    
    SOLVE_FOR3(omega, n, nvi);

    if (electromagnetic || FiniteElMass) {
      SOLVE_FOR(Ajpar);
    }
    
    /////////////////////////////////////////////////////////////////////
    // Laplacian inversions
    
    // Create Laplacian inversion objects for potentials
    phiSolver = LaplaceXZ::create(mesh, Options::getRoot()->getSection("phisolver"));
    phi = psi = Apar = Jpar = 0.0; // Initial value
    phiSolver->setCoefs(Field3D(1.0),Field3D(0.0));

    if (electromagnetic && FiniteElMass) {
      // Need to solve a Helmholtz problem to get Apar
      psiSolver = LaplaceXZ::create(mesh, Options::getRoot()->getSection("psisolver"));
      Field3D A = 1.0;
      Field3D D = -1.0/(beta * mi_me);
      psiSolver->setCoefs(A,D);
  
      SAVE_REPEAT(Apar);
    }

    // // Laplacian inversions
    
    // // Create Laplacian inversion objects for potentials
    // phiSolver = Laplacian::create(Options::getRoot()->getSection("phisolver"));
    // phi = psi = Apar = Jpar = 0.0; // Initial value

    // if (electromagnetic && FiniteElMass) {
    //   // Need to solve a Helmholtz problem to get Apar
    //   psiSolver = Laplacian::create(Options::getRoot()->getSection("psisolver"));
    //   psiSolver->setCoefA(1.0); // density n
    //   psiSolver->setCoefD(-1.0/(beta * mi_me));

    //   SAVE_REPEAT(Apar);
    // }
    
    // Additional outputs
    SAVE_REPEAT3(phi, Jpar, psi);
    
    SAVE_REPEAT2(ddt(nvi), ddt(n));
    
    n.splitYupYdown();
    SAVE_REPEAT2(n.yup(), n.ydown());
    
    SAVE_REPEAT(vac_mask);

    SAVE_REPEAT3(a,b,c);
    SAVE_REPEAT2(d,f);
    return 0;
  }
  
  int rhs(BoutReal t) {
    printf("TIME = %e\r", t);

    Coordinates *coord = mesh->coordinates();

    ddt(omega) = 0.0;
    mesh->communicate(Apar, omega, n, nvi);
    
    // Apply Dirichlet boundary conditions in z
    for(int i=0;i<mesh->LocalNx;i++) {
      for(int j=0;j<mesh->LocalNy;j++) {
        Apar(i,j,0) = -Apar(i,j,1);
        Apar(i,j,mesh->LocalNz-1) = -Apar(i,j,mesh->LocalNz-2);
        
        omega(i,j,0) = -omega(i,j,1);
        omega(i,j,mesh->LocalNz-1) = -omega(i,j,mesh->LocalNz-2);
      }
    }
    
    ////////////////////////////////////
    // Get phi from vorticity
    phi = phiSolver->solve(omega, 0.0);
    mesh->communicate(phi);
    
    phi -= Ti * n; // Ion diamagnetic drift
    
    Field3D ntot = floor(n + N0, 1e-8);    // Total density
    Field3D jtot = Jpar + J0; // Total parallel current
    Field3D phitot = phi - Ti * N0; // Total electric potential, assuming omega0 = 0
    
    // Parallel flow
    //Field3D vi = nvi / ntot;
    Field3D vi = nvi;
    
    Field3D ptot = (Te + Ti) * ntot; // Total pressure

    // 1 in the vacuum, 0 where there is plasma
    vac_mask = (1.0 - tanh( (ntot - vacuum_density) / vacuum_trans )) / 2.0;
    
    ////////////////////////////////////
    // Ohm's law

    if (electromagnetic) {
      if (!FiniteElMass) {
        // Electromagnetic + zero electron mass
        // Get J from Apar
        Apar = Ajpar;
	const Field3D &beta_inv = -(1./beta);
        Jpar = coord->Div_Perp_Lap_FV(beta_inv, Apar, false); //-(1./beta) * Delp2(Apar);
      } else {
        // Electromagnetic + finite electron mass
        Apar = psiSolver->solve(ntot*Ajpar + nvi/mi_me, 0.0);
        Jpar = ntot*mi_me*(Ajpar - Apar) + nvi;
      }
    } else {
      // Electrostatic
      if (FiniteElMass) {
        // Ajpar = -me_mi v_||e
        Jpar = ntot * Ajpar*mi_me + nvi;
      } else {
        // Electrostatic + zero electron mass
        Jpar = Grad_parP(Te*ntot - phitot) / (resistivity*(1.0 + vac_mask*vacuum_mult));
      }
    }
    mesh->communicate(Jpar);
    Jpar.applyBoundary("neumann");
    for(int i=0;i<mesh->LocalNx;i++) {
      for(int j=0;j<mesh->LocalNy;j++) {
        Jpar(i,j,0) = Jpar(i,j,1);
        Jpar(i,j,mesh->LocalNz-1) = Jpar(i,j,mesh->LocalNz-2);
      }
    }
    
    // Poloidal flux
    psi = Apar * R;
    
    ////////////////////////////////////
    
    // Vorticity
    {
      TRACE("ddt(omega)");
      ddt(omega) = 
        poisson(omega, phitot)   // ExB advection
        - Div_parP2(omega, vi)   // Parallel advection
        + Div_parP(jtot)          // Parallel current
        + curvature(ptot)         // Diamagnetic current (ballooning drive)
        + (1.0 + vac_mask*vacuum_mult)*coord->Div_Perp_Lap_FV(viscosity,omega, false)  //viscosity*(1.0 + vac_mask*vacuum_mult)*Delp2(omega)  // Viscosity
        - vacuum_damp*vac_mask*omega // Damping in vacuum region
        ;
      
      ddt(omega) *= 1.0 - vac_mask;

      c = - curvature(ptot); // Debugging variable

    }
    
    // Perturbed vector potential
    if (electromagnetic || FiniteElMass) {
      TRACE("ddt(Ajpar)");
      
      ddt(Ajpar) = 
        Grad_parP(Te*ntot - phitot)
        - resistivity*(1.0 + vac_mask*vacuum_mult) * Jpar // Resistivity
        ;
    }
    
    // Perturbed density
    {
      TRACE("ddt(n)");
      
      ddt(n) = 
        poisson(ntot, phitot) // ExB advection
        - Div_parP2(ntot, vi) // Parallel advection
        + Div_parP(jtot)
        + curvature(ntot * Te)  // Electron curvature drift
        + coord->Div_Perp_Lap_FV(viscosity,n, false)//viscosity*Delp2(n)
        ;
      f =   curvature(ntot * Te);

    }
    
    // Parallel momentum
    {
      TRACE("ddt(nvi)");
      
      ddt(nvi) =
        poisson(nvi, phitot) // ExB advection
        - Div_parP2(nvi, vi) // Parallel advection
        - curvature(nvi * Ti)  // Ion curvature drift
        - Grad_parP(ptot)    // pressure gradient
        + coord->Div_Perp_Lap_FV(viscosity,nvi, false) //viscosity*Delp2(nvi) // Perpendicular viscosity
        + viscosity_par * Diffusion_parP(vi) // Parallel viscosity
        ;
    }

    // Parallel Diffusion (for blob initial conditions)
    if (Diffusion){
      ddt(n) =
	1e4*Diffusion_parP(n);
      ddt(nvi) = 0.0;
      ddt(Ajpar) = 0.0;
      ddt(omega) = 0.0;
    }
    return 0;
  }

  /// Magnetic curvature term C(f)
  Field3D curvature(const Field3D &f) {
    return 2 * bracket(logB, f, BRACKET_ARAKAWA) * bracket_factor;
  }

  /// Poisson bracket
  /// Note that the Bxyz factor is here because the coordinates are not flux coordinates
  Field3D poisson(const Field3D f, const Field3D &g) {
    return bracket(f, g, BRACKET_ARAKAWA) * bracket_factor;
  }

  // Parallel gradient
  Field3D Grad_parP(const Field3D &f) {
    TRACE("Grad_parP");
    Field3D yup, ydown;
    
    if (!f.hasYupYdown()) {
      // Communicate to get yup and ydown fields
      Field3D fcom = f;
      mesh->communicate(fcom);
      fcom.applyParallelBoundary("parallel_neumann");
      
      yup = fcom.yup();
      ydown = fcom.ydown();
    } else {
      yup = f.yup();
      ydown = f.ydown();
    }

    Field3D result;
    result.allocate();
    
    for(auto &i : result.region(RGN_NOBNDRY)) {
      result[i] = (yup[i.yp()] - ydown[i.ym()]) * 0.5 * inv_dy[i];
    }
    
    return result;
  }

  // Parallel diffusion
  Field3D Diffusion_parP(const Field3D &f) {
    TRACE("Diffusion_parP");
    Field3D yup, ydown;
    
    if (!f.hasYupYdown()) {
      // Communicate to get yup and ydown fields
      Field3D fcom = f;
      mesh->communicate(fcom);
      fcom.applyParallelBoundary("parallel_neumann");
      
      yup = fcom.yup();
      ydown = fcom.ydown();
    } else {
      yup = f.yup();
      ydown = f.ydown();
    }
    
    Field3D Bup = Bxyz.yup();
    Field3D Bdown = Bxyz.ydown();
    
    Field3D result;
    result.allocate();
    
    for(auto &i : result.region(RGN_NOBNDRY)) {
      auto yp = i.yp();
      auto ym = i.ym();

      // Magnetic field half-way between cells
      // This is because the parallel gradient doesn't include
      // gradients of B
      BoutReal Bp = 0.5*(Bup[yp] + Bxyz[i]);
      BoutReal Bm = 0.5*(Bup[ym] + Bxyz[i]);
      
      result[i] = ( (yup[yp]/Bp) - f[i]*(1./Bp + 1./Bm) + (ydown[ym]/Bm)) * SQ(inv_dy[i]) * Bxyz[i];
    }
    
    return result;
  }
  
  // Parallel divergence
  Field3D Div_parP(const Field3D &f) {
    TRACE("Div_parP");

    Field3D yup, ydown;
    
    if (!f.hasYupYdown()) {
      // Communicate to get yup and ydown fields
      Field3D fcom = f;
      mesh->communicate(fcom);
      fcom.applyParallelBoundary("parallel_neumann");
      
      yup = fcom.yup();
      ydown = fcom.ydown();
    } else {
      yup = f.yup();
      ydown = f.ydown();
    }
    
    Field3D result;
    result.allocate();
    Field3D Bup = Bxyz.yup();
    Field3D Bdown = Bxyz.ydown();
    
    for(auto &i : result.region(RGN_NOBNDRY)) {
      auto yp = i.yp();
      auto ym = i.ym();
      result[i] = ( (yup[yp] / Bup[yp]) - (ydown[ym] / Bdown[ym]) )  * 0.5 * inv_dy[i] * Bxyz[i];
    }

    return result;
  }
  
  // Parallel divergence of a product Div_par(fv)
  // Written in skew-symmetric form to suppress zigzag modes
  Field3D Div_parP2(const Field3D f, const Field3D v) {
    return 0.5*(Div_parP(f*v) + v*Grad_parP(f) + f*Div_parP(v));
  }
  
private:
  // Evolving variables
  Field3D Ajpar, omega, n, nvi;  // Electromagnetic potential, vorticity, density perturbation, parallel momentum

  Field3D Apar;
  Field3D Jpar;   // Perturbed parallel current density
  Field3D phi;    // Electrostatic potential
  Field3D psi;    // Poloidal flux

  Field3D Bxyz; // 3D Magnetic field [T]
  Field3D R;    // Major radius  [Normalised]
  Field3D J0;   // Equilibrium parallel current [Normalised]
  Field3D N0;   // Equilibrium density [Normalised]
  
  Field3D bracket_factor; // sqrt(g_yy) / (JB) factor appearing in ExB advection
  Field3D logB; // Log(Bxyz) used in curvature

  Field3D dRdx, dRdz;
  
  LaplaceXZ *phiSolver; // Solver for potential phi from vorticity
  LaplaceXZ *psiSolver; // Solver for potential Apar(psi) from Ajpar
  
  BoutReal resistivity; // [Normalised]
  BoutReal viscosity;   // Perpendicular viscosity [Normalised]
  BoutReal viscosity_par; // Parallel viscosity [Normalised]
  
  BoutReal Nnorm; // Normalisation density [m^-3]
  BoutReal Tnorm; // Normalisation temperature [eV]
  BoutReal Bnorm; // Normalisation magnetic field [T]
  BoutReal AA; // Atomic mass number
  
  BoutReal Omega_ci; // Normalisation frequency [s^-1]
  BoutReal Cs, rho_s; // Normalisation sound speed [m/s], length [m]
  BoutReal beta;     // Normalisation beta
  BoutReal mi_me;    // Ion mass / electron mass
  
  BoutReal Te, Ti;  // Electron and ion temperatures [Normalised]
  
  Field3D inv_dy; // 1 / (sqrt(g_22) * dy) for parallel gradients

  Field3D vac_mask; // 1 in vacuum, 0 in plasma
  BoutReal vacuum_density, vacuum_trans; // Determines the transition from 0 to 1
  BoutReal vacuum_mult; // Multiply dissipation in vacuum
  BoutReal vacuum_damp; // Damping in vacuum region

  bool electromagnetic; // Include electromagnetic effects?
  bool FiniteElMass;    // Finite electron mass?

  bool Diffusion;   //Simulate only parallel diffusion?
  bool curvilinear;   //Switch for curvilinear poloidal grids

  Field3D a,b,c,d,e,f;  //Debugging Variables
};

BOUTMAIN(MergingFlux);
