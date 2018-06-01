
#include "bout/physicsmodel.hxx"
#include "bout/constants.hxx"
#include "invert_laplace.hxx"

#include "parallel_boundary_op.hxx"
#include "boundary_factory.hxx"
#include "field_factory.hxx"

class ParDirichletMidpoint : public BoundaryOpPar {
public:
  ParDirichletMidpoint() :
    BoundaryOpPar(NULL, 0.) {}
  ParDirichletMidpoint(BoundaryRegionPar *region) :
    BoundaryOpPar(region, 0.) {}
  ParDirichletMidpoint(BoundaryRegionPar *region, std::shared_ptr<FieldGenerator>  value) :
    BoundaryOpPar(region, value) {}
  ParDirichletMidpoint(BoundaryRegionPar *region, Field3D* value) :
    BoundaryOpPar(region, value) {}
  ParDirichletMidpoint(BoundaryRegionPar *region, BoutReal value) :
    BoundaryOpPar(region, value) {}
  
  BoundaryOpPar* clone(BoundaryRegionPar *region, const list<string> &args);
  BoundaryOpPar* clone(BoundaryRegionPar *region, Field3D *f);

  using BoundaryOpPar::apply;
  void apply(Field3D &f) override {return apply(f, 0);}
  void apply(Field3D &f, BoutReal t) override;

};

const Field3D Div_a_Laplace_xz(const Field3D &a, const Field3D &f);
const Field3D Div_a_Laplace_xz(const Field3D &f);

//Field3D a,b,c,d;

class MergingFlux : public PhysicsModel {
protected:
  int init(bool restarting) {

    // SAVE_REPEAT4(a,b,c,d);
    
    // Add custom boundary condition
    BoundaryFactory::getInstance()->add(new ParDirichletMidpoint(), "parallel_dirichlet_midpoint");
    
    // Read options
    Options *opt = Options::getRoot()->getSection("model");
    
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
    OPTION(opt, diffusion, viscosity); // [m^2/s] . Normalised later
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
    diffusion *= Bnorm / Tnorm;
    
    output.write("\tNormalised viscosity = %e, parallel viscosity = %e, resistivity = %e\n", viscosity, viscosity_par, resistivity);

    OPTION(opt, hyperresist, 1.0);
    
    /////////////////////////////////////////////////////////////////////
    // Coordinates and mesh
    
    Coordinates *coord = mesh->coordinates();
    
    mesh->get(Bxyz, "B");
    SAVE_ONCE(Bxyz);

    Bxyz.applyBoundary("neumann");
    ASSERT1( min(Bxyz) > 0.0 );
    
    mesh->communicate(Bxyz); // To get yup/ydown fields

    // Note: A Neumann condition simplifies boundary conditions on fluxes
    // where the condition e.g. on J should be on flux (J/B)
    Bxyz.applyParallelBoundary("parallel_neumann");
    
    logB = log(Bxyz);
    
    // Metric factor appearing in front of the bracket operator
    // For Clebsch coordinates this is 1
    bracket_factor = sqrt(coord->g_22) / (coord->J * Bxyz);
    
    // Normalise by leaving metric as-is (SI units), but dividing dx,dy and dz by rho_s
    coord->dx /= rho_s;
    coord->dy /= rho_s;
    coord->dz /= rho_s;
    
    R = sqrt(coord->g_22) / rho_s;
    
    inv_dy = 1. / (sqrt(coord->g_22) * coord->dy);
    
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

    OPTION(opt, drifts, true); // Include drifts, electric fields?
    OPTION(opt, electromagnetic, true);
    OPTION(opt, FiniteElMass, false);

    OPTION(opt, n_div_integrate, false);
    OPTION(opt, nvi_div_integrate, false);
    
    SOLVE_FOR2(n, nvi);
    
    if (drifts) {
      SOLVE_FOR(omega);
      if (electromagnetic || FiniteElMass) {
        SOLVE_FOR(Ajpar);
      }
      SAVE_REPEAT3(phi, Jpar, psi);
    } else {
      omega = 0.0;
    }
    
    /////////////////////////////////////////////////////////////////////
    // Laplacian inversions
    
    // Create Laplacian inversion objects for potentials
    phiSolver = Laplacian::create(Options::getRoot()->getSection("phisolver"));
    phi = psi = Apar = Jpar = 0.0; // Initial value

    if (electromagnetic && FiniteElMass) {
      // Need to solve a Helmholtz problem to get Apar
      psiSolver = Laplacian::create(Options::getRoot()->getSection("psisolver"));
      psiSolver->setCoefA(1.0); // density n
      psiSolver->setCoefD(-1.0/(beta * mi_me));

      SAVE_REPEAT(Apar);
    }

    ntot.setBoundary("ntot");
    jtot.setBoundary("jtot");
    vi.setBoundary("vi");
    
    return 0;
  }
  
  int rhs(BoutReal t) {

    if (drifts) {
      mesh->communicate(Apar, omega, n, nvi);
      Apar.applyParallelBoundary();
      omega.applyParallelBoundary();
      
    } else {
      mesh->communicate(n, nvi);
    }
    
    n.applyParallelBoundary();
    // Note: Boundary applied to vi rather than nvi
    
    // Apply Dirichlet boundary conditions in z
    for(int i=0;i<mesh->LocalNx;i++) {
      for(int j=0;j<mesh->LocalNy;j++) {
        Apar(i,j,0) = -Apar(i,j,1);
        Apar(i,j,mesh->LocalNz-1) = -Apar(i,j,mesh->LocalNz-2);
        
        omega(i,j,0) = -omega(i,j,1);
        omega(i,j,mesh->LocalNz-1) = -omega(i,j,mesh->LocalNz-2);
      }
    }

    if (drifts) {
      ////////////////////////////////////
      // Get phi from vorticity
      phi = phiSolver->solve(omega, 0.0);
      mesh->communicate(phi);
      
      phi -= Ti * n; // Ion diamagnetic drift
    } else {
      phi = 0.0;
    }
    
    ntot = floor(n + N0, 0.0);    // Total density
      
    Field3D ntot_lim = floor(ntot, 1e-4);
    
    Field3D phitot = phi - Ti * N0; // Total electric potential, assuming omega0 = 0
    
    // Parallel flow
    vi = nvi / ntot_lim;
    Field3D momflux = nvi * vi; // Momentum flux
    
    // Apply boundary conditions to v
    vi.splitYupYdown();
    vi.yup() = nvi.yup() / ntot_lim;
    vi.ydown() = nvi.ydown() / ntot_lim;
    vi.applyParallelBoundary();

    // Ensure that boundary conditions are consistent
    // between v, nv and momentum flux
    
    momflux.splitYupYdown();
    for (const auto &reg : mesh->getBoundariesPar()) {
      // Using the values of density and velocity on the boundary
      const Field3D &n_next = n.ynext(reg->dir);
      const Field3D &vi_next = vi.ynext(reg->dir);

      // Set the momentum and momentum flux
      Field3D &nvi_next = nvi.ynext(reg->dir);
      Field3D &momflux_next = momflux.ynext(reg->dir);
      momflux_next.allocate();
        
      for (reg->first(); !reg->isDone(); reg->next()) {
        // Density at the boundary
        BoutReal n_b = 0.5*(n_next(reg->x, reg->y+reg->dir, reg->z) +
                            n(reg->x, reg->y, reg->z));
        if (n_b < 0.0) {
          n_b = 0.0;
        }
        // Velocity at the boundary
        BoutReal vi_b = 0.5*(vi_next(reg->x, reg->y+reg->dir, reg->z) +
                             vi(reg->x, reg->y, reg->z));
        
        nvi_next(reg->x, reg->y+reg->dir, reg->z) = 
          2.*n_b*vi_b - nvi(reg->x, reg->y, reg->z);

        momflux_next(reg->x, reg->y+reg->dir, reg->z) = 
          2.*n_b*vi_b*vi_b - momflux(reg->x, reg->y, reg->z);
      }
    }
    
    Field3D ptot = (Te + Ti) * ntot; // Total pressure
    
    ////////////////////////////////////
    // Ohm's law

    if (drifts) {
      if (electromagnetic) {
        if (!FiniteElMass) {
          // Electromagnetic + zero electron mass
          // Get J from Apar
          Apar = Ajpar;
          Jpar = -(1./beta) * Div_a_Laplace_xz(Apar);
        } else {
          // Electromagnetic + finite electron mass
          Apar = psiSolver->solve(Ajpar + nvi/mi_me, 0.0);
          Jpar = mi_me*ntot*(Ajpar - Apar) + nvi;
        }
      } else {
        // Electrostatic
        if (FiniteElMass) {
          // Ajpar = -me_mi v_||e
          Jpar = ntot * Ajpar*mi_me + nvi;
        } else {
          // Electrostatic + zero electron mass
          Jpar = Grad_parP(Te*ntot - phitot) / resistivity;
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
    } else {
      Jpar = Apar = Ajpar = 0.0;

      phitot = 0.0;
    }
    
    // Poloidal flux
    psi = Apar * R;

    jtot = Jpar + J0; // Total parallel current
    
    ntot.applyBoundary();
    jtot.applyBoundary();

    //jtot *= Bxyz; // Note:Applying boundary on flux (j/B)
    mesh->communicate(ntot, jtot);
    ntot.applyParallelBoundary();
    
    jtot.applyParallelBoundary(0.0);
    // jtot /= Bxyz;
    // jtot.yup() /= Bxyz.yup();
    // jtot.ydown() /= Bxyz.ydown();

    if (!drifts) {
      jtot = 0.0;
    }
    
    ////////////////////////////////////

    if (drifts) {
      // Vorticity
      {
        TRACE("ddt(omega)");
        ddt(omega) = 
          poisson(omega, phitot)   // ExB advection
          - Div_parP2(omega, vi)   // Parallel advection
          + Div_parP(jtot)          // Parallel current
          + curvature(ptot)         // Diamagnetic current (ballooning drive)
          + viscosity*Div_a_Laplace_xz(omega)  // Viscosity
          ;
      }
    
      // Perturbed vector potential
      if (electromagnetic || FiniteElMass) {
        TRACE("ddt(Ajpar)");
        
        Field3D Jpar2 = Div_a_Laplace_xz(Jpar);
        
        ddt(Ajpar) = 
          Grad_parP(Te*ntot - phitot)
          - resistivity * Jpar // Resistivity
          + hyperresist * Jpar2  // Hyper-resistivity
          ;
      }
    }
    
    // Perturbed density
    {
      TRACE("ddt(n)");
      
      ddt(n) = diffusion*Div_a_Laplace_xz(n)
        ;

      if (n_div_integrate) {
        ddt(n) -= Div_par_integrate(nvi);
        
        if (drifts) {
          /*
          Field3D jflux = jtot/ntot_lim;

          jflux.splitYupYdown();
          for (const auto &reg : mesh->getBoundariesPar()) {
            // Using the values of density and velocity on the boundary
            const Field3D &ntot_next = ntot_lim.ynext(reg->dir);
            const Field3D &jtot_next = jtot.ynext(reg->dir);
            
            // Set the momentum and momentum flux
            Field3D &nvi_next = nvi.ynext(reg->dir);
            Field3D &momflux_next = momflux.ynext(reg->dir);
            momflux_next.allocate();
            
            for (reg->first(); !reg->isDone(); reg->next()) {
              // Density at the boundary
              BoutReal n_b = 0.5*(n_next(reg->x, reg->y+reg->dir, reg->z) +
                                  n(reg->x, reg->y, reg->z));
              // Velocity at the boundary
              BoutReal vi_b = 0.5*(vi_next(reg->x, reg->y+reg->dir, reg->z) +
                                   vi(reg->x, reg->y, reg->z));
              
              nvi_next(reg->x, reg->y+reg->dir, reg->z) = 
                2.*n_b*vi_b - nvi(reg->x, reg->y, reg->z);
              
              momflux_next(reg->x, reg->y+reg->dir, reg->z) = 
                2.*n_b*vi_b*vi_b - momflux(reg->x, reg->y, reg->z);
            }
          }
          */
          ddt(n) -= Div_par_U1(ntot, -jtot/ntot_lim);
        }
        
      } else {
        ddt(n) -= Div_par_U1(ntot, vi - jtot/ntot_lim);
      }
      
      if (drifts) {
        ddt(n) +=
          poisson(ntot, phitot) // ExB advection
          + curvature(ntot * Te)  // Electron curvature drift
          ;
      }
    }
    
    // Parallel momentum
    {
      TRACE("ddt(nvi)");
      
      ddt(nvi) =
        - (Te + Ti) * Grad_par(n)
        //- Grad_parP((Te + Ti) * (n + N0)) // Pressure gradient, no flooring at 0
        
        + diffusion*Div_a_Laplace_xz(vi, n)    // Perpendicular diffusion
        + viscosity*Div_a_Laplace_xz(ntot,vi)     // Perpendicular viscosity
        + viscosity_par * Diffusion_parP(ntot, vi) // Parallel viscosity
        //+ viscosity_par * Grad2_par2(nvi)
        ;
      
      // Parallel advection
      if (nvi_div_integrate) {
        ddt(nvi) -=  Div_par_integrate(momflux);
      } else {
        ddt(nvi) -= Div_par_U1(nvi, vi);
      }
      
      if (drifts) {
        ddt(nvi) +=
          poisson(nvi, phitot) // ExB advection
          - curvature(nvi * Ti)  // Ion curvature drift
          ;
      }
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
      //fcom.applyParallelBoundary("parallel_neumann");
      fcom.applyParallelBoundary("parallel_dirichlet_midpoint");
      
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
  Field3D Diffusion_parP(const Field3D &D, const Field3D &f) {
    TRACE("Diffusion_parP");
    Field3D yup, ydown;
    Field3D Dup, Ddown;
    
    if (!f.hasYupYdown()) {
      // Communicate to get yup and ydown fields
      Field3D fcom = f;
      mesh->communicate(fcom);
      fcom.applyParallelBoundary("parallel_neumann");
      //fcom.applyParallelBoundary("parallel_dirichlet_midpoint");
      
      yup = fcom.yup();
      ydown = fcom.ydown();
    } else {
      yup = f.yup();
      ydown = f.ydown();
    }
    
    if (!D.hasYupYdown()) {
      // Communicate to get yup and ydown fields
      Field3D Dcom = D;
      mesh->communicate(Dcom);
      Dcom.applyParallelBoundary("parallel_neumann");
      
      Dup = Dcom.yup();
      Ddown = Dcom.ydown();
    } else {
      Dup = D.yup();
      Ddown = D.ydown();
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
      BoutReal Bm = 0.5*(Bdown[ym] + Bxyz[i]);

      BoutReal Dp = 0.5*(Dup[yp] + D[i]);
      BoutReal Dm = 0.5*(Ddown[ym] + D[i]);
      
      result[i] = ( (yup[yp]*Dp/Bp) - f[i]*(Dp/Bp + Dm/Bm) + (ydown[ym]*Dm/Bm)) * SQ(inv_dy[i]) * Bxyz[i];

      ASSERT3(finite(result[i]));
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
      //fcom.applyParallelBoundary("parallel_neumann");
      fcom.applyParallelBoundary("parallel_dirichlet_midpoint");
      
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
  
  /// Parallel divergence of a product Div_par(fv)
  /// Written in skew-symmetric form to suppress zigzag modes
  Field3D Div_parP2(const Field3D f, const Field3D v) {
    return 0.5*(Div_parP(f*v) + v*Grad_parP(f) + f*Div_parP(v));
  }

  /// First order upwinding method
  Field3D Div_par_U1(const Field3D f, const Field3D v) {
    Field3D result;
    result.allocate();

    Field3D fyup, fydown;
    
    if (!f.hasYupYdown()) {
      // Communicate to get yup and ydown fields
      Field3D fcom = f;
      mesh->communicate(fcom);
      //fcom.applyParallelBoundary("parallel_neumann");
      fcom.applyParallelBoundary("parallel_dirichlet_midpoint");
      
      fyup = fcom.yup();
      fydown = fcom.ydown();
    } else {
      fyup = f.yup();
      fydown = f.ydown();
    }

    Field3D vyup, vydown;
    
    if (!v.hasYupYdown()) {
      // Communicate to get yup and ydown fields
      Field3D vcom = v;
      mesh->communicate(vcom);
      vcom.applyParallelBoundary("parallel_dirichlet_midpoint");
      
      vyup = vcom.yup();
      vydown = vcom.ydown();
    } else {
      vyup = v.yup();
      vydown = v.ydown();
    }
    
    Field3D Bup = Bxyz.yup();
    Field3D Bdown = Bxyz.ydown();
    
    for(auto &i : result.region(RGN_NOBNDRY)) {
      auto yp = i.yp();
      auto ym = i.ym();
      
      // Lower boundary
      BoutReal vm = 0.5*(vydown[ym] + v[i]); // Velocity at lower boundary
      BoutReal fm = (vm > 0.0) ? fydown[ym] : f[i]; // Field at lower boundary
      BoutReal Bm = 0.5*(Bdown[ym] + Bxyz[i]);  // Magnetic field
      
      // Upper boundary
      BoutReal vp = 0.5*(vyup[yp] + v[i]);
      BoutReal fp = (vp > 0.0) ? f[i] : fyup[yp];
      BoutReal Bp = 0.5*(Bup[yp] + Bxyz[i]);  // Magnetic field
      
      result[i] = (vp * fp / Bp - vm * fm / Bm)*inv_dy[i] * Bxyz[i];
    }
    return result;
  }

  /// Parallel divergence, using integration over projected cells
  Field3D Div_par_integrate(const Field3D &f) {
    Field3D f_B = f / Bxyz;
    
    f_B.splitYupYdown();
    mesh->getParallelTransform().integrateYUpDown(f_B);

    // integrateYUpDown replaces all yup/down points, so the boundary conditions
    // now need to be applied. If Bxyz has neumann parallel boundary conditions
    // then the boundary condition is simpler since f = 0 gives f_B=0 boundary condition.

    /// Loop over the mesh boundary regions
    for (const auto &reg : mesh->getBoundariesPar()) {
      Field3D &f_B_next = f_B.ynext(reg->dir);
      const Field3D &f_next = f.ynext(reg->dir);
      const Field3D &B_next = Bxyz.ynext(reg->dir);
      
      for (reg->first(); !reg->isDone(); reg->next()) {
        f_B_next(reg->x, reg->y+reg->dir, reg->z) =
          f_next(reg->x, reg->y+reg->dir, reg->z) / B_next(reg->x, reg->y+reg->dir, reg->z);
      }
    }
    
    Field3D result;
    result.allocate();
    
    Coordinates *coord = mesh->coordinates();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      result[i] = Bxyz[i] * (f_B.yup()[i.yp()] - f_B.ydown()[i.ym()]) * 0.5 * inv_dy[i];
    }
    
    return result;
  }
  
private:
  // Evolving variables
  Field3D Ajpar, omega, n, nvi;  // Electromagnetic potential, vorticity, density perturbation, parallel momentum

  Field3D Apar;
  Field3D Jpar;   // Perturbed parallel current density
  Field3D phi;    // Electrostatic potential
  Field3D psi;    // Poloidal flux
  Field3D vi;     // Parallel flow velocity

  Field3D jtot; // Total current
  Field3D ntot; // Total density
  
  Field3D Bxyz; // 3D Magnetic field [T]
  Field2D R;    // Major radius  [Normalised]
  Field3D J0;   // Equilibrium parallel current [Normalised]
  Field3D N0;   // Equilibrium density [Normalised]
  
  Field3D bracket_factor; // sqrt(g_yy) / (JB) factor appearing in ExB advection
  Field3D logB; // Log(Bxyz) used in curvature

  Laplacian *phiSolver; // Solver for potential phi from vorticity
  Laplacian *psiSolver; // Solver for potential Apar(psi) from Ajpar
  
  BoutReal resistivity; // [Normalised]
  BoutReal hyperresist; // Hyper-resistivity
  BoutReal diffusion;   // Density diffusion [Normalised]
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
  
  Field2D inv_dy; // 1 / (sqrt(g_22) * dy) for parallel gradients

  bool drifts;          // Include currents and electric fields?
  bool electromagnetic; // Include electromagnetic effects?
  bool FiniteElMass;    // Finite electron mass?

  bool n_div_integrate; // Density parallel flux method
  bool nvi_div_integrate; // Parallel momentum flux method
 
};

///////////////////////////////////////////////////////////
// Parallel boundary condition
// using a midpoint approximation for Dirichlet condition

BoundaryOpPar* ParDirichletMidpoint::clone(BoundaryRegionPar *region, const list<string> &args) {
  if(!args.empty()) {
    try {
      real_value = stringToReal(args.front());
      return new ParDirichletMidpoint(region, real_value);
    } catch (BoutException& e) {
      std::shared_ptr<FieldGenerator>  newgen = 0;
      // First argument should be an expression
      newgen = FieldFactory::get()->parse(args.front());
      return new ParDirichletMidpoint(region, newgen);
    }
  }
  return new ParDirichletMidpoint(region);
}

BoundaryOpPar* ParDirichletMidpoint::clone(BoundaryRegionPar *region, Field3D *f) {
  return new ParDirichletMidpoint(region, f);
}

void ParDirichletMidpoint::apply(Field3D &f, BoutReal t) {

  Field3D& f_next = f.ynext(bndry->dir);

  Coordinates& coord = *(mesh->coordinates());

  // Loop over grid points If point is in boundary, then fill in
  // f_next such that the field would be VALUE on the boundary
  for (bndry->first(); !bndry->isDone(); bndry->next()) {
    // temp variables for convenience
    int x = bndry->x; int y = bndry->y; int z = bndry->z;

    // Generate the boundary value
    BoutReal value = getValue(*bndry, t);
    
    f_next(x, y+bndry->dir, z) = 2.*value - f(x,y,z);
  }
}

BOUTMAIN(MergingFlux);

// Div ( a Laplace_xz(f) )  -- Vorticity
const Field3D Div_a_Laplace_xz(const Field3D &a, const Field3D &f) {
  Mesh *mesh = a.getMesh();

  Field3D result(mesh);
  result = 0.0;

  Coordinates *coord = mesh->coordinates();

  // Flux in x

  int xs = mesh->xstart - 1;
  int xe = mesh->xend;
  
  for (int i = xs; i <= xe; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Calculate flux from i to i+1

        BoutReal fout =
            (coord->J(i, j) * a(i, j, k) * coord->g11(i, j) +
             coord->J(i + 1, j) * a(i + 1, j, k) * coord->g11(i + 1, j)) *
            (f(i + 1, j, k) - f(i, j, k)) /
            (coord->dx(i, j) + coord->dx(i + 1, j));

        result(i, j, k) += fout / (coord->dx(i, j) * coord->J(i, j));
        result(i + 1, j, k) -= fout / (coord->dx(i + 1, j) * coord->J(i + 1, j));
      }
    }
  
  // Z flux
  // Easier since all metrics constant in Z

  for (int i = mesh->xstart; i <= mesh->xend; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      // Coefficient in front of df/dy term
      BoutReal coef =
          coord->g_23(i, j) /
          (coord->dy(i, j + 1) + 2. * coord->dy(i, j) + coord->dy(i, j - 1)) /
          SQ(coord->J(i, j) * coord->Bxy(i, j));

      for (int k = 0; k < mesh->LocalNz; k++) {
        // Calculate flux between k and k+1
        int kp = (k + 1) % (mesh->LocalNz);

        BoutReal fout =
            0.5 * (a(i, j, k) + a(i, j, kp)) * coord->g33(i, j) *
            (
                // df/dz
                (f(i, j, kp) - f(i, j, k)) / coord->dz

                // - g_yz * df/dy / SQ(J*B)
                - coef * (f.yup()(i, j + 1, k) + f.yup()(i, j + 1, kp) -
                          f.ydown()(i, j - 1, k) - f.ydown()(i, j - 1, kp)));

        result(i, j, k) += fout / coord->dz;
        result(i, j, kp) -= fout / coord->dz;
      }
    }
  
  return result;
}

// Div ( Laplace_xz(f) )
const Field3D Div_a_Laplace_xz(const Field3D &f) {
  Mesh *mesh = f.getMesh();

  Field3D result(mesh);
  result = 0.0;

  Coordinates *coord = mesh->coordinates();

  // Flux in x

  int xs = mesh->xstart - 1;
  int xe = mesh->xend;
  
  for (int i = xs; i <= xe; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        // Calculate flux from i to i+1

        BoutReal fout =
            (coord->J(i, j) * coord->g11(i, j) +
             coord->J(i + 1, j) * coord->g11(i + 1, j)) *
            (f(i + 1, j, k) - f(i, j, k)) /
            (coord->dx(i, j) + coord->dx(i + 1, j));

        result(i, j, k) += fout / (coord->dx(i, j) * coord->J(i, j));
        result(i + 1, j, k) -= fout / (coord->dx(i + 1, j) * coord->J(i + 1, j));
      }
    }
  
  // Z flux
  // Easier since all metrics constant in Z

  for (int i = mesh->xstart; i <= mesh->xend; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      // Coefficient in front of df/dy term
      BoutReal coef =
          coord->g_23(i, j) /
          (coord->dy(i, j + 1) + 2. * coord->dy(i, j) + coord->dy(i, j - 1)) /
          SQ(coord->J(i, j) * coord->Bxy(i, j));

      for (int k = 0; k < mesh->LocalNz; k++) {
        // Calculate flux between k and k+1
        int kp = (k + 1) % (mesh->LocalNz);

        BoutReal fout =
            coord->g33(i, j) *
            (
                // df/dz
                (f(i, j, kp) - f(i, j, k)) / coord->dz

                // - g_yz * df/dy / SQ(J*B)
                - coef * (f.yup()(i, j + 1, k) + f.yup()(i, j + 1, kp) -
                          f.ydown()(i, j - 1, k) - f.ydown()(i, j - 1, kp)));

        result(i, j, k) += fout / coord->dz;
        result(i, j, kp) -= fout / coord->dz;
      }
    }
  
  return result;
}
