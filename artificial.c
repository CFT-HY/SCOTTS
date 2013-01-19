/* artificial.c
 *
 * Artificial viscosity.
 */
#include "hydro.h"


void artificial_viscosity(hydro_fields f, hydro_params p) {

  int x, y, z;
  
  double ***gradu = make_field(p);

  double Ux, Uxp, Uy, Uyp, Uz, Uzp;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {


	Ux = 0.25*(f.U[0][x][y][z] + f.U[0][x][y+1][z] + f.U[0][x][y][(z+1)%p.Lz]
		   + f.U[0][x][y+1][(z+1)%p.Lz]);

	Uxp = 0.25*(f.U[0][x+1][y][z] + f.U[0][x+1][y+1][z] + f.U[0][x+1][y][(z+1)%p.Lz]
		    + f.U[0][x+1][y+1][(z+1)%p.Lz]);


	Uy = 0.25*(f.U[1][x][y][z] + f.U[1][x+1][y][z] + f.U[1][x][y][(z+1)%p.Lz]
		   + f.U[1][x+1][y][(z+1)%p.Lz]);

	Uyp = 0.25*(f.U[1][x][y+1][z] + f.U[1][x+1][y+1][z] + f.U[1][x][y+1][(z+1)%p.Lz]
		    + f.U[1][x+1][y+1][(z+1)%p.Lz]);


	Uz = 0.25*(f.U[2][x][y][z] + f.U[2][x+1][y][z] + f.U[2][x][y+1][z] + f.U[2][x+1][y+1][z]);

	Uzp = 0.25*(f.U[2][x][y][(z+1)%p.Lz] f.U[2][x+1][y][(z+1)%p.Lz]
		    + f.U[2][x][y+1][(z+1)%p.Lz] + f.U[2][x+1][y+1][(z+1)%p.Lz]);

	gradu[x][y][z] = (Uxp - Ux + Uyp - Uy + Uzp - Uz)/p.dx;
	
      }
    }
  }

  halo_field(dxphi[0], p);
  halo_field(dxphi[1], p);
  halo_field(dxphi[2], p);

  Vdpot(p, f.T[0][0], phiav[0][0], Vdmid[0][0]);


  halo_field(Vdmid, p);



  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	vdnb = 0.125*(Vdmid[x-1][y][z]
		      + Vdmid[x][y][z]
		      + Vdmid[x][y-1][z]
		      + Vdmid[x-1][y-1][z]
		      + Vdmid[x][y][((z-1+p.Lz)%p.Lz)]
		      + Vdmid[x][y-1][((z+p.Lz-1)%p.Lz)]
		      + Vdmid[x-1][y][((z+p.Lz-1)%p.Lz)]
		      + Vdmid[x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	
	pinb = 0.125*(f.pi[x-1][y][z]
		      + f.pi[x][y][z]
		      + f.pi[x][y-1][z]
		      + f.pi[x-1][y-1][z]
		      + f.pi[x][y][((z-1+p.Lz)%p.Lz)]
		      + f.pi[x][y-1][((z+p.Lz-1)%p.Lz)]
		      + f.pi[x-1][y][((z+p.Lz-1)%p.Lz)]
		      + f.pi[x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	
	wnb = 0.125*(f.W[x-1][y][z]
		     + f.W[x][y][z]
		     + f.W[x][y-1][z]
		     + f.W[x-1][y-1][z]
		     + f.W[x][y][((z-1+p.Lz)%p.Lz)]
		     + f.W[x][y-1][((z+p.Lz-1)%p.Lz)]

		     + f.W[x-1][y][((z+p.Lz-1)%p.Lz)]
		     + f.W[x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	
	dxphinb0 = 0.25*(dxphi[0][x-1][y][z]
			 + dxphi[0][x-1][y-1][z]
			 + dxphi[0][x-1][y][((z+p.Lz-1)%p.Lz)]
			 + dxphi[0][x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	dxphinb1 = 0.25*(dxphi[1][x][y-1][z]
			 + dxphi[1][x-1][y-1][z]
			 + dxphi[1][x][y-1][((z+p.Lz-1)%p.Lz)]
			 + dxphi[1][x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	dxphinb2 = 0.25*(dxphi[2][x][y][((z-1+p.Lz)%p.Lz)]
			 + dxphi[2][x-1][y][((z+p.Lz-1)%p.Lz)]
			 + dxphi[2][x][y-1][((z+p.Lz-1)%p.Lz)]
			 + dxphi[2][x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
       
    f.Z[0][x][y][z] = f.Z[0][x][y][z] - p.dt*p.C*wnb*pinb*dxphinb0;
    f.Z[1][x][y][z] = f.Z[1][x][y][z] - p.dt*p.C*wnb*pinb*dxphinb1;
    f.Z[2][x][y][z] = f.Z[2][x][y][z] - p.dt*p.C*wnb*pinb*dxphinb2;
    
    f.Z[1][x][y][z] = f.Z[1][x][y][z] 
      - p.dt*p.C*wnb*(f.V[1][x][y][z]*dxphinb1
		      + f.V[0][x][y][z]*dxphinb0
		      + f.V[2][x][y][z]*dxphinb2)*dxphinb1;

    f.Z[0][x][y][z] = f.Z[0][x][y][z] 
      - p.dt*p.C*wnb*(f.V[1][x][y][z]*dxphinb1
		      + f.V[0][x][y][z]*dxphinb0
		      + f.V[2][x][y][z]*dxphinb2)*dxphinb0;    

    f.Z[2][x][y][z] = f.Z[2][x][y][z]
      - p.dt*p.C*wnb*(f.V[1][x][y][z]*dxphinb1
		      + f.V[0][x][y][z]*dxphinb0
		      + f.V[2][x][y][z]*dxphinb2)*dxphinb2;

    f.Z[0][x][y][z] = f.Z[0][x][y][z] - p.dt*vdnb*dxphinb0;     
    f.Z[1][x][y][z] = f.Z[1][x][y][z] - p.dt*vdnb*dxphinb1;
    f.Z[2][x][y][z] = f.Z[2][x][y][z] - p.dt*vdnb*dxphinb2;
				       
    gpi = f.W[x][y][z]*(f.pi[x][y][z]
			+ 0.5*(f.V[0][x][y][z]
			       *dxphi[0][x-1][y][z]
			       + f.V[0][x+1][y][z]
			       *dxphi[0][x][y][z])
			+ 0.5*(f.V[1][x][y][z]
                               *dxphi[1][x][y-1][z]
                               + f.V[1][x][y+1][z]
                               *dxphi[1][x][y][z])
			+ 0.5*(f.V[2][x][y][z]
                               *dxphi[2][x][y][((z-1+p.Lz)%p.Lz)]
                               + f.V[2][x][y][((z+1)%p.Lz)]
                               *dxphi[2][x][y][z]));
      
    f.E[x][y][z] = f.E[x][y][z] + p.dt*(p.C*gpi*gpi + gpi*Vdmid[x][y][z]);

    
      }
    }
  }


  // halo_field(f.Z[0], p);
  // halo_field(f.Z[1], p);
  // halo_field(f.Z[2], p);
  //  halo_field(f.E, p);



  // Pressure acceleration
  // W&M sec 2.4.12, 3.5.1, DONE
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// equation (3.68)
	// repeated calculation of these unnecessary, so split it out later
	p_bar_x_plus = (f.p[x][y][z]
			+ f.p[x][y-1][z] 
			+ f.p[x][y][((z-1+p.Lz)%p.Lz)]
			+ f.p[x][y-1][((z+p.Lz-1)%p.Lz)]
			)/4.0;
	
	p_bar_x_minus = (f.p[x-1][y][z] 
			 + f.p[x-1][y-1][z]
			 + f.p[x-1][y][((z+p.Lz-1)%p.Lz)]
			 + f.p[x-1][y-1][((z+p.Lz-1)%p.Lz)]
			 )/4.0;
	
	p_bar_y_plus = (f.p[x][y][z]
			+ f.p[x-1][y][z] 
			+ f.p[x][y][((z-1+p.Lz)%p.Lz)]
			+ f.p[x-1][y][((z+p.Lz-1)%p.Lz)]
			)/4.0;
	
	
	p_bar_y_minus = (f.p[x][y-1][z]
			 + f.p[x-1][y-1][z] 
			 + f.p[x][y-1][((z+p.Lz-1)%p.Lz)]
			 + f.p[x-1][y-1][((z+p.Lz-1)%p.Lz)]
			 )/4.0;
	
	
	p_bar_z_plus = (f.p[x][y][z]
			+ f.p[x-1][y][z] 
			+ f.p[x][y-1][z]
			+ f.p[x-1][y-1][z]
		    )/4.0;
	
	
	p_bar_z_minus = (f.p[x][y][((z-1+p.Lz)%p.Lz)]
			 + f.p[x-1][y][((z+p.Lz-1)%p.Lz)] 
			 + f.p[x][y-1][((z+p.Lz-1)%p.Lz)]
			 + f.p[x-1][y-1][((z+p.Lz-1)%p.Lz)]
			 )/4.0;
	
	
	
	
	f.Z[0][x][y][z] = f.Z[0][x][y][z] - (p_bar_x_plus
					     - p_bar_x_minus)*p.dt/p.dx;
	
	f.Z[1][x][y][z] = f.Z[1][x][y][z] - (p_bar_y_plus
					     - p_bar_y_minus)*p.dt/p.dx;
	
	f.Z[2][x][y][z] = f.Z[2][x][y][z] - (p_bar_z_plus
					     - p_bar_z_minus)*p.dt/p.dx;

      }
    }
  }

  //  halo_field(f.Z[0], p);
  //  halo_field(f.Z[1], p);
  //  halo_field(f.Z[2], p);



  // update velocity v; denominator is W&M eq (2.85)
  // but note grid is Eulerian and D=0
  // then gv is the four-velocity (2.84)
  // (what we call kappa they call sigma, ish?)
  //
  // Section 3.4.5, equations 3.5.7, 3.5.8

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

    /*
    if(x == 0)
      fprintf(stderr,"old Ux = %lf\n", f.U[0][1]);
    */

	sigmabar = (f.kappa[x][y][z]
		    *f.E[x][y][z] 
		    + f.kappa[x-1][y][z]
		    *f.E[x-1][y][z]
		    + f.kappa[x-1][y-1][z]
		    *f.E[x-1][y-1][z]
		    + f.kappa[x-1][y-1][((z+p.Lz-1)%p.Lz)]
		    *f.E[x-1][y-1][((z+p.Lz-1)%p.Lz)]
		    + f.kappa[x][y-1][((z+p.Lz-1)%p.Lz)]
		    *f.E[x][y-1][((z+p.Lz-1)%p.Lz)]
		    + f.kappa[x-1][y][((z+p.Lz-1)%p.Lz)]
		    *f.E[x-1][y][((z+p.Lz-1)%p.Lz)]
		    + f.kappa[x][y-1][z]
		    *f.E[x][y-1][z]
		    + f.kappa[x][y][((z-1+p.Lz)%p.Lz)]
		    *f.E[x][y][((z-1+p.Lz)%p.Lz)]
		)/8.0;

   
    f.U[0][x][y][z] = f.Z[0][x][y][z]/sigmabar;
    f.U[1][x][y][z] = f.Z[1][x][y][z]/sigmabar;
    f.U[2][x][y][z] = f.Z[2][x][y][z]/sigmabar;
    
    //    fprintf(stderr,"%lf %lf\n", f.Z[0][x], sigmabar);
    /*
      if(x == 0)
      fprintf(stderr,"Ux = %lf\n", f.U[0][1]);
    */
      }
    }
  }

  halo_field(f.U[0], p);
  halo_field(f.U[1], p);
  halo_field(f.U[2], p);


  // section 3.4.5 continued
  // face-centred x-cpt of 3-velocity
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	ubarx = (f.U[0][x][y][z]
		 + f.U[0][x][y+1][z] 
		 + f.U[0][x][y][((z+1)%p.Lz)]
		 + f.U[0][x][y+1][((z+1)%p.Lz)]
		 )/4.0;
	
	ubary = (f.U[1][x][y][z]
		 + f.U[1][x][y+1][z]
		 + f.U[1][x][y][((z+1)%p.Lz)]
		 + f.U[1][x][y+1][((z+1)%p.Lz)]
		 )/4.0;
	
	ubarz = (f.U[2][x][y][z]
		 + f.U[2][x][y+1][z]
		 + f.U[2][x][y][((z+1)%p.Lz)]
		 + f.U[2][x][y+1][((z+1)%p.Lz)]
		 )/4.0;
	
	Wfacex[x][y][z] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);
	
	f.V[0][x][y][z] = ubarx/Wfacex[x][y][z];
      }
    }
  }
  



  // y-cpt
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	ubarx = (f.U[0][x][y][z]
		 + f.U[0][x+1][y][z] 
		 + f.U[0][x][y][((z+1)%p.Lz)]
		 + f.U[0][x+1][y][((z+1)%p.Lz)]
		 )/4.0;

	ubary = (f.U[1][x][y][z]
		 + f.U[1][x+1][y][z] 
		 + f.U[1][x][y][((z+1)%p.Lz)]
		 + f.U[1][x+1][y][((z+1)%p.Lz)]
		 )/4.0;

	ubarz = (f.U[2][x][y][z]
		 + f.U[2][x+1][y][z] 
		 + f.U[2][x][y][((z+1)%p.Lz)]
		 + f.U[2][x+1][y][((z+1)%p.Lz)]
		 )/4.0;

	Wfacey[x][y][z] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.V[1][x][y][z] = ubary/Wfacey[x][y][z];
	
      }
    }
  }

  // z-cpt
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	ubarx = (f.U[0][x][y][z]
		 + f.U[0][x][y+1][z] 
		 + f.U[0][x+1][y][z]
		 + f.U[0][x+1][y+1][z]
		 )/4.0;
	
	ubary = (f.U[1][x][y][z]
		 + f.U[1][x][y+1][z] 
		 + f.U[1][x+1][y][z]
		 + f.U[1][x+1][y+1][z]
		 )/4.0;
	
	ubarz = (f.U[2][x][y][z]
		 + f.U[2][x][y+1][z] 
		 + f.U[2][x+1][y][z]
		 + f.U[2][x+1][y+1][z]
		 )/4.0;
	
	Wfacez[x][y][z] = sqrt(1.0 
			       + ubarx*ubarx
			       + ubary*ubary
			       + ubarz*ubarz);
    
	f.V[2][x][y][z] = ubarz/Wfacez[x][y][z];
      }
    }
  }
  

  halo_field(f.V[0], p);
  halo_field(f.V[1], p);
  halo_field(f.V[2], p);



  // End section 3.4.5

  //  fprintf(stderr,"Vx[1] = %lf\n", f.V[0][1]);
  //  getchar();

  // Section 3.4.6
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	
	utildex = (f.U[0][x][y][z] 
		   + f.U[0][x+1][y][z]
		   + f.U[0][x][y+1][z]
		   + f.U[0][x][y][((z+1)%p.Lz)]
		   + f.U[0][x+1][y][((z+1)%p.Lz)]
		   + f.U[0][x][y+1][((z+1)%p.Lz)]
		   + f.U[0][x+1][y+1][z]
		   + f.U[0][x+1][y+1][((z+1)%p.Lz)]
		   )/8.0;
	
	utildey = (f.U[1][x][y][z] 
		   + f.U[1][x+1][y][z]
		   + f.U[1][x][y+1][z]
		   + f.U[1][x][y][((z+1)%p.Lz)]
		   + f.U[1][x+1][y][((z+1)%p.Lz)]
		   + f.U[1][x][y+1][((z+1)%p.Lz)]
		   + f.U[1][x+1][y+1][z]
		   + f.U[1][x+1][y+1][((z+1)%p.Lz)]
		   )/8.0;
    
	utildez = (f.U[2][x][y][z] 
		   + f.U[2][x+1][y][z]
		   + f.U[2][x][y+1][z]
		   + f.U[2][x][y][((z+1)%p.Lz)]
		   + f.U[2][x+1][y][((z+1)%p.Lz)]
		   + f.U[2][x][y+1][((z+1)%p.Lz)]
		   + f.U[2][x+1][y+1][z]
		   + f.U[2][x+1][y+1][((z+1)%p.Lz)]
		   )/8.0;

	
	//    if(x < 3 || x > p.N-3)
	//      fprintf(stderr,"utildex %d = %lf\n", x, utildex);
     
	Wold[x][y][z] = f.W[x][y][z];
	f.W[x][y][z] = sqrt(1.0 
			    + utildex*utildex 
			    + utildey*utildey 
			    + utildez*utildez);
	
	
    //	s = (f.kappa[iix(x,y,z,p)] - 1)*(f.W[iix(x,y,z,p)] - Wold[iix(x,y,z,p)])
    //	  /(f.W[iix(x,y,z,p)] + Wold[iix(x,y,z,p)]);
    //	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
    // sort of E*exp(-2.0*s)
    
    
    // Top of p90 ch 3
	f.E[x][y][z] = f.E[x][y][z]*pow(Wold[x][y][z]/f.W[x][y][z],
					f.kappa[x][y][z]-1.0);



      }
    }
  }



  //  halo_field(f.W, p);

  halo_field(Wfacex, p);
  halo_field(Wfacey, p);
  halo_field(Wfacez, p);

  // Section 3.4.4 -- pressure terms
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	qx = (f.V[0][x+1][y][z] - f.V[0][x][y][z])/p.dx;
	qy = (f.V[1][x][y+1][z] - f.V[1][x][y][z])/p.dx;
	qz = (f.V[2][x][y][((z+1)%p.Lz)] - f.V[2][x][y][z])/p.dx;
	
	divv = qx + qy + qz;

		
	// s = 0.5*(f.kappa[x] - 1.0)*divv*p.dt/2.0;
	// f.E[x] = f.E[x]*(1-s)/(1+s);
	

	// Is it 1.0 or 0.5? -- not the major problem I think
	// also note that dx=1.0 still get problems
	// so it is not some factor of dx somewhere
	// nor is it this particular expression:
	f.E[x][y][z] = f.E[x][y][z]
	  *exp(-1.0*(f.kappa[x][y][z] - 1.0)*divv*p.dt);

      }
    }
  }



  //  halo_field(f.E, p);



  // Last bit of 3.4.6
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	qx = (f.V[0][x][y][z] + f.V[0][x+1][y][z])
	  *(Wfacex[x+1][y][z] - Wfacex[x][y][z])*p.dt/(2.0*p.dx);
	
	qy = (f.V[1][x][y][z] + f.V[1][x][y+1][z])
	  *(Wfacey[x][y+1][z] - Wfacey[x][y][z])*p.dt/(2.0*p.dx);
	
	qz = (f.V[2][x][y][z] + f.V[2][x][y][((z+1)%p.Lz)])
	  *(Wfacez[x][y][((z+1)%p.Lz)] - Wfacez[x][y][z])*p.dt/(2.0*p.dx);
	
	gradv = (qx+qy+qz)/f.W[x][y][z];
	
	//	s = (f.kappa[x] - 1.0)*gradv/2.0;
	//	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
	f.E[x][y][z] = f.E[x][y][z]*exp(-1.0*(f.kappa[x][y][z]-1.0)*gradv);
      }
    }      
  }



  halo_field(f.E, p);
  
  
  
  // Clean up memory. Surely we don't need all these temporary arrays?

  
  free_field(p, Vdmid);
  free_field(p, Wold);

  free_field(p, phiav);

  /*
  free_field(p, dxphi[0]);
  free_field(p, dxphi[1]);
  free_field(p, dxphi[2]);
  free(dxphi);
  */

  free_vector(p, dxphi);

  free_field(p, Wfacex);
  free_field(p, Wfacey);
  free_field(p, Wfacez);



}

