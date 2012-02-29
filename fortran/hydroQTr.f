*  hydroQTr.f  30.5.95  9.1.96   
      
      Program HydroQTr

*   **  Special relativistic spherical hydro + a scalar field phi with
*       potential V(phi,T) coupled to it.

*   --  Transport = Donor Cell  14.8.95

*   --  Log of changes
* 9.1.96   Use reflective core at origin (Ncore) to fix stability
* 2.11.95  List of slices to store (lprint = -1)
* 17.10.95 Tracer particles
* 11.10.95 pi = -vdefl*dphi/dt;  use ML's static wall profile
* 27.9.95  Raster images
* 15.9.95  init = 8:  read droplet profile
* 95/8/21  V0
* 95/8/18  init = 5:  droplet    8/22:  Tq, Th
* 95/8/15  init = 4:  gcen
* 93/5/18  Added init=3: alpha, gamma, lambda params
C 93/5/22  Added option to use Gaussian initial bubble 
C           (parameters istdbu (= 0 or 1), cstre, Rtensp)
C 93/5/26  T(eps,phi) is now solved analytically (FindTa instead of FindT)
C           => program is now 11.3 times as fast as it used to be
C 93/5/26  vwall and vshock are now calculated from the movement of xwall or 
C           xshock during Delta t = dt*sprint (WaVeDt instead of WallVelo)
C 93/5/29  Added option to use one-dimensional extremum bubble  (istdbu=2)
C           (best to use slightly oversized bu: cstre, Rtens = 1.01 ... 1.05)

*   **  Main variables

      integer M, N, bc, j, Ncycle, sprint, lprint, cycle, init, nres1,
     .        nres2, jV, jQ, mwall, mshock, rprint, xcycle, test, Ntb,
     .        r2first, r2last, r2left, r2right, Nwprf, Nprint, rres,
     .        Ntr, Ntrmax, trprint, Nlist, Mlist, sprf, sprl, ssprint,
     .        trprf, trprl, wprint, Ng, Ncore
      parameter (M = 65539, NTb = 1001, Nwprf = 2000, Ntrmax = 101,
     .   Mlist = 30)
      integer slclist(Mlist)

      real dt, dx, Lgrid, C, Cav, Courant, mD, kappa(M),
     .     fi(M), pi(M), D(M), E(M), Z(M), v(M), gb(M), p(M),
     .     Q(M), gbold(M), fiold(M), pifull(M),
c    .     pix(M), fix(M), Zx(M), Dxx(M), Ex(M), gbx(M),
     .     xe(M), xc(M),
     .     Tc, T(M), T0, alpha, gamma, lambda, gdeg, a, pii, 
     .     Lheat, sigma, lcorr, gq, gh, V0, aq, ah,
     .     Tplus, ThTb(NTb), ehTb(NTb), fir(2*M), er(2*M), vr(2*M), 
     .     pr(2*M), Tr(2*M), ehscale,
     .     xw(Nwprf), fiw(Nwprf), Tw(Nwprf), vw(Nwprf), vdefl,
     .     xtr(Ntrmax), Ltr, dtr,
*   --  parameters for initial data
     .     rE, xxwall, Tbar, Tq, Th,
*   --  workspace for subroutines
     .     Vdold(M), Delta(M), F(M), dxfi(M), gbv(M),
     .     fiav(M), Vdmid(M), Vnew(M),
*   --  energy balance
     .     Etot, restE, kinE, kinphi, grdphi, VE, QoverE,
*   --  wall velocity
     .     vwall, vshock, xwall, xshock, vwallfi, xwallfi, p1, p2, p3,
*   --  for raster plots
     .      emin, emax, Tmin, Tmax, fimin, fimax, vmin, vmax

C      Parameters for the initial bubble:
C      istdbu = 
C        0: Gaussian RLapl-sized bubble (cstre, Rtensp are dummy params)
C        1: Gaussian customized bubble 
C           (cstre = relat. central strength, Rtensp/rMboso = abs. width)
C        2: One-d extremum bubble (true extremum bu, if  cstre, Rtensp = 1.)
C        3: Tanh-bubble (cstre, Rtensp dummy).
C      If ".in"-file ends before the params istdbu, cstre, Rtensp 
C      are read in, then Tanh-bubble is used.
     
      real cstre, Rtensp
      integer istdbu 
                
      character*66 name

*     M:    maximum number of grid points
*     N:    actual number of grid points (physical = N-2)
*     bc:   boundary conditions    0 = periodic,  1 = reflective
*     rprint:  frequency of output for raster images
*     sprint:  frequency of output for t plots
*     lprint:  frequency of output for x plots
*     wprint:  frequency of output for wall motion
*     init:    type of initial data   1 = shock tube
*                 2 = V(phi,T) wall, with L, l, sigma params (gdeg = g_q)
*                 3 = V(phi,T) wall, with alpha, gamma, lambda params 
*                 4 = V(phi,T) wall, with L, l, sigma params (gdeg = gcen)

*     Lgrid:   physical length of grid = (N-2)*dx
*     C:       strength of field-fluid interaction
*     Cav:     artificial viscosity parameter
*     mD:      dynamical weight of D-density
*     kappa:   'equation of state'

*   --  The dynamical variables:
*     fi    the scalar field (order parameter, Higgs field)   at xc
*     pi    its conjugate momentum partial_t fi                  xc
*     D     dust density in lab frame                            xc
*     E     energy density in lab frame                          xc
*     Z                                                          xe
*     v     fluid flow velocity                                  xe
*     gb    the boost factor                                     xc
*     p     pressure (in fluid rest frame)                       xc
*     Q                                                          xc



*  **  Main program

      read(5,11) name
      write(6,11) name
 11   format (a66)
      
      xcycle = -1
      test = 0
      V0 = 0.0
      read(5,*) Ng, bc, Nprint, nres1, nres2, rres
      write(6,*) ' Ng ', Ng, ' bc ', bc, ' Nprint ', Nprint,
     .   ' nres ', nres1, nres2, ' rres ', rres
      N = Ng+2
      if (Nprint.eq.Ng) Nprint = N
      read(5,*) Ncycle, rprint, sprint, lprint, wprint
      write(6,*) ' Ncycle ', Ncycle,
     .   ' rprint ',  rprint,
     .   ' sprint ',  sprint, ' lprint ', lprint, ' wprint ', wprint
      if (rres.ge.0) then 
         read(5,*) emin, emax, Tmin, Tmax, fimin, fimax, vmin, vmax
      end if
      r2first = 2
      r2last = 1
      r2left = 1
      r2right = 1
      if (rres.eq.0) then
         read(5,*) r2first, r2last, r2left, r2right
         write(6,*) ' r2first ', r2first, ' r2last ', r2last,
     .      ' r2left ', r2left, ' r2right ', r2right
      end if
      if (lprint.eq.-1) then
         read(5,*) Nlist
         read(5,*) (slclist(i), i = 1, Nlist)
         write(6,*) ' Nlist ', Nlist, (slclist(i), i = 1, Nlist)
      end if
*   --  Deploy tracer particles
      read(5,*) Ntr, trprint, trprf, trprl
      write(6,*) ' Ntr ', Ntr, ' trprint ', trprint,
     .   ' trprf ', trprf, ' trprl ', trprl
      if (Ntr.ge.1) then
         read(5,*) (xtr(j), j = 1, Ntr)
         write(6,*) (xtr(j), j = 1, Ntr)
      end if
c     dtr = Ltr/Ntr
c     do 13 j = 1, Ntr
c        xtr(j) = j*dtr
c13   continue
      read(5,*) C, Cav, Courant, mD
      write(6,*) ' C ', C, ' Cav ', Cav, ' Courant ', Courant,
     .   ' mD ', mD
      read(5,*) init
      write(6,*) ' init ', init
      if (init.eq.1) then
         read(5,*) rE
         write(6,*) ' rE ', rE
      else if (init.eq.2) then
         read(5,*) Lgrid, xxwall, Ncore, mwall, mshock
         write(6,*) ' Lgrid ', Lgrid, ' xxwall ', xxwall,
     .      ' Ncore ', Ncore, ' mwall ', mwall, ' mshock ', mshock
         N = Ng-Ncore+2
         Nprint = Nprint-Ncore
         read(5,*) Tbar, gdeg
         write(6,*) ' Tbar ', Tbar, ' gdeg ', gdeg
         read(5,*) Lheat, sigma, lcorr
         write(6,*) ' Lheat ', Lheat, ' sigma ', sigma,
     .       ' lcorr ', lcorr
      else if (init.eq.3) then
         read(5,*) Lgrid, xxwall, Ncore, mwall, mshock
         write(6,*) ' Lgrid ', Lgrid, ' xxwall ', xxwall,
     .      ' Ncore ', Ncore, ' mwall ', mwall, ' mshock ', mshock
         N = Ng-Ncore+2
         Nprint = Nprint-Ncore
         read(5,*) Tbar, gdeg
         write(6,*) ' Tbar ', Tbar, ' gdeg ', gdeg
         read(5,*) gamma, alpha, lambda
         write(6,*) ' gamma ', gamma, ' alpha ', alpha,
     .       ' lambda ', lambda
      else if (init.eq.4) then
         read(5,*) Lgrid, xxwall, mwall, mshock
         write(6,*) ' Lgrid ', Lgrid, ' xxwall ', xxwall,
     .      ' Ncore ', Ncore, ' mwall ', mwall, ' mshock ', mshock
         N = Ng-Ncore+2
         Nprint = Nprint-Ncore
         read(5,*) Tbar, gdeg
         write(6,*) ' Tbar ', Tbar, ' gcen ', gdeg
         read(5,*) Lheat, sigma, lcorr
         write(6,*) ' Lheat ', Lheat, ' sigma ', sigma,
     .       ' lcorr ', lcorr
      else if (init.eq.5) then
         read(5,*) Lgrid, xxwall, mwall, mshock
         write(6,*) ' Lgrid ', Lgrid, ' xxwall ', xxwall,
     .      ' mwall ', mwall, ' mshock ', mshock
         read(5,*) Tq, Th, gq, gh, V0
         write(6,*) ' Tq ', Tq, ' Th ', Th, ' gq ', gq, ' gh ', gh,
     .      ' V0 ', V0 
         read(5,*) sigma, lcorr
         write(6,*) ' sigma ', sigma, ' lcorr ', lcorr
      else if (init.eq.6) then
         read(5,*) Lgrid, xxwall, mwall, mshock
         write(6,*) ' Lgrid ', Lgrid, ' xxwall ', xxwall,
     .      ' mwall ', mwall, ' mshock ', mshock
         read(5,*) Tq, Th, gdeg, V0
         write(6,*) ' Tq ', Tq, ' Th ', Th, ' gdeg ', gdeg,
     .      ' V0 ', V0
         read(5,*) Lheat, sigma, lcorr
         write(6,*) ' Lheat ', Lheat, ' sigma ', sigma,
     .      ' lcorr ', lcorr
      else if (init.eq.8) then
         read(5,*) sprf, sprl
         write(6,*) ' sprf ', sprf, ' sprl ', sprl
         read(5,*) Lgrid, mwall, mshock, vdefl
         write(6,*) ' Lgrid ', Lgrid,
     .      ' mwall ', mwall, ' mshock ', mshock, ' vdefl ', vdefl
         read(5,*) gdeg, V0
         write(6,*) ' gdeg ', gdeg, ' V0 ', V0
         read(5,*) Lheat, sigma, lcorr
         write(6,*) ' Lheat ', Lheat, ' sigma ', sigma,
     .      ' lcorr ', lcorr
      end if
      
      read(5,*,end=175) istdbu, cstre, Rtensp
      goto 185
175   istdbu = 3 
185   continue
       

      call OpenFile (name, init, rres)
      write(7,112) Ncycle, lprint, Nlist, N, Nprint, nres1, nres2
 112  format (7i6)
      if (rres.ge.1) then
      write(11,112) Ncycle, rprint, N, rres, Ncycle/rprint+1,
     .    (N-1)/rres+1
      write(12,112) Ncycle, rprint, N, rres, Ncycle/rprint+1,
     .    (N-1)/rres+1
      write(13,112) Ncycle, rprint, N, rres, Ncycle/rprint+1,
     .    (N-1)/rres+1
      write(14,112) Ncycle, rprint, N, rres, Ncycle/rprint+1, 
     .   (N-1)/rres+1
      write(21,112) Ncycle, rprint, N, nres1, r2last-r2first+1, 
     .   r2right-r2left+1
      write(22,112) Ncycle, rprint, N, nres1, r2last-r2first+1,
     .   r2right-r2left+1
      write(23,112) Ncycle, rprint, N, nres1, r2last-r2first+1,
     .   r2right-r2left+1
      write(24,112) Ncycle, rprint, N, nres1, r2last-r2first+1,
     .   r2right-r2left+1
      end if
      
      Tc = 1.
      pii = 4*atan(1.)
      cycle = 0

      if (init.eq.1) then
         call ShckTube (rE, N,   dx, Lgrid, xe, xc, fi, pifull,
     .       D, E, Z, v, gb)
      else if (init.eq.2 .or. init.eq.3 .or. init.eq.4) then
        if (istdbu.eq.0) then
          write(6,*) 'Initial bubble: standard Gaussian;  istdbu ', 
     .                 istdbu
          call Gauss (Lgrid, xxwall, Tbar, Tc, pii, gdeg,
     .       Lheat, sigma, lcorr, N, init,
     .       T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb, istdbu, cstre, Rtensp)
        else if (istdbu.eq.1) then
          write(6,*) 'Initial bubble: customized Gaussian;  istdbu ',
     .               istdbu, '  cstre ', cstre, '  Rtensp ', Rtensp
          call Gauss (Lgrid, xxwall, Tbar, Tc, pii, gdeg,
     .       Lheat, sigma, lcorr, N, init,
     .       T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb, istdbu, cstre, Rtensp)
        else if (istdbu.eq.2) then
          write(6,*) 'Initial bubble: One-d;  istdbu ',
     .               istdbu, '  cstre ', cstre, '  Rtensp ', Rtensp
          call Onedbu (Lgrid, xxwall, Tbar, Tc, pii, gdeg,
     .       Lheat, sigma, lcorr, N, init,
     .       T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb, istdbu, cstre, Rtensp)
        else if (istdbu.eq.3) then
          write(6,*) 'Initial bubble: Tanh;  istdbu ', 
     .               istdbu
          call Wall (Lgrid, xxwall, Tbar, Tc, pii, gdeg,
     .       Lheat, sigma, lcorr, N, Ncore, init,
     .       T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb)
        else
           write(6,*) 'No such initial bubble!'
        end if         
      else if (init.eq.5) then
         write(6,*) 'Initial droplet: Tanh'
         aq = gq*pii**2/90
         ah = gh*pii**2/90
         Lheat = 4*(aq-ah)*Tc**4
         a = aq
         call Droplet (Lgrid, xxwall, Tq, Th, Tc,
     .       Lheat, sigma, lcorr, N, init,
     .       T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb)
      else if (init.eq.6) then
         write(6,*) 'Initial droplet: Tanh'
         a = gdeg*pii**2/90
         call Droplet (Lgrid, xxwall, Tq, Th, Tc,
     .       Lheat, sigma, lcorr, N, init,
     .       T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb)
      else if (init.eq.8) then
         write(6,*) 'Initial similarity droplet: Tanh'
c        aq = gq*pii**2/90
c        ah = gh*pii**2/90
c        Lheat = 4*(aq-ah)*Tc**4
         a = gdeg*pii**2/90
         do 19 j = 1, 2*N
            read(4,*) fir(j), vr(j), Tr(j), er(j), pr(j)
 19      continue
         do 195 j = 1, Nwprf
            read(3,*) xw(j), fiw(j), Tw(j), vw(j)
 195     continue
         call SmVSDrop (Lgrid, fir, Tr, er, pr, vr, N,
     .      xw, fiw, Tw, vw, Nwprf, vdefl,
     .      Lheat, sigma, lcorr, a, Tc,
     .         alpha, gamma, lambda, T0, T, dx, xe, xc,
     .      fi, pifull, D, E, Z, v, gb)
      else
         write(6,*) 'No such initial data'
      end if

      dt = Courant*dx

      call BackStep (dt, dx, bc, N, fi, pifull, xe, xc, 
     .   alpha, gamma, lambda, T, T0,   pi)

      do 20 cycle = 0, Ncycle

*      --  Quantities that are needed from the previous timestep to
*          do a new one:  pi, fi, Z, v, D, E, gb

         if (cycle.eq.xcycle) test = 1
         call Tracer (dt, dx, bc, N, v, xe, cycle, Ntr,   xtr)
         if (cycle.ne.0) then

            if (cycle.eq.xcycle) write(6,*) 'Calling Field'
            call Field (dt, dx, C, bc, N, gb, v, xe, xc,   pi, fi,
     .            Vdold, alpha, gamma, lambda, T, T0, fiold, pifull)

            if (cycle.eq.xcycle) write(6,*) 'Calling Donorr'
            call Donorr  (dt, dx, bc, N, v, xe, xc,   D,   Delta, F)
            call Donorr  (dt, dx, bc, N, v, xe, xc,   E,   Delta, F)
            if (cycle.eq.xcycle) write(6,*) 'Calling DonorZr'
            call DonorZr (dt, dx, bc, N, v, xe, xc,   Z,   Delta, F)

            if (cycle.eq.xcycle) write(6,*) 'Calling EqState'
            call EqState (a, bc, N, E, gb, fi, alpha, gamma,
     .         lambda, V0, T0, test,   T,   p, kappa, Vnew)
            if (test.eq.2) then
               write(6,*) 'Stopping at cycle ', cycle
               stop
            end if

            if (cycle.eq.xcycle) write(6,*) 'Calling Hydro'
            call Hydro (dt, dx, kappa, C, Cav, mD, bc, N,
     .         fiold, fi, pi, p, D, xe, xc,
     .            alpha, gamma, lambda, T0, a,
     .            E, Z, v, gb, T,
     .            Q, gbold, dxfi, fiav, gbv, Vdmid)

C            do 12 j = 2, N-1
C               call FindT(E(j), gb(j), fi(j), a, alpha, gamma,
C     .            lambda, V0, T0, j,   T(j)) 
C 12         continue
            if (cycle.eq.xcycle) write(6,*) 'Calling FindTa'
            call FindTa( E, gb, fi, a, alpha, gamma,
     .         lambda, V0, T0, N, test,   T ) 
  
            if (cycle.eq.xcycle) write(6,*) 'Calling WrapEven'
            call WrapEven (bc, N,   T)
         end if

         if (mod(cycle,sprint).eq.0
     .      .or. cycle.ge.sprf .and. cycle.le.sprl) then
            if (cycle.gt.sprf .and. cycle.le.sprl) then
               ssprint = 1
            else
               ssprint = sprint
            end if
            call Ebalance (dx, N, kappa, E, gb, xe, xc,
     .         fi, pifull, Q, alpha, gamma, lambda, V0, T, T0,
     .            Etot, restE, kinE, kinphi, grdphi, VE, QoverE)
C            call WallVelo (xc, N, fi, Q, dt, cycle,
C     .         alpha, gamma, lambda, V0, T, T0,
C     .            vwall, vshock, jV, jQ, xwall, xshock) 
            call WaVeDt (xc, N, fi, Q, dt, cycle,
     .         alpha, gamma, lambda, V0, T, T0,
     .            vwall, vshock, jV, jQ, xwall, xshock, ssprint) 
         end if
            if (cycle.eq.xcycle) write(6,*) 'Calling Output'
         call Output (init, cycle, sprint, lprint, dt, dx, N, nres1,
     .      nres2, Nprint, xc, xe, fi, D, E, Q, Z, v, gb, kappa,
     .      alpha, gamma, lambda, V0, T, T0,
     .      Etot, restE, kinE, kinphi, grdphi, VE, QoverE, 
     .      vwall, vshock, jV, jQ, xwall, xshock, mwall, mshock,
     .      slclist, Nlist, sprf, sprl)
         if (mod(cycle,wprint).eq.0) then
            call WaVefi (xe, xc, N, gb, v, pi, fi, Q, E, kappa,
     .         dt, dx, cycle, C, alpha, gamma, lambda, V0, T, T0,
     .            vwallfi, jV, xwallfi, p1, p2, p3) 
            write(18,198) cycle*dt, xwallfi, vwallfi, p1, p2, p3
         end if
 198     format(6f14.6)
         call Trout (cycle, dt, trprint, trprf, trprl, Ntr, xtr)
         if (rres.ge.1)
     .   call Raster (cycle, rprint, r2first, r2last, r2left, r2right,
     .      N, rres, E, gb, T, fi, v,
     .      emin, emax, Tmin, Tmax, fimin, fimax, vmin, vmax)

 20   continue

      stop
      end
*     of Main program



      subroutine OpenFile (name, init, rres)
      integer init, rres, lblank
      character*66 name, fname(24)

      lblank = index (name, ' ')
      fname(3) = 'initial.dat'
      fname(4) = name(1:lblank-1)//':4'
      fname(7) = name(1:lblank-1)//'.7'
      fname(8) = name(1:lblank-1)//'.8'
      fname(9) = name(1:lblank-1)//'.9'
      fname(11) = name(1:lblank-1)//'.11'
      fname(12) = name(1:lblank-1)//'.12'
      fname(13) = name(1:lblank-1)//'.13'
      fname(14) = name(1:lblank-1)//'.14'
      fname(18) = name(1:lblank-1)//'.18'
      fname(19) = name(1:lblank-1)//'.19'
      fname(21) = name(1:lblank-1)//'.21'
      fname(22) = name(1:lblank-1)//'.22'
      fname(23) = name(1:lblank-1)//'.23'
      fname(24) = name(1:lblank-1)//'.24'
      open (unit = 3, file = fname(3), status = 'unknown')
      open (unit = 4, file = fname(4), status = 'unknown')
      open (unit = 7, file = fname(7), status = 'unknown')
      open (unit = 8, file = fname(8), status = 'unknown')
      open (unit = 9, file = fname(9), status = 'unknown')
      open (unit = 18, file = fname(18), status = 'unknown')
      open (unit = 19, file = fname(19), status = 'unknown')
      if (rres.ge.1) then
      open (unit = 11, file = fname(11), status = 'unknown')
      open (unit = 12, file = fname(12), status = 'unknown')
      open (unit = 13, file = fname(13), status = 'unknown')
      open (unit = 14, file = fname(14), status = 'unknown')
      open (unit = 21, file = fname(21), status = 'unknown')
      open (unit = 22, file = fname(22), status = 'unknown')
      open (unit = 23, file = fname(23), status = 'unknown')
      open (unit = 24, file = fname(24), status = 'unknown')
      write(11,12) name
      write(12,12) name
      write(13,12) name
      write(14,12) name
      write(21,12) name
      write(22,12) name
      write(23,12) name
      write(24,12) name
      end if
      write(6,12) name
      write(7,12) name
      write(8,12) name
      write(9,12) name
 12   format (a66)

      return
      end
*     of OpenFile



      subroutine ehTable (a, alpha, gamma, lambda, V0, T0, Tc, NTb,
     .      Tplus, ThTb, ehTb)
      integer Ntb, i
      real a, alpha, gamma, lambda, V0, T0, Tc, Tplus, ThTb(NTb),
     .     ehTb(NTb), dT, fimin, Vf, VTf

      write(6,*) ' ehTable '
      Tplus = sqrt(8*T0**2/(9*(T0/Tc)**2-1))
      dT = (Tplus-T0)/Ntb
      write(6,*) ' T0 ', T0, ' Tplus ', Tplus
      do 20 i = 1, Ntb
         T = T0 + (i-0.5)*dT
         ThTb(i) = T
         fimin = ( alpha*T
     .      + sqrt((alpha*T)**2-4*lambda*gamma*(T**2-T0**2)) )
     .      / (2*lambda)
         ehTb(i) = 3*a*T**4
     .      + Vf(alpha, gamma, lambda, V0, T, T0, fimin)
     .      - T*VTf(alpha, gamma, lambda, T, T0, fimin)
 20   continue
    
      return
      end
*     of ehTable


      REAL FUNCTION Theh(eh, ThTb, ehTb, NTb)
      INTEGER NTb, jlo
      REAL eh, ThTb(NTb), ehTb(NTb), Th, dTh
      call hunt (ehTb, NTb, eh, jlo)
      if (jlo.eq.0 .or. jlo.eq.NTb) then
         WRITE(6,*) ' Theh:  eh = ', eh, ' out of range,  jlo = ', jlo
         STOP
      end if
      if (jlo.eq.1) jlo = 2
      if (jlo.eq.NTb-1) jlo = NTb-2
      call polint (ehTb(jlo-1), ThTb(jlo-1), 4, eh, Th, dTh)
      Theh = Th
      RETURN
      END
*     off Theh



      subroutine ShckTube (rE, N,   dx, Lgrid, xe, xc, fi, pifull,
     .   D, E, Z, v, gb)
      integer N, j, Nhalf
      real rE, dx, Lgrid, xe(N), xc(N), fi(N), pifull(N),
     .     D(N), E(N), Z(N), v(N), gb(N), Er, El
*   --  Initial data for the relativistic shock tube

      Lgrid = 1.
      dx = Lgrid/(N-2)
      do 11 j = 1, N
         xe(j) = (j-1)*dx
         xc(j) = (j-1.5)*dx
 11   continue
      Nhalf = N/2
      Er = 1./(1.+rE)
      El = rE*Er
      do 12 j = 1, Nhalf
         E(j) = El
 12   continue
      do 14 j = Nhalf+1, N
         E(j) = Er
 14   continue
      do 16 j = 1, N
         fi(j) = 0.
         pifull(j) = 0.
         D(j) = 1.
         Z(j) = 0.
         v(j) = 0.
         gb(j) = 1.
 16   continue

      return
      end
*     of ShckTube



      subroutine Wall (Lgrid, xxwall, Tbar, Tc, pii, gdeg,
     .   Lheat, sigma, lcorr, N, Ncore, init,
     .      T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .   fi, pifull, D, E, Z, v, gb)
      integer N, j, init
      real Lgrid, xxwall, Tbar, Tc, gdeg, Lheat, sigma, lcorr,
     .     Tconst, T(N), T0, alpha, gamma, lambda, V0, a, dx,
     .     xe(N), xc(N), fi(N), pifull(N), D(N), E(N), Z(N),
     .     v(N), gb(N), fimin, pii, Vf, VTf, Lcore
*   --  Initial data for a phase transition wall
C       phi-field is initially a tanh-function, which has been cut off
    
      Tconst = Tbar*Tc
      if (init.eq.2 .or. init.eq.4) then
         T0 = Tc/sqrt(1+6*sigma/Lheat/lcorr)
         alpha = 1/sqrt(2*sigma*lcorr**5*Tc**2/3)
         gamma = (Lheat+6*sigma/lcorr)/(6*sigma*lcorr*Tc**2)
         lambda = 1/(3*sigma*lcorr**3)
      else if (init.eq.3) then
         T0 = sqrt(1.-2*alpha**2/9/gamma/lambda)*Tc
         Lheat = (4./9.)*(alpha**2*gamma/lambda**2)*(T0*Tc)**2
         lcorr = 1./sqrt(2./9./lambda)/(alpha*Tc)
         sigma = 2*sqrt(2.)*(alpha*Tc)**3/(81*lambda**2.5)
      else
         stop
      end if
      if (init.eq.4) then
         a = gdeg*pii**2/90 + Lheat/8/Tc**4
      else
         a = gdeg*pii**2/90
      end if
      fimin = ( alpha*Tconst
     .   + sqrt((alpha*Tconst)**2-4*lambda*gamma*(Tconst**2-T0**2)) )
     .   / (2*lambda)
      write(6,12) Tconst/Tc, T0/Tc, Tc
      write(6,14) a
      write(6,16) alpha, gamma, lambda
      write(6,17) fimin
      write(6,175) Lheat, sigma, lcorr
 12   format (' T =', f9.6, 'Tc   T0 =', f9.6, 'Tc   Tc =', f9.6)
 14   format (' a =', f9.6)
 16   format (' alpha =', f9.6, '   gamma =', f9.6,
     .   '   lambda =', f9.6)
 17   format (' phi min =', f9.6)
 175  format (' Lheat =', f9.5, '   sigma =', f9.6, '   lcorr =', f9.6)

      dx = Lgrid/(N+Ncore-2)
      Lcore = Ncore*dx
      write(6,177) Lcore
 177  format (' Lcore = ', f9.5)
      do 18 j = 1, N
         xe(j) = Lcore+(j-1)*dx
         xc(j) = Lcore+(j-1.5)*dx
         fi(j) = 0.5*fimin*(1 -
     .      tanh(sqrt(lambda/2)*fimin/2*(xc(j)-xxwall)))
         pifull(j) = 0.
         T(j) = Tconst
         D(j) = 1.
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .      - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
         Z(j) = 0.
         v(j) = 0.
         gb(j) = 1.
 18   continue

      return
      end
*     of Wall



      subroutine Gauss (Lgrid, xxwall, Tbar, Tc, pii, gdeg,
     .   Lheat, sigma, lcorr, N, init,
     .      T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .   fi, pifull, D, E, Z, v, gb, istdbu, cstre, Rtensp)
      integer N, j, init
      real Lgrid, xxwall, Tbar, Tc, gdeg, Lheat, sigma, lcorr,
     .     Tconst, T(N), T0, alpha, gamma, lambda, V0, a, dx,
     .     xe(N), xc(N), fi(N), pifull(N), D(N), E(N), Z(N),
     .     v(N), gb(N), fimin, pii, Vf, VTf
      real cstre, Rtensp, rMboso, cstrab, Rtenab, sigmlo, Rlapab
      integer istdbu
*   --  Initial data for a phase transition wall
C       A Gaussian bubble for the phi-field
    
      Tconst = Tbar*Tc
      if (init.eq.2) then
         T0 = Tc/sqrt(1+6*sigma/Lheat/lcorr)
         alpha = 1/sqrt(2*sigma*lcorr**5*Tc**2/3)
         gamma = (Lheat+6*sigma/lcorr)/(6*sigma*lcorr*Tc**2)
         lambda = 1/(3*sigma*lcorr**3)
         
         sigmlo = sigma
      else if (init.eq.3) then
         T0 = sqrt(1.-2*alpha**2/9/gamma/lambda)*Tc

         sigmlo = 2.*sqrt(2.)/81. * alpha**3 /lambda**2.5 *Tc**3
      else
         stop
      end if
      a = gdeg*pii**2/90
      fimin = ( alpha*Tconst
     .   + sqrt((alpha*Tconst)**2-4*lambda*gamma*(Tconst**2-T0**2)) )
     .   / (2*lambda)
      write(6,12) Tconst/Tc, T0/Tc, Tc
      write(6,14) a
      write(6,16) alpha, gamma, lambda
      write(6,17) fimin
 12   format (' T =', f9.6, 'Tc   T0 =', f9.6, 'Tc   Tc =', f9.6)
 14   format (' a =', f9.6)          
 16   format (' alpha =', f9.6, '   gamma =', f9.6,
     .   '   lambda =', f9.6)
 17   format (' phi min =', f9.6)

C      rt = 3.37923 / 0.0576937
C      stren = 0.93295 * 0.969143
      if (istdbu.eq.0) then
         cstrab = 1. * fimin
         Rlapab = 2.*sigmlo
     .      /(-Vf(alpha,gamma,lambda,V0,Tconst,T0,fimin))
         Rtenab = Rlapab
C         write(6,'(3f12.7)') sigmlo, 
C     .      Vf(alpha,gamma,lambda,V0,Tconst,T0,fimin), Rlapab
      else if (istdbu.eq.1) then
         rMboso = sqrt( gamma*(Tconst**2 - T0**2) )
         cstrab = cstre * fimin
         Rtenab = Rtensp / rMboso
      else
         write(6,*) 'Gauss: Error;  istdbu = ', istdbu
      end if
      write(6,195) cstrab, Rtenab
195   format( ' cstrab = ', f12.7, '   Rtenab =', f12.7 )

      dx = Lgrid/(N-2)
      do 18 j = 1, N
         xe(j) = -xxwall+(j-1)*dx
         xc(j) = -xxwall+(j-1.5)*dx
C         fi(j) = 0.5*fimin*(1-tanh(sqrt(lambda/2)*fimin/2*xc(j)))
C         fi(j) = stren * exp( - ( xc(j)+xxwall+0.5*dx )**2 /2./rt**2 )
         fi(j) = cstrab * 
     .            exp( - ( xc(j)+xxwall+0.5*dx )**2 /2./Rtenab**2 )
         pifull(j) = 0.
         T(j) = Tconst
         D(j) = 1.
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .      - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
         Z(j) = 0.
         v(j) = 0.
         gb(j) = 1.
 18   continue

      return
      end
*     of Gauss



      subroutine Onedbu (Lgrid, xxwall, Tbar, Tc, pii, gdeg,
     .   Lheat, sigma, lcorr, N, init,
     .      T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .   fi, pifull, D, E, Z, v, gb, istdbu, cstre, Rtensp)
      integer N, j, init
      real Lgrid, xxwall, Tbar, Tc, gdeg, Lheat, sigma, lcorr,
     .     Tconst, T(N), T0, alpha, gamma, lambda, V0, a, dx,
     .     xe(N), xc(N), fi(N), pifull(N), D(N), E(N), Z(N),
     .     v(N), gb(N), fimin, pii, psibar, Vf, VTf
      real cstre, Rtensp, rMboso, deofV, lbar
      integer istdbu
*   --  Initial data for a phase transition wall
C       An one-d extremum bubble for the phi-field.
C       cstre = relative central strength, Rtensp = relative radius 
C       (true extremum bubble, if  cstre, Rtensp = 1.)
    
      Tconst = Tbar*Tc
      if (init.eq.2) then
         T0 = Tc/sqrt(1+6*sigma/Lheat/lcorr)
         alpha = 1/sqrt(2*sigma*lcorr**5*Tc**2/3)
         gamma = (Lheat+6*sigma/lcorr)/(6*sigma*lcorr*Tc**2)
         lambda = 1/(3*sigma*lcorr**3)
      else if (init.eq.3) then
         T0 = sqrt(1.-2*alpha**2/9/gamma/lambda)*Tc
      else
         stop
      end if
      a = gdeg*pii**2/90
      fimin = ( alpha*Tconst
     .   + sqrt((alpha*Tconst)**2-4*lambda*gamma*(Tconst**2-T0**2)) )
     .   / (2*lambda)
      write(6,12) Tconst/Tc, T0/Tc, Tc
      write(6,14) a
      write(6,16) alpha, gamma, lambda
      write(6,17) fimin
 12   format (' T =', f9.6, 'Tc   T0 =', f9.6, 'Tc   Tc =', f9.6)
 14   format (' a =', f9.6)
 16   format (' alpha =', f9.6, '   gamma =', f9.6,
     .   '   lambda =', f9.6)
 17   format (' phi min =', f9.6)

      rMboso = sqrt( gamma*(Tconst**2 - T0**2) )
      deofV = alpha*Tconst
      lbar = 9./2. * lambda*rMboso**2 / deofV**2
      write(6,'('' lbar = '',f12.7)') lbar


      dx = Lgrid/(N-2)
      do 18 j = 1, N
         xe(j) = -xxwall+(j-1)*dx
         xc(j) = -xxwall+(j-1.5)*dx
         fi(j) = cstre *3.*rMboso**2/2./deofV 
     .    * psibar( rMboso*( xc(j)+xxwall+0.5*dx )/Rtensp , lbar)
         pifull(j) = 0.
         T(j) = Tconst
         D(j) = 1.
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .      - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
         Z(j) = 0.
         v(j) = 0.
         gb(j) = 1.
 18   continue

      return
      end
*     of Onedbu



      subroutine Droplet (Lgrid, xxwall, Tq, Th, Tc,
     .   Lheat, sigma, lcorr, N, init,
     .      T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .   fi, pifull, D, E, Z, v, gb)
      integer N, j, init
      real Lgrid, xxwall, Tq, Th, Tc, gdeg, Lheat, sigma, lcorr,
     .     Tconst, Tqc, T(N), T0, alpha, gamma, lambda, V0, a, dx,
     .     xe(N), xc(N), fi(N), pifull(N), D(N), E(N), Z(N),
     .     v(N), gb(N), fimin, Vf, VTf
*   --  Initial data for a phase transition wall: Droplqt
C       phi-field is initially a tanh-function, which has been cut off

      Tconst = Th*Tc
      Tqc = Tq*Tc
      T0 = Tc/sqrt(1+6*sigma/Lheat/lcorr)
      alpha = 1/sqrt(2*sigma*lcorr**5*Tc**2/3)
      gamma = (Lheat+6*sigma/lcorr)/(6*sigma*lcorr*Tc**2)
      lambda = 1/(3*sigma*lcorr**3)
      fimin = ( alpha*Tconst
     .   + sqrt((alpha*Tconst)**2-4*lambda*gamma*(Tconst**2-T0**2)) )
     .   / (2*lambda)
      write(6,12) Tconst/Tc, T0/Tc, Tc
      write(6,14) a
      write(6,16) alpha, gamma, lambda
      write(6,17) fimin
 12   format (' T =', f9.6, 'Tc   T0 =', f9.6, 'Tc   Tc =', f9.6)
 14   format (' a =', f9.6)
 16   format (' alpha =', f9.6, '   gamma =', f9.6,
     .   '   lambda =', f9.6)
 17   format (' phi min =', f9.6)

      dx = Lgrid/(N-2)
      do 18 j = 1, N
         xe(j) = (j-1)*dx
         xc(j) = (j-1.5)*dx
         fi(j) = 0.5*fimin*(1 -
     .      tanh(sqrt(lambda/2)*fimin/2*(xxwall-xc(j))))
         pifull(j) = 0.
         T(j) = 0.5*(Tconst+Tqc) + 0.5*(Tqc-Tconst)
     .      * tanh(sqrt(lambda/2)*fimin/2*(xxwall-xc(j)))
         D(j) = 1.
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .      - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
         Z(j) = 0.
         v(j) = 0.
         gb(j) = 1.
 18   continue

      return
      end
*     of Droplet


      SUBROUTINE SmVSDrop (Lgrid, fir, Tr, er, pr, vr, N,
     .   xw, fiw, Tw, vw, Nwprf, vdefl,
     .   Lheat, sigma, lcorr, a, Tc,
     .       alpha, gamma, lambda, T0, T, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb)
      INTEGER N, j, jwall, jw, Nwprf
      REAL Lgrid, fir(2*N), Tr(2*N), er(2*N), pr(2*N), vr(2*N),
     .     xw(Nwprf), fiw(Nwprf), Tw(Nwprf), vw(Nwprf),
     .     Lheat, sigma, lcorr, a, Tc, alpha, gamma, lambda, T0,
     .     T(N), xe(N), xc(N),
     .     fi(N), pifull(N), D(N), E(N), Z(N), v(N), gb(N), vc,
     .     Tf, Tq, fiminf, vvf, flf, frg, Tlf, Trg, vlf, vrg,
     .     fix, Tx, vx, vdefl, w
      
      T0 = Tc/sqrt(1+6*sigma/Lheat/lcorr)
      alpha = 1/sqrt(2*sigma*lcorr**5*Tc**2/3)
      gamma = (Lheat+6*sigma/lcorr)/(6*sigma*lcorr*Tc**2)
      lambda = 1/(3*sigma*lcorr**3)
      write(6,11) T0/Tc, Tc
      write(6,12) a
      write(6,13) alpha, gamma, lambda
 11   format (' T0 =', f9.6, 'Tc   Tc =', f9.6)
 12   format (' a =', f9.6)
 13   format (' alpha =', f9.6, '   gamma =', f9.6,
     .   '   lambda =', f9.6)
      dx = Lgrid/(N-2)
      do 14 j = 1, N
         xe(j) = (j-1)*dx
         xc(j) = (j-1.5)*dx
 14   continue

*   --  Find the phase boundary from input profile:  xxwall, vvf, Tfc, fiminf
      do 16 j = 2, 2*N
         if (fir(j).ne.fir(j-1)) then
            if (mod(j,2).eq.0) then
               xxwall = 0.5*(xc(j/2)+xe(j/2))
            else
               xxwall = 0.5*(xe(j/2)+xc(j/2+1))
            end if
            jwall = j
            vvf = vr(j)
            Tq = Tr(j-1)
            Tf = Tr(j)
            fiminf = fir(j)
            go to 17
         end if
 16   continue
      write(6,*) ' SmVSDrop:  did not find phase boundary '
      STOP
 17   continue
      write(6,*) ' SmVSDrop:  Phase boundary at xxwall = ', xxwall,
     .   ' vvf ', vvf, ' Tf ', Tf, ' fiminf ', fiminf

c     GO TO 27
*   --  Find the phase boundary from wall profile:
      do 175 j = 1, Nwprf
         if(xw(j-1)*xw(j).lt.0.) then
            jw = j
            go to 176
         end if
 175  continue
      write(6,*) ' SimSDrop:  did not find center of wall profile '
      STOP
 176  continue
      write(6,*) ' Wall ', xw(jw-1), xw(jw), jw
      flf = fiw(jw-100)
      frg = fiw(jw+99)
      Tlf = Tw(jw-100)
      Trg = Tw(jw+99)
      vlf = vw(jw-100)
      vrg = vw(jw+99)
      write(6,*) flf, frg, fiminf
      write(6,*) Tlf, Trg, Tq, Tf
      write(6,*) vlf, vrg, vvf
      write(6,*) xw(jw+1)-xw(jw-1), xe(j/2+1)-xe(j/2)

*   --  Apply profile to similarity solution:
c     do 179 j = -120, 120
c        write(6,1785) xw(jw+j), fir(jwall+j), Tr(jwall+j), vr(jwall+j)
c1785    format (4f12.6)
c179  continue
      do 178 j = -100, 99
         fix = fiminf
         Tx = Tf-Tq
         vx = vvf
         if (j.ge.0) then
            fix = fir(jwall+j)
            Tx = Tr(jwall+j)-Tq
            vx = vr(jwall+j)
         end if
         fir(jwall+j) = fix/frg*fiw(jw+j)
         Tr(jwall+j) = Tq+Tx/(Trg-Tlf)*(Tw(jw+j)-Tlf)
         vr(jwall+j) = vx/vrg*vw(jw+j)
 178  continue

      pifull(1) = 0.
      do 22 j = 1, N
*      --  edge quantities v, Z
         fi(j) = fir(2*j)
         v(j) = vr(2*j)
         gb(j) = 1./sqrt(1.-v(j)**2)
         T(j) = Tr(2*j)
         Z(j) = (4*a*T(j)**4
     .      - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j)))
     .      *gb(j)**2*v(j)
*      --  cell quantities  fi, T, gb, E
         fi(j) = fir(2*j-1)
         if (j.ge.2)
     .      pifull(j) = -vdefl*(fir(2*j)-fir(2*j-2))/(xe(j)-xe(j-1))
         vc = vr(2*j-1)
         gb(j) = 1./sqrt(1.-vc**2)
         T(j) = Tr(2*j-1)
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
         w = 4*a*T(j)**4 - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
         E(j) = E(j)*gb(j)
         D(j) = w/T(j)
 22   continue
          
      RETURN

 27   CONTINUE
*   --  Tanh-profile
      pifull(1) = 0.
      do 28 j = 1, N
         D(j) = 1.
c        pifull(j) = 0.
*      --  edge quantities v, Z
         if (xe(j).lt.xxwall) then
            fi(j) = 0.5*fiminf*(1 +
     .         tanh(sqrt(lambda/2)*fiminf/2*(xe(j)-xxwall)))
            v(j) = 0.5*vvf*(1+
     .         tanh(sqrt(lambda/2)*fiminf/2*(xe(j)-xxwall)))
            gb(j) = 1./sqrt(1.-v(j)**2)
            T(j) = 0.5*(Tf+Tq) + 0.5*(Tf-Tq)
     .         * tanh(sqrt(lambda/2)*fiminf/2*(xe(j)-xxwall))
            Z(j) = (4*a*T(j)**4
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j)))
     .         *gb(j)**2*v(j)
         else
            fi(j) = 0.5*fir(2*j)*(1 +
     .         tanh(sqrt(lambda/2)*fiminf/2*(xe(j)-xxwall)))
            v(j) = 0.5*vr(2*j)*(1+
     .         tanh(sqrt(lambda/2)*fiminf/2*(xe(j)-xxwall)))
            gb(j) = 1./sqrt(1.-v(j)**2)
            T(j) = 0.5*(Tr(2*j)+Tq) + 0.5*(Tr(2*j)-Tq)
     .         * tanh(sqrt(lambda/2)*fiminf/2*(xe(j)-xxwall))
            Z(j) = (4*a*T(j)**4
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j)))
     .         *gb(j)**2*v(j)
         end if
*      --  cell quantities  fi, T, gb, E
         if (xc(j).lt.xxwall) then
            if (j.ge.2)
     .        pifull(j) = -vdefl*(fi(j)-fi(j-1))/(xe(j)-xe(j-1))
            fi(j) = 0.5*fiminf*(1 +
     .         tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall)))
            vc = 0.5*vvf*(1+
     .         tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall)))
            gb(j) = 1./sqrt(1.-vc**2)
            T(j) = 0.5*(Tf+Tq) + 0.5*(Tf-Tq)
     .         * tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall))
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
            E(j) = E(j)*gb(j)
         else
            fi(j) = 0.5*fir(2*j-1)*(1 +
     .         tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall)))
            vc = 0.5*vr(2*j-1)*(1+
     .         tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall)))
            gb(j) = 1./sqrt(1.-vc**2)
            T(j) = 0.5*(Tr(2*j-1)+Tq) + 0.5*(Tr(2*j-1)-Tq)
     .         * tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall))
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
            E(j) = E(j)*gb(j)
         end if

 28   continue

      RETURN
      END
*     of SmVSDrop
         


      
      SUBROUTINE SimSDrop (Lgrid, Tq, Tc, fir, er, vr, ThTb, ehTb,
     .       NTb, ehscale, Lheat, sigma, lcorr, N,
     .       T, T0, alpha, gamma, lambda, V0, a, dx, xe, xc,
     .       fi, pifull, D, E, Z, v, gb)
      INTEGER N, j, NTb
      REAL Lgrid, Tq, Tc, fir(2*N), er(2*N), vr(2*N),
     .       ThTb(NTb), ehTb(NTb), ehscale, 
     .       Lheat, sigma, lcorr, 
     .       T(N), T0, alpha, gamma, lambda, V0, a, dx, xe(N), xc(N),
     .       fi(N), pifull(N), D(N), E(N), Z(N), v(N), gb(N), Theh,
     .       Tplus, fiminf, fimin, Tqc, Tfc, vvf, vc
      EXTERNAL Theh

      Tqc = Tq*Tc
      T0 = Tc/sqrt(1+6*sigma/Lheat/lcorr)
      alpha = 1/sqrt(2*sigma*lcorr**5*Tc**2/3)
      gamma = (Lheat+6*sigma/lcorr)/(6*sigma*lcorr*Tc**2)
      lambda = 1/(3*sigma*lcorr**3)

      call ehTable (a, alpha, gamma, lambda, V0, T0, Tc, NTb,
     .      Tplus, ThTb, ehTb)

      dx = Lgrid/(N-2)
      do 14 j = 1, N
         xe(j) = (j-1)*dx
         xc(j) = (j-1.5)*dx
 14   continue

*   --  Find the phase boundary from input profile:  xxwall, vvf, Tfc, fiminf
      do 16 j = 2, 2*N
         if (fir(j).ne.fir(j-1)) then
            if (mod(j,2).eq.0) then
               xxwall = 0.5*(xc(j/2)+xe(j/2))
            else
               xxwall = 0.5*(xe(j/2)+xc(j/2+1))
            end if
            vvf = vr(j)
           Tfc = Theh(ehscale*er(j)/sqrt(1.-vr(j)**2), ThTb, ehTb, NTb)
            fiminf = ( alpha*Tfc
     .         + sqrt((alpha*Tfc)**2-4*lambda*gamma*(Tfc**2-T0**2)) )
     .         / (2*lambda)
            go to 17
         end if
 16   continue
      write(6,*) ' SimSDrop:  didn not find phase boundary '
      STOP
 17   continue
      write(6,*) ' SimSDrop:  Phase boundary at xxwall = ', xxwall,
     .   ' vvf ', vvf, ' Tfc ', Tfc, ' fiminf ', fiminf

      do 18 j = 1, N
         D(j) = 1.
         pifull(j) = 0.
         if (fir(2*j-1).eq.0.) then
*         --  Inside the quark droplet
c           v(j) = 0.
c           gb(j) = 1.
c           E(j) = er(2*j-1)
c           T(j) = Tqc
            fi(j) = 0.5*fiminf*(1 +
     .         tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall)))
            T(j) = 0.5*(Tfc+Tqc) + 0.5*(Tfc-Tqc)
     .         * tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall))
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
            v(j) = 0.5*vvf*(1+
     .         tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall)))
            gb(j) = 1./sqrt(1.-v(j)**2)
            Z(j) = (4*a*T(j)**4
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j)))
     .         *gb(j)**2*v(j)
         else
*         --  Outside the quark droplet: v(j), Z(j)
*         --    First find the reference values  v(j), fimin, T(j)
            v(j) = vr(2*j)
            gb(j) = 1./sqrt(1.-v(j)**2)
            E(j) = ehscale*er(2*j)*gb(j)
            T(j) = Theh(ehscale*er(2*j), ThTb, ehTb, NTb)
            fimin = ( alpha*T(j)
     .         + sqrt((alpha*T(j))**2-4*lambda*gamma*(T(j)**2-T0**2)) )
     .         / (2*lambda)
*         --    Then impose the tanh profile
            fi(j) = 0.5*fimin*(1 -
     .         tanh(sqrt(lambda/2)*fimin/2*(xc(j)-xxwall)))
            T(j) = 0.5*(T(j)+Tqc) + 0.5*(T(j)-Tqc)
     .         * tanh(sqrt(lambda/2)*fimin/2*(xc(j)-xxwall))
            v(j) = 0.5*v(j)*(1+
     .         tanh(sqrt(lambda/2)*fiminf/2*(xc(j)-xxwall)))
            gb(j) = 1./sqrt(1.-v(j)**2)
            Z(j) = (4*a*T(j)**4
     .              - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j)))
     .             *gb(j)**2*v(j)

*         --  Outside the quark droplet: fi(j), T(j), gb(j), E(j) 
*         --    First find the reference values  fimin, T(j), vc
            vc = vr(2*j-1)
            gb(j) = 1./sqrt(1.-vc**2)
            E(j) = ehscale*er(2*j-1)*gb(j)
            T(j) = Theh(ehscale*er(2*j-1), ThTb, ehTb, NTb)
            fimin = ( alpha*T(j)
     .         + sqrt((alpha*T(j))**2-4*lambda*gamma*(T(j)**2-T0**2)) )
     .         / (2*lambda)
*         --    Then impose the tanh profile
            fi(j) = 0.5*fimin*(1 +
     .         tanh(sqrt(lambda/2)*fimin/2*(xc(j)-xxwall)))
            T(j) = 0.5*(T(j)+Tqc) + 0.5*(T(j)-Tqc)
     .         * tanh(sqrt(lambda/2)*fimin/2*(xc(j)-xxwall))
         E(j) = 3*a*T(j)**4 + Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
     .         - T(j)*VTf(alpha,gamma,lambda,T(j),T0,fi(j))
            vc =  0.5*vc*(1+
     .         tanh(sqrt(lambda/2)*fimin/2*(xc(j)-xxwall))) 
            gb(j) = 1./sqrt(1.-vc**2)
            E(j) = E(j)*gb(j)
         end if
 18   continue
      RETURN
      END
*     of SimSDrop



      subroutine BackStep (dt, dx, bc, N, fi, pifull, xe, xc,
     .   alpha, gamma, lambda, T, T0,   pi)
      integer bc, N, j
      real dt, dx, fi(N), pifull(N), pi(N), xe(N), xc(N), Vdf,
     .     alpha, gamma, lambda, T, T0
*   --  The half-time quantity pi(N) evolved backwards to  -dt/2
*       for proper initial configuration.
      
      do 12 j = 2, N-1
         pi(j) = pifull(j) - 0.5*dt/(xc(j)*dx)**2
     .      * (xe(j)**2*(fi(j+1)-fi(j)) - xe(j-1)**2*(fi(j)-fi(j-1)))
     .      + 0.125*(dt/dx/xc(j))**2
     .      * (xe(j)**2*(pifull(j+1)-pifull(j))
     .         - xe(j-1)**2*(pifull(j)-pifull(j-1)))
     .      + 0.5*dt*Vdf(alpha,gamma,lambda,T,T0,
     .                   fi(j)-0.25*dt*pifull(j))
 12   continue
      call WrapEven (bc, N,   pi)

      return
      end
*     of BackStep



      subroutine Ebalance (dx, N, kappa, E, gb, xe, xc,
     .   fi, pifull, Q, alpha, gamma, lambda, V0, T, T0,
     .      Etot, restE, kinE, kinphi, grdphi, VE, QoverE)
      integer N, j
      real dx, kappa(N), E(N), gb(N), xe(N), xc(N), fi(N), pifull(N),
     .     Q(N), alpha, gamma, lambda, V0, T(N), T0,
     .     Etot, restE, kinE, kinphi, grdphi, VE, QoverE, vol, Vf
      Etot = 0.
      restE = 0.
      kinE = 0.
      kinphi = 0.
      grdphi = 0.
      VE = 0.
      QoverE = 0.
      do 12 j = 2, N-1
         vol = xe(j)**3-xe(j-1)**3
         restE = restE + E(j)/gb(j)*vol
         kinE = kinE + kappa(j)*E(j)/gb(j)*(gb(j)**2-1)*vol
         kinphi = kinphi + 0.5*pifull(j)**2*vol
         grdphi = grdphi + 0.125*((fi(j+1)-fi(j-1))/dx)**2*vol
         VE = VE + Vf(alpha, gamma, lambda, V0, T, T0, fi(j))*vol
         QoverE = QoverE + Q(j)*vol
 12   continue
c     Etot = restE + kinE + kinphi + grdphi + VE
      Etot = restE + kinE + kinphi + grdphi
      vol = xe(N-1)**3
      restE = restE/vol
      kinE = kinE/vol
      kinphi = kinphi/vol
      grdphi = grdphi/vol
      VE = VE/vol
      QoverE = QoverE/Etot
      Etot = Etot/vol

      return
      end
*     of Ebalance



      subroutine WaVefi (xe, xc, N, gb, v, pi, fi, Q, E, kappa,
     .   dt, dx, cycle, C, alpha, gamma, lambda, V0, T, T0,
     .      vwall, jV, xwall, p1, p2, p3) 
      integer N, cycle, j, jV
      real xe(N), xc(N), gb(N), v(N), pi(N), fi(N), Q(N),
     .   E(N), kappa(N),
     .   dt, dx, vwall, xwall, p1, p2, p3,
     .   C, alpha, gamma, lambda, V0, T(N), T0,
     .   Vmax, Qmax, Vj, Vjp, Vjm, time, Vf, parabtop
      real xwaold, xshold, Deltat, s, p, p0

C      WaVefi measures velocities from phi derivatives

      jV = 2
      Vmax = Vf(alpha,gamma,lambda,V0,T(jV),T0,fi(jV))
      p1 = (kappa(1)-1.)*E(1)/gb(1)
      p2 = (kappa(2)-1.)*E(2)/gb(2)
      p3 = max(p1,p2)
      p = min(p1,p2)
      p0 = p
      p1 = p
      p2 = p3
      do 12 j = 3, N-1
         Vj = Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
         if (Vj.gt.Vmax) then
            jV = j
            Vmax = Vj
         end if
         p = (kappa(j)-1.)*E(j)/gb(j)
         if (p.gt.p2) then
            p1 = p0
            p2 = p
            p3 = p
         end if
         if (p.lt.p0) p0 = p
         if (p.lt.p3) p3 = p
 12   continue
      Vjm = Vf(alpha,gamma,lambda,V0,T(jV-1),T0,fi(jV-1))
      Vjp = Vf(alpha,gamma,lambda,V0,T(jV+1),T0,fi(jV+1))
      xwall = parabtop(xc(jV-1),xc(jV),xc(jV+1),Vjm,Vmax,Vjp)
      j = jV
      s = -0.5*dt*C*gb(j)
      piplus = ( (1+s)*pi(j) + dt*(
     +      (xe(j)**2*(fi(j+1)-fi(j)) - xe(j-1)**2*(fi(j)-fi(j-1)))
     +      /(xc(j)*dx)**2
     .      - Vj - 0.5*C*gb(j)*
     .      ( v(j-1)*(fi(j)-fi(j-1)) + v(j)*(fi(j+1)-fi(j)) )/dx
     .      ) ) / (1-s)
      vwall = -(pi(j)+piplus)/(fi(j+1)-fi(j-1))*dx
  
      return
      end
*     of WaVefi



      subroutine WaVeDt (xc, N, fi, Q, dt, cycle,
     .   alpha, gamma, lambda, V0, T, T0,
     .      vwall, vshock, jV, jQ, xwall, xshock, sprint) 
      integer N, cycle, j, jV, jQ
      real xc(N), fi(N), Q(N), dt, vwall, vshock, xwall, xshock,
     .   alpha, gamma, lambda, V0, T(N), T0,
     .   Vmax, Qmax, Vj, Vjp, Vjm, time, Vf, parabtop
      integer sprint
      real xwaold, xshold, Deltat

C      WaVeDt measures velocities using sprint*dt as Delta t

      xwaold = xwall
      xshold = xshock

      jV = 2
      jQ = 2
      Vmax = Vf(alpha,gamma,lambda,V0,T(jV),T0,fi(jV))
      Qmax = Q(jQ)
      do 12 j = 3, N-1
         Vj = Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
         if (Vj.gt.Vmax) then
            jV = j
            Vmax = Vj
         end if
         if (Q(j).gt.Qmax) then
            jQ = j
            Qmax = Q(j)
         end if
 12   continue
      Vjm = Vf(alpha,gamma,lambda,V0,T(jV-1),T0,fi(jV-1))
      Vjp = Vf(alpha,gamma,lambda,V0,T(jV+1),T0,fi(jV+1))
      xwall = parabtop(xc(jV-1),xc(jV),xc(jV+1),Vjm,Vmax,Vjp)
      xshock = parabtop(xc(jQ-1),xc(jQ),xc(jQ+1),Q(jQ-1),Qmax,Q(jQ+1))
      time = cycle*dt
      Deltat = real(sprint)*dt
      if (time.eq.0.) then
         vwall = 0.
         vshock = 0.
      else if (cycle.eq.sprint) then
         vwall = (xwall-xwaold)/Deltat
         vshock = 0.
      else
         vwall = (xwall-xwaold)/Deltat
         vshock = (xshock-xshold)/Deltat
      end if
  
      return
      end
*     of WaVeDt



      subroutine WallVelo (xc, N, fi, Q, dt, cycle,
     .   alpha, gamma, lambda, V0, T, T0,
     .      vwall, vshock, jV, jQ, xwall, xshock) 
      integer N, cycle, j, jV, jQ
      real xc(N), fi(N), Q(N), dt, vwall, vshock, xwall, xshock,
     .   alpha, gamma, lambda, V0, T(N), T0,
     .   Vmax, Qmax, Vj, Vjp, Vjm, time, Vf, parabtop

C      WallVelo measures velocities using t=0 as starting point
C       -  WaVeDt is more accurate

      jV = 2
      jQ = 2
      Vmax = Vf(alpha,gamma,lambda,V0,T(jV),T0,fi(jV))
      Qmax = Q(jQ)
      do 12 j = 3, N-1
         Vj = Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))
         if (Vj.gt.Vmax) then
            jV = j
            Vmax = Vj
         end if
         if (Q(j).gt.Qmax) then
            jQ = j
            Qmax = Q(j)
         end if
 12   continue
      Vjm = Vf(alpha,gamma,lambda,V0,T(jV-1),T0,fi(jV-1))
      Vjp = Vf(alpha,gamma,lambda,V0,T(jV+1),T0,fi(jV+1))
      xwall = parabtop(xc(jV-1),xc(jV),xc(jV+1),Vjm,Vmax,Vjp)
      xshock = parabtop(xc(jQ-1),xc(jQ),xc(jQ+1),Q(jQ-1),Qmax,Q(jQ+1))
      time = cycle*dt
      if (time.eq.0.) then
         vwall = 0.
         vshock = 0.
      else
         vwall = xwall/time
         vshock = xshock/time
      end if
  
      return
      end
*     of WallVelo



      subroutine Output (init, cycle, sprint, lprint, dt, dx, N, nres1,
     .   nres2, Nprint, xc, xe, fi, D, E, Q, Z, v, gb, kappa,
     .   alpha, gamma, lambda, V0, T, T0,
     .   Etot, restE, kinE, kinphi, grdphi, VE, QoverE,
     .   vwall, vshock, jV, jQ, xwall, xshock, mwall, mshock,
     .   slclist, Nlist, sprf, sprl)
      integer init, cycle, sprint, lprint, N, j, nres1, nres2, mwall,
     .         mshock, jj, jV, jQ, Nprint, Nlist, slclist(Nlist),
     .        sprf, sprl
      real dt, dx, xc(N), xe(N), fi(N), D(N), E(N), Q(N), Z(N),
     .     gb(N), kappa(N), v(N), time0, Vf, alpha, gamma, lambda, V0,
     .     T(N), T0, Etot, restE, kinE, kinphi, grdphi, VE, QoverE,
     .     vwall, vshock, xwall, xshock, truncsm
*   --  Output as a function of time to file.6
*       and as a function of x to file.7
*   --  Wall detail to file.8
*       Shock detail to file.9
      logical Listed
      external Listed

c     write(6,*) 'Output'
      time0 = 0.
      if (cycle.eq.0) then
         write(6,101)
 101     format (/' cycle      time      Etot       restE   ',
     .      '   kinE      kinphi     grdphi       VE   ',
     .      '    QoverE   vwall   vshock   xwall   xshock '/)
         write(8,11) 2*mwall+1
         write(9,11) 2*mshock+1
 11      format (i6)
      end if
      if (mod(cycle,sprint).eq.0.
     .   .or. cycle.ge.sprf .and. cycle.le.sprl) then
         write(6,12) cycle, time0+cycle*dt, 
     .      Etot, restE, kinE, kinphi, grdphi,
     .      VE, QoverE, vwall, vshock, xwall, xshock
 12      format (i6,  f13.5, 6f11.6, e10.2, 2f9.5, 2f8.2)
      end if
      
      if (lprint.eq.-1 .and. Listed(cycle,slclist,Nlist)
     .   .or. lprint.gt.0 .and. mod(cycle,lprint).eq.0) then
         write(7,14) cycle, time0+cycle*dt
 14      format (i6, f14.6)
         do 18 j = 1, Nprint, nres1
            write(7,16) xc(j), truncsm(fi(j)),
     .         truncsm(Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))),
     .         E(j)/gb(j), E(j), truncsm(Q(j)), T(j),
     .         xe(j), truncsm(v(j)),
c    .         D(j)/gb(j), D(j), 
     .         (kappa(j)-1.)*E(j)/gb(j)
 16         format (11e17.9)
 18      continue
         do 19 j = Nprint+nres2, N, nres2
            write(7,16) xc(j), truncsm(fi(j)),
     .         truncsm(Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))),
     .         E(j)/gb(j), E(j), truncsm(Q(j)), T(j),
     .         xe(j), truncsm(v(j)),
c    .         D(j)/gb(j), D(j),
     .         (kappa(j)-1.)*E(j)/gb(j)
 19      continue
         write(8,20) cycle, time0+cycle*dt, xwall
 20      format (i6, f12.6, e16.8)
         do 22 jj = jV-mwall, jV+mwall
            if (jj.lt.1) then
               j = 1
            else if (jj.gt.N) then
               j = N
            else
               j = jj
            end if
            write(8,16) xc(j), truncsm(fi(j)),
     .         truncsm(Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))),
     .         E(j)/gb(j), E(j), truncsm(Q(j)), T(j), 
     .         xe(j), truncsm(v(j))
 22      continue
         write(9,20) cycle, time0+cycle*dt, xshock
         do 24 jj = jQ-mshock, jQ+mshock
            if (jj.lt.1) then
               j = 1
            else if (jj.gt.N) then
               j = N
            else
               j = jj
            end if
            write(9,16) xc(j), truncsm(fi(j)),
     .         truncsm(Vf(alpha,gamma,lambda,V0,T(j),T0,fi(j))),
     .         E(j)/gb(j), E(j), truncsm(Q(j)), T(j),
     .         xe(j), truncsm(v(j))
 24      continue
      end if

      return
      end
*     of Output



      logical function Listed (c, list, Nlist)
      integer c, Nlist, list(Nlist), i
      Listed = .false.
      do 20 i = 1, Nlist
         if (c.eq.list(i)) Listed = .true.
 20   continue
      RETURN
      END
*     off Listed
     



      subroutine Raster (cycle, rprint, r2first, r2last,
     .      r2left, r2right, N, rres, E, gb, T, fi, v,
     .      emin, emax, Tmin, Tmax, fimin, fimax, vmin, vmax)
      integer cycle, rprint, r2first, r2left, r2right, r2last, N, rres,
     .        ncol, cshft, i, j
      parameter (ncol = 250, cshft = 5)
      real E(N), gb(N), T(N), fi(N), v(N),
     .     emin, emax, edif, Tmin, Tmax, Tdif,
     .     fimin, fimax, fidif, vmin, vmax, vdif
c     parameter (emin = 36.04, emax = 36.2,
c    .   Tmin = 0.99985, Tmax = 1.00005,
c    .   fimin = 0.0, fimax = 0.4, vmin = -0.0003, vmax = 0.0004)

      edif = emax-emin
      Tdif = Tmax-Tmin
      fidif = fimax-fimin
      vdif = vmax-vmin
      if (mod(cycle,rprint).eq.0) then
         write(11,110) (
     .      mod( int(ncol*(E(j)/gb(j)-emin)/edif), ncol )+cshft,
     .      j = 1, N, rres)
         write(12,110) (
     .      mod( int(ncol*(T(j)-Tmin)/Tdif), ncol )+cshft,
     .      j = 1, N, rres)
         write(13,110) (
     .      mod( int(ncol*(fi(j)-fimin)/fidif), ncol )+cshft,
     .      j = 1, N, rres)
         write(14,110) (
     .      mod( int(ncol*(v(j)-vmin)/vdif), ncol )+cshft,
     .      j = 1, N, rres)
 110     format (20i4)
      end if

      if (cycle.ge.r2first .and. cycle.le.r2last) then
         write(21,110) (
     .      mod( int(ncol*(E(j)/gb(j)-emin)/edif), ncol )+cshft,
     .      j = r2left, r2right)
         write(22,110) (
     .      mod( int(ncol*(T(j)-Tmin)/Tdif), ncol )+cshft,
     .      j = r2left, r2right)
         write(23,110) (
     .      mod( int(ncol*(fi(j)-fimin)/fidif), ncol )+cshft,
     .      j = r2left, r2right)
         write(24,110) (
     .      mod( int(ncol*(v(j)-vmin)/vdif), ncol )+cshft,
     .      j = r2left, r2right)
      end if

      return
      end
*     of Raster



      SUBROUTINE Tracer (dt0, dx, bc, N, v, xe, cycle, Ntr,   xtr)
      INTEGER bc, N, Ntr, cycle
      REAL dt0, dx, v(N), xe(N), xtr(Ntr)
*   --  Advance tracer particles
      INTEGER i, jl
      REAL dt
      dt = dt0
      if (cycle.eq.0) dt = 0.5*dt0
*   --  Should add boundary conditions!!
      do 20 i = 1, Ntr
         jl = int((xtr(i)-xe(1))/dx)+1
         v0 = v(jl)
         v1 = v(jl+1)
         xtr(i) = (xtr(i)+v0*dt+(v1-v0)*(0.5*xtr(i)-xe(jl))*dt/dx) /
     .             (1 - (v1-v0)*0.5*dt/dx)
 20   continue
      RETURN
      END
*     of Tracer



      SUBROUTINE Trout (cycle, dt, trprint, trprf, trprl, Ntr, xtr)
      INTEGER cycle, Ntr, trprint, trprf, trprl,   i
      REAL dt, xtr(Ntr)
      if (mod(cycle,trprint).eq.0 
     .   .or. cycle.ge.trprf .and. cycle.le.trprl)
     .   write(19,21) cycle*dt, (xtr(i), i = 1, Ntr)
 21   format (102f10.4)
      RETURN
      END
*     of Trout
      


    
      subroutine Field (dt, dx, C, bc, N, gb, v, xe, xc,   pi, fi,
     .      Vdold, alpha, gamma, lambda, T, T0, fiold, pifull)
      integer bc, N, j
      real dt, dx, C, gb(N), v(N), xe(N), xc(N), pi(N), fi(N), Vdold(N), 
     .     fiold(N), pifull(N), piold, s,
     .     alpha, gamma, lambda, T(N), T0
*   --  Evolve the scalar field:  fi to the new time level,
*       pi half-way

c     write(6,*) 'Field'
      call Vdpot (bc, N, alpha, gamma, lambda, T, T0, fi,   Vdold)
      do 12 j = 2, N-1
         piold = pi(j)
         s = -0.5*dt*C*gb(j)
         pi(j) = ( (1+s)*pi(j) + dt*(
     +      (xe(j)**2*(fi(j+1)-fi(j)) - xe(j-1)**2*(fi(j)-fi(j-1)))
     +      /(xc(j)*dx)**2
     .      - Vdold(j) - 0.5*C*gb(j)*
     .      ( v(j-1)*(fi(j)-fi(j-1)) + v(j)*(fi(j+1)-fi(j)) )/dx
     .      ) ) / (1-s)
         pifull(j) = piold + 1.5*(pi(j)-piold)
 12   continue
      call WrapEven (bc, N,   pi)
      call WrapEven (bc, N,   pifull)

      do 14 j = 2, N-1
         fiold(j) = fi(j)
         fi(j) = fi(j) + dt*pi(j)
 14   continue
      call WrapEven (bc, N,   fi)
      call WrapEven (bc, N,   fiold)

      return
      end
*     of Field



      subroutine Trnsport (dt, dx, bc, N, v, xe, xc,   D,   Delta, F)
      integer bc, N, j
      real dt, dx, v(N), xe(N), xc(N), D(N), Delta(N), F(N), s, r
*   --  van Leer monotonic transport (2. order only)
*   --  see notes NR 23 (90/7/10..18) and QH 18.6 (93/3/10)

c     write(6,*) 'Trnsport'
      s = dt/dx
*   --  obtain gradients Delta (monotonicity applied here)
      do 12 j = 2, N-1
         r = (D(j)-D(j-1))*(D(j+1)-D(j))
         if (r.gt.0.) then
            Delta(j) = 2*r/(D(j+1)-D(j-1))
         else
            Delta(j) = 0.
         end if
 12   continue
      if (bc.eq.0) then
         Delta(1) = Delta(N-1)
         Delta(N) = Delta(2)
      else if (bc.eq.1) then
         Delta(1) = 0
         Delta(2) = 0
         Delta(N-1) = 0
         Delta(N) = 0
      else
         write(6,*) 'Need to have b.c. for Trnsport Delta'
      end if

*   --  obtain fluxes F
      do 14 j = 2, N-1
         if (v(j).gt.0.) then
            F(j) = v(j)*xe(j)**2 * (D(j) + 0.5*(1.-v(j)*s)*Delta(j))
         else
            F(j) = v(j)*xe(j)**2 * (D(j+1) - 0.5*(1.+v(j)*s)*Delta(j+1))
         end if
 14   continue
      call WrapOdd (bc, N,   F)

*   --  advect D
      do 16 j = 2, N-1
         D(j) = D(j) - s*(F(j)-F(j-1))/xc(j)**2
 16   continue
      call WrapEven (bc, N,   D)

      return
      end
*     of Trnsport



      subroutine Donorr (dt, dx, bc, N, v, xe, xc,   D,   Delta, F)
      integer bc, N, j
      real dt, dx, v(N), xe(N), xc(N), D(N), Delta(N), F(N), s, r
*   --  Donor Cell transport

c     write(6,*) ' Donorr '
      s = dt/dx
      do 12 j = 2, N-1
         if (v(j).ge.0.) then
            F(j) = v(j)*xe(j)**2*D(j)
         else
            F(j) = v(j)*xe(j)**2*D(j+1)
         end if
 12   continue
      call WrapOdd (bc, N,   F)
      do 16 j = 2, N-1
         D(j) = D(j) - s*(F(j)-F(j-1))/xc(j)**2
 16   continue
      call WrapEven (bc, N,   D)
      
      return
      end
*     of Donorr




      subroutine TransZ (dt, dx, bc, N, v, xe, xc,   Z,   Delta, F)
      integer bc, N, j
      real dt, dx, v(N), xe(N), xc(N), Z(N), Delta(N), F(N), s, r, vc
*   --  Transport for Z.  Like Trnsport, but shifted half dx,
*       and Z treated like momentum for boundary conditions
*       (changes sign for reflective bc)

c     write(6,*) 'TransZ'
      s = dt/dx
*   --  obtain gradients Delta (monotonicity applied here)
      do 12 j = 2, N-1
         r = (Z(j)-Z(j-1))*(Z(j+1)-Z(j))
         if (r.gt.0.) then
            Delta(j) = 2*r/(Z(j+1)-Z(j-1))
         else
            Delta(j) = 0.
         end if
 12   continue
      if (bc.eq.0) then
         Delta(1) = Delta(N-1)
         Delta(N) = Delta(2)
      else if (bc.eq.1) then
         Delta(1) = Z(2)
*      --  (sic)
         Delta(N) = Delta(N-2)
      else
         write(6,*) 'Need to have b.c. for TransZ Delta'
      end if

*   --  obtain fluxes F
      do 14 j = 2, N-1
         vc = 0.5*(v(j-1)+v(j))
         if (vc.gt.0.) then
            F(j) = vc*xc(j)**2 * (Z(j-1) + 0.5*(1.-vc*s)*Delta(j-1))
         else
            F(j) = vc*xc(j)**2 * (Z(j) - 0.5*(1.+vc*s)*Delta(j))
         end if
 14   continue
      call WrapEven (bc, N,   F)

*   --  advect Z
      do 16 j = 2, N-1
         Z(j) = Z(j) - s*(F(j+1)-F(j))/xe(j)**2
 16   continue
      call WrapOdd (bc, N,   Z)

      return
      end
*     of TransZ


      subroutine DonorZr (dt, dx, bc, N, v, xe, xc,   Z,   Delta, F)
      integer bc, N, j
      real dt, dx, v(N), xe(N), xc(N), Z(N), Delta(N), F(N), s, r, vc
*   --  Donor Cell transport

c     write(6,*) ' DonorZr '
      s = dt/dx
      do 12 j = 2, N-1
         vc = 0.5*(v(j-1)+v(j))
         if (vc.ge.0.) then
            F(j) = vc*xc(j)**2 * Z(j-1)
         else
            F(j) = vc*xc(j)**2 * Z(j)
         end if
 12   continue
      call WrapEven (bc, N,   F)
      do 16 j = 2, N-1
         Z(j) = Z(j) - s*(F(j+1)-F(j))/xe(j)**2
 16   continue
      call WrapOdd (bc, N,   Z)

      return
      end
*     of DonorZr




      subroutine EqState (a, bc, N, E, gb, fi, alpha, gamma,
     .    lambda, V0, T0, test,   T,   p, kappa, Vnew)
      integer bc, N, j, test
      real E(N), gb(N), fi(N), T(N), p(N), kappa(N),
     .     a, alpha, gamma, lambda, V0, T0, tolE, Vnew(N)
      parameter (tolE = 1.e-3)

*   --  Obtain pressure, temperature, and  kappa-1 = gb*p/E  from E.
*   --  This pressure is used only for pressure gradient acceleration of
*       Z, everywhere else use (kappa-1)E/gb.

c      write(6,*) 'EqState'
C      do 12 j = 2, N-1
C         call FindT(E(j), gb(j), fi(j), a, alpha, gamma,
C     .       lambda, V0, T0, j,   T(j)) 
C 12   continue
      call FindTa( E, gb, fi, a, alpha, gamma,
     .    lambda, V0, T0, N, test,   T ) 
      call Vpot (bc, N, alpha, gamma, lambda, V0, T, T0, fi,   Vnew)
      do 14 j = 2, N-1
         if (E(j).lt.tolE*gb(j)*3.*a*T(j)**4) then
            write(6,*) ' E getting dangerously small due to negative',
     .         '  V  contribution '
            write(6,125) j, E(j), T(j), Vnew(j)
 125        format (' j', i6, ' E', e12.4, ' T', f9.6, ' V', e12.4)
            test = 2
            return
         end if
         p(j) = a*T(j)**4-Vnew(j)
         kappa(j) = 1.+gb(j)*p(j)/E(j)
 14   continue
      call WrapEven (bc, N,   T)
      call WrapEven (bc, N,   p)
      call WrapEven (bc, N,   kappa)

      return
      end
*     of EqState




      subroutine FindTa (E, gb, fi, a, alpha, gamma, lambda, V0, T0, N,
     .   test,   T) 
      integer N, j, test
      real E(N), gb(N), fi(N), T(N),
     .     a, alpha, gamma, lambda, V0, T0, Tfix
C      real Ttst

*   --  Given the vectors E and fi, find the temperature vector T
C       by solving the quartic equation

C      Ttst = T(100)

      if (test.eq.1) write(6,*) 'FindTa'
      do 232 j=2,N-1
c        if (test.eq.1) write(6,100) j, fi(j), E(j), gb(j)
         Tfix = gamma**2/4.*fi(j)**4 - 12.*a*(lambda/4.*fi(j)**4
     .       - gamma/2.*T0**2*fi(j)**2 - E(j)/gb(j) + V0)
c      --  This if-write stuff ruins the vectorization: code takes 5*cpu time
c        if (Tfix.lt.0.) then
c           write(6,*) 'FindTa:  Tfix = ', Tfix,
c    .         ' set to zero at  j = ', j
c           Tfix = 0.0
c        end if
         T(j) = sqrt( 1./6./a * (gamma/2.*fi(j)**2 + sqrt(Tfix)) ) 
c        if (test.eq.1 .and. Tfix.eq.0.0)
c    .      write(6,100) j, fi(j), E(j), gb(j), T(j)
c100     format (i6, 4e12.4)

C         Tcons = lambda/4.*fi(j)**4 
C     .            - gamma/2.*T0**2*fi(j)**2 - E(j)/gb(j) + V0 
C         Tsti = gamma**2/4.*fi(j)**4 - 12.*a*Tcons
C         Tso = gamma/2.*fi(j)**2 + sqrt(Tsti)
C         Tapu = sqrt( 1./6./a * Tso )
C         T(j) = Tapu
C         write(6,*) Tcons, Tsti, Tso, Tapu

C         if (j.eq.100) then
C            call FindT(E(j), gb(j), fi(j), a, alpha, gamma,
C     .         lambda, V0, T0, j,   Ttst) 
C            write(6,*) j, ':  T_anal = ', T(j), '   T_num = ', Ttst
C         end if
232   continue

      return   
      end
*     of FindTa




      subroutine FindT (E, gb, fi, a, alpha, gamma, lambda, V0, T0, j,
     .      T) 
      integer j, k, kmax
      parameter (kmax = 20)
      real E, gb, fi, a, alpha, gamma, lambda, V0, T0, T, ET, ETD,
     . Tin, dT, tolflat, bound, acc, Vf, VTf, VTTf
      parameter (tolflat = 0.01, bound = 0.05, acc = 1.e-12)

*   --  Given E and fi, find the temperature T, by searching
*       the root of E(T)-E, using the Newton-Raphson method
C     Unnecessarily complicated, nowadays FindTa is used instead

         Tin = T
         do 12 k = 1, kmax
            ET = gb*(3*a*T**4 + Vf(alpha,gamma,lambda,V0,T,T0,fi)
     .           - T*VTf(alpha,gamma,lambda,T,T0,fi))
            ETd = gb*(12*a*T**3 - T*VTTf(alpha,gamma,lambda,T,T0,fi))
            if (abs(ETd/ET).lt.tolflat/T) then
               write(6,102) j, T, ET, ETd
 102           format (' FindT: E(T) is too flat at j =', i6,
     .            ' T ', f9.3, ' E ', f9.3, ' dE/dT', e10.2)
               stop
            end if
            dT = (ET-E)/ETd
            T = T-dT
            if (abs((T-Tin)/Tin).gt.bound) then
               write(6,104) bound, j, Tin, E, T, ET, ETd
 104           format (' FindT: T changing by more than', f6.3,
     .            ' at j =', i6, ' Tin', f9.3, ' E', f9.3, ' T', f9.3,
     .            ' E(T)', f9.3, ' dE/dT', f9.3)
               stop
            end if
            if (abs(dT/Tin).lt.acc) return
 12      continue
         write(6,122) j, Tin, E, T, ET, dT, ETd
 122     format (' FindT exceeding max iterations at j =', i6,
     .      ' Tin', e20.12, ' E    ', e20.12 /
     .      ' T  ', e20.12, ' E(T) ', e20.12 /
     .      ' dT ', e20.12, ' dE/dT', e20.12) 

         stop
         end
*        of FindT


      
      subroutine Hydro (dt, dx, kappa, C, Cav, mD, bc, N,
     .   fiold, fi, pi, p, D, xe, xc, alpha, gamma, lambda, T0, a,
     .      E, Z, v, gb, T,   Q, gbold, dxfi, fiav, gbv, Vdmid) 
      integer bc, N, j
      real dt, dx, kappa(N), C, Cav, mD, fiold(N), fi(N), pi(N), p(N),
     .     D(N), xe(N), xc(N), alpha, gamma, lambda, T0, a,
     .     E(N), Z(N), v(N), gb(N), T(N), Q(N), gbold(N),
     .     dxfi(N), gbv(N), dv, gv, s, fiav(N), Vdmid(N), gpi
*   --  Does hydro except transport.  E  and  Z evolved.  v  and gb
*       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.  
*   --  The order chosen is not unique, but is what I used in
*       dissertation.  Should experiment with different orders,
*       especially with position of field-fluid interaction,
*       now placed first (otherwise order is from CW).

c     write(6,*) 'Hydro'
      do 12 j = 2, N-1
         dxfi(j) = 0.5*(fiold(j+1)+fi(j+1)-fiold(j)-fi(j))/dx
         fiav(j) = 0.5*(fiold(j)+fi(j))
 12   continue
      call WrapOdd (bc, N,   dxfi)
      call WrapEven (bc, N,   fiav)
      call Vdpot (bc, N, alpha, gamma, lambda, T, T0, fiav,   Vdmid)
     
*   --  The field-fluid interaction + art.viscosity on  E.
      do 14 j = 2, N-1
         Z(j) = Z(j) - dt*0.5* ( C*(gb(j)+gb(j+1)) *
     .      ( 0.5*(pi(j)+pi(j+1)) + v(j)*dxfi(j) )
     .      +(Vdmid(j)+Vdmid(j+1)) ) * dxfi(j)
         gpi = gb(j) * ( pi(j) + 0.5*(v(j-1)*dxfi(j-1)+v(j)*dxfi(j)) )
         E(j) = E(j) + dt * (C*gpi**2 + Vdmid(j)*gpi)
         dv = v(j)-v(j-1)
         if (dv.lt.0.) then
            Q(j) = Cav*(mD*D(j)+E(j))*dv*dv
         else
            Q(j) = 0.
         end if
         E(j) = E(j) - dt*Q(j)*gb(j)*dv/dx
 14   continue
      call WrapEven (bc, N,   Q)
      call WrapEven (bc, N,   E)

*   --  Pressure acceleration + art. viscosity on  Z.  Obtain  v.
*   --  The pressure used is advected old pressure
      do 16 j = 2, N-1
         Z(j) = Z(j) - dt*(p(j+1)-p(j)+Q(j+1)-Q(j))/dx
         gv = 2*Z(j) / ( mD*(D(j)+D(j+1)) + kappa(j)*(E(j)+E(j+1)) )
         v(j) = gv/sqrt(1.+gv*gv)
 16   continue
      call WrapOdd (bc, N,   Z)
      call WrapOdd (bc, N,   v)
         
*   --  Obtain boost factor  gb,  pressure work on boost.
      do 18 j = 2, N-1
         gbold(j) = gb(j)
         gb(j) = sqrt( 1. +
     .      ( 0.5*(Z(j-1)+Z(j)) / (mD*D(j)+kappa(j)*E(j)) )**2 )
         s = (kappa(j)-1.)*(gb(j)-gbold(j))/(gbold(j)+gb(j))
         E(j) = E(j)*(1.-s)/(1.+s)
 18   continue
      call WrapEven (bc, N,   gb)
      call WrapEven (bc, N,   gbold)

*   --  Pressure work on divergence.
*   --  Like CW, use time-averaged  gb,  but new  v.
      do 20 j = 2, N-1
         gbv(j) = v(j)*xe(j)**2*(gbold(j)+gb(j)+gbold(j+1)+gb(j+1))/4
 20   continue
      call WrapOdd (bc, N,   gbv)
      do 22 j = 2, N-1
         s = (kappa(j)-1.)*dt/dx*(gbv(j)-gbv(j-1))
     .      /(gbold(j)+gb(j))/xc(j)**2
         E(j) = E(j)*(1.-s)/(1.+s)
 22   continue
      call WrapEven (bc, N,   E)

      return
      end
*     of Hydro



      subroutine WrapEven (bc, N,   D)
      integer bc, N
      real D(N)
*   --  Applies boundary condition to even quantities that lie at zone
*       centers.

      if (bc.eq.0) then
         D(1) = D(N-1)
         D(N) = D(2)
      else if (bc.eq.1) then
         D(1) = D(2)
         D(N) = D(N-1)
      else
         write(6,*) 'WrapEven: No such boundary condition.'
         stop
      end if

      return
      end
*     of WrapEven



      subroutine WrapOdd (bc, N,   v)
      integer bc, N
      real v(N)
*   --  Applies boundary condition to odd quantities that lie at zone
*       edges.

      if (bc.eq.0) then
         v(1) = v(N-1)
         v(N) = v(2)     
      else if (bc.eq.1) then
         v(1) = 0.   
         v(N) = -v(N-2)
      else
         write(6,*) 'WrapOdd: No such boundary condition.'
         stop
      end if

      return
      end
*     of WrapOdd  



      real function Vf(alpha, gamma, lambda, V0, T, T0, fi)
      real alpha, gamma, lambda, V0, T, T0, fi
*   --  The effective potential.  Remember to modify Vdf and Vdpot too!

      Vf = V0 + 0.5*gamma*(T**2-T0**2)*fi**2 - alpha*T*fi**3/3
     .     + 0.25*lambda*fi**4

      return
      end
*     of Vf



      real function Vdf(alpha, gamma, lambda, T, T0, fi)
      real alpha, gamma, lambda, T, T0, fi
*   --  The fi-derivative of the effective potential.

      Vdf = gamma*(T**2-T0**2)*fi - alpha*T*fi**2 + lambda*fi**3

      return
      end
*     of Vdf



      real function VTf(alpha, gamma, lambda, T, T0, fi)
      real alpha, gamma, lambda, T, T0, fi
*   --  The T-derivative of the effective potential.

      VTf = gamma*T*fi**2 - alpha*fi**3/3

      return
      end
*     of VTf



      real function VTTf(alpha, gamma, lambda, T, T0, fi)
      real alpha, gamma, lambda, T, T0, fi
*   --  The second T-derivative of the effective potential.

      VTTf = gamma*fi**2

      return
      end
*     of VTTf



      subroutine Vpot (bc, N, alpha, gamma, lambda, V0, T, T0, fi,
     .      Vnew)
      integer bc, N, j
      real alpha, gamma, lambda, V0, T(N), T0, fi(N), Vnew(N)
*   --  The values of the effective potential 
*       stored in Vnew(N).  (Vectorizes).

c     write(6,*) 'Vpot'
      do 12 j = 2, N-1
         Vnew(j) = V0 + 0.5*gamma*(T(j)**2-T0**2)*fi(j)**2 
     .             - alpha*T(j)*fi(j)**3/3
     .             + 0.25*lambda*fi(j)**4
 12   continue
      call WrapEven (bc, N,   Vnew)

      return
      end
*     of Vdpot



      subroutine Vdpot (bc, N, alpha, gamma, lambda, T, T0, fi,   Vd)
      integer bc, N, j
      real alpha, gamma, lambda, T(N), T0, fi(N), Vd(N)
*   --  The values of the fi-derivative of the effective potential 
*       stored in Vd(N).  (Vectorizes).

      do 12 j = 2, N-1
         Vd(j) = gamma*(T(j)**2-T0**2)*fi(j) - alpha*T(j)*fi(j)**2
     .           + lambda*fi(j)**3
 12   continue
      call WrapEven (bc, N,   Vd)

      return
      end
*     of Vdpot



      real function truncsm(x)
      real x, tiny
      parameter (tiny = 1.e-30)
*   --  Truncate small numbers (from Cray) to zero, so small machines
*       can read them
      
      if (abs(x).lt.tiny) then
         truncsm = 0.
      else
         truncsm = x
      end if
      return
      end
*     of truncsm



      real function parabtop (x1, x2, x3, y1, y2, y3)
      real x1, x2, x3, y1, y2, y3, C1, C2, C3
      
      C1 = y1/(x1-x2)/(x1-x3)
      C2 = y2/(x2-x1)/(x2-x3)
      C3 = y3/(x3-x1)/(x3-x2)
      if (y1.eq.y2 .and. y2.eq.y3) then
         parabtop = x2  
      else
         parabtop = 0.5*(C1*(x2+x3)+C2*(x1+x3)+C3*(x1+x2))/(C1+C2+C3)
      end if

      return
      end
*     of parabtop



      real function psibar( x, lbar )
      real x, lbar, reps, thx2, fx
      parameter( reps=1.E-8 )
C   --  The one-d extremum bubble in dimensionless units

      thx2 = (tanh(x))**2
      if (abs(thx2).lt.reps) then
        psibar = 2./lbar *(1. - sqrt(1.-lbar))
      else if (abs(thx2-1).lt.reps) then
        psibar = 0.
      else
        fx = (1. - lbar/thx2)/(1. - 1./thx2)
        psibar = 2./fx *(1. - sqrt(1.-fx))
      end if

      return
      end
C     of psibar
     


      real function arcoth (x)
      real x

      if (abs(x).gt.1) then
        arcoth = log((x+1)/(x-1)) / 2
      else
        write(6,*) 'Wrong type of argument to arcoth;  x = ', x
      end if

      return
      end
C     of arcoth


      real function cthj (x)
      real x, rmaxst, y
      parameter( rmaxst=1.E16 )

      y = tanh(x)
      if (y.eq.0) then
        cthj = rmaxst
      else
        cthj = 1./y
      end if
      write(6,*) 'cthj = ', cthj

      return
      end
C     of cthj

*     End of Program HydroQTr
