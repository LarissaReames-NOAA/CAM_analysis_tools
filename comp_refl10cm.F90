module cmpref_mod

     use share_mod

     IMPLICIT NONE


     interface cmpref
        module procedure calcrefl10cm
        module procedure calcrefl10cm1d
        module procedure rayleighsoakwetgraupel 
        module procedure get_m_mix_nested
        module procedure get_m_mix
        module procedure mcomplexmaxwellgarnett
        module procedure make_rainnumber
     end interface cmpref

     CONTAINS

subroutine calcrefl10cm (qv, qc, qr, qs, qg, &
                t, p, dBZ, kte, ite, jte)

        IMPLICIT NONE

        INTEGER, INTENT(IN)                         :: kte, ite, jte
        REAL, INTENT(IN), DIMENSION(kte,jte,ite)    :: qv, qc, qr, qs, qg, t, p
        REAL, INTENT(OUT), DIMENSION(kte,jte,ite) :: dBZ

        INTEGER :: i, j, k


        do j = 1,jte
        do i = 1,ite
                call calcrefl10cm1d(qv(:,j,i), qc(:,j,i), qr(:,j,i), &
                         qs(:,j,i), qg(:,j,i), t(:,j,i), &
                        p(:,j,i), dBZ(:,j,i), 0.0, 1, kte)
        enddo
        enddo

end subroutine calcrefl10cm

!!! ASSUMES TOP->BOTTOM LEVEL ORDER
subroutine calcrefl10cm1d (qv1d, qc1d, qr1d, qs1d, qg1d, &
                t1d, p1d, dBZ, rand1, kts, kte)
      IMPLICIT NONE

      LOGICAL, PARAMETER :: iiwarm = .false.

!..Densities of rain, snow, graupel, and cloud ice.
      REAL, PARAMETER :: rhow = 1000.0
      REAL, PARAMETER :: rhos = 100.0
      REAL, PARAMETER :: rhog = 500.0
      REAL, PARAMETER :: rhoi = 890.0

!..Generalized gamma distributions for rain, graupel and cloud ice.
!.. N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.
      REAL, PARAMETER :: mur = 0.0
      REAL, PARAMETER :: mug = 0.0
      REAL, PARAMETER :: mui = 0.0

!..Sum of two gamma distrib for snow (Field et al. 2005).
!.. N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
!..    + Kap1*(M2/M3)**mus * D**mus * exp(-M2*Lam1*D/M3)]
!.. M2 and M3 are the (bms)th and (bms+1)th moments respectively
!.. calculated as function of ice water content and temperature.
      REAL, PARAMETER :: mus = 0.6357
      REAL, PARAMETER :: Kap0 = 490.6
      REAL, PARAMETER :: Kap1 = 17.46
      REAL, PARAMETER :: Lam0 = 20.78
      REAL, PARAMETER :: Lam1 = 3.29

!..Y-intercept parameter for graupel is not constant and depends on
!.. mixing ratio.  Also, when mug is non-zero, these become equiv
!.. y-intercept for an exponential distrib and proper values are
!.. computed based on same mixing ratio and total number concentration.
      REAL, PARAMETER :: gonvmin = 1.E2
      REAL, PARAMETER :: gonvmax = 1.E6

!..Mass power law relations:  mass = am*D**bm
!.. Snow from Field et al. (2005), others assume spherical form.
      REAL amr, amg, ami 
      REAL, PARAMETER :: bmr = 3.0
      REAL, PARAMETER :: ams = 0.069
      REAL, PARAMETER :: bms = 2.0
      REAL, PARAMETER :: bmg = 3.0
      REAL, PARAMETER :: bmi = 3.0

!..Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
!.. Rain from Ferrier (1994), ice, snow, and graupel from
!.. Thompson et al (2008). Coefficient fv is zero for graupel/ice.
      REAL, PARAMETER :: avr = 4854.0
      REAL, PARAMETER :: bvr = 1.0
      REAL, PARAMETER :: fvr = 195.0
      REAL, PARAMETER :: avs = 40.0
      REAL, PARAMETER :: bvs = 0.55
      REAL, PARAMETER :: fvs = 100.0
      REAL, PARAMETER :: avg = 442.0
      REAL, PARAMETER :: bvg = 0.89
      REAL, PARAMETER :: avi = 1493.9
      REAL, PARAMETER :: bvi = 1.0
      REAL, PARAMETER :: avc = 0.316946E8
      REAL, PARAMETER :: bvc = 2.0

!..Minimum microphys values
!.. R1 value, 1.E-12, cannot be set lower because of numerical
!.. problems with Paul Field's moments and should not be set larger
!.. because of truncation problems in snow/ice growth.
      REAL, PARAMETER :: R1 = 1.E-12
      REAL, PARAMETER :: R2 = 1.E-6
      REAL, PARAMETER :: eps = 1.E-15

!..Rhonot used in fallspeed relations (rhonot/rho)**.5 adjustment.
      REAL, PARAMETER :: rhonot = 101325.0/(287.05*298.0)

!..Water vapor and air gas constants at constant pressure
      REAL, PARAMETER :: Rv = 461.5
      REAL, PARAMETER :: oRv = 1./Rv
      REAL, PARAMETER :: R = 287.04

!> For snow moments conversions (from Field et al. 2005)
      REAL, DIMENSION(10) :: sa, sb

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte 
      REAL, INTENT(IN):: rand1
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                          qv1d, qc1d, qr1d, qs1d, qg1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT) :: dBZ

!..Local variables
      LOGICAL :: allowwetgraupel
      LOGICAL :: allowwetsnow
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho, rhof, nr1d
      REAL, DIMENSION(kts:kte):: rc, rr, nr, rs, rg

      REAL, DIMENSION(kts:kte):: ilamr, ilamg, N0r, N0g
      REAL, DIMENSION(kts:kte):: mvdr
      REAL, DIMENSION(kts:kte):: smob, smo2, smoc, smoz
      REAL:: oM3, M0, Mrat, slam1, slam2, xDs

      REAL, DIMENSION(kts:kte):: zerain, zesnow, zegraupel

      REAL:: N0exp, N0min, lamexp, lamr, lamg
      REAL:: a, b, loga, tc0
      DOUBLE PRECISION:: fmelts, fmeltg

      INTEGER:: i, k, k0, kbot, n
      LOGICAL, DIMENSION(kts:kte):: Lqr, Lqs, Lqg

      DOUBLE PRECISION :: cback, x, eta, fd
      REAL:: xslw1, ygra1, zans1

      REAL:: oams, obms, ocms
      REAL, DIMENSION(5,15) :: cce, ccg
      REAL, DIMENSION(15) ::  ocg1, ocg2
      REAL, DIMENSION(7) :: cie, cig
      REAL :: oig1, oig2, obmi
      REAL, DIMENSION(13) :: cre, crg
      REAL :: ore1, org1, org2, org3, obmr
      REAL, DIMENSION(18) :: cse, csg
      REAL, DIMENSION(12) :: cge, cgg
      REAL :: oge1, ogg1, ogg2, ogg3, oamg, obmg, ocmg
      !..Ice initiates with this mass (kg), corresponding diameter calc.
      !..Min diameters and mass of cloud, rain, snow, and graupel (m, kg).
      REAL, PARAMETER :: xm0i = 1.E-12
      REAL, PARAMETER :: D0c = 1.E-6
      REAL, PARAMETER :: D0r = 50.E-6
      REAL, PARAMETER :: D0s = 300.E-6
      REAL, PARAMETER :: D0g = 350.E-6
      REAL :: D0i, xm0s, xm0g

      REAL, PARAMETER :: Ntc = 100.E6
      REAL, PARAMETER :: Ntcmax = 1999.E6
      REAL, PARAMETER :: Sc = 0.632
      REAL :: Sc3

      sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
              0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
      sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
              0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

      amr = PI*rhow/6.0
      amg = PI*rhog/6.0
      ami = PI*rhoi/6.0

      allowwetsnow = .true.
      allowwetgraupel = .true.

      do k = kts, kte
         dBZ(k) = -35.0
      enddo
      call radar_init()

!> - Compute Schmidt number to one-third used numerous times
      Sc3 = Sc**(1./3.)

!> - Compute minimum ice diam from mass, min snow/graupel mass from diam
      D0i = (xm0i/ami)**(1./bmi)
      xm0s = ams * D0s**bms
      xm0g = amg * D0g**bmg

!> - Compute constants various exponents and gamma() associated with cloud,
!! rain, snow, and graupel
      do n = 1, 15
         cce(1,n) = n + 1.
         cce(2,n) = bmr + n + 1.
         cce(3,n) = bmr + n + 4.
         cce(4,n) = n + bvc + 1.
         cce(5,n) = bmr + n + bvc + 1.
         call WGAMMA(ccg(1,n), cce(1,n))
         call WGAMMA(ccg(3,n), cce(2,n))
         call WGAMMA(ccg(3,n), cce(3,n))
         call WGAMMA(ccg(4,n), cce(4,n))
         call WGAMMA(ccg(5,n), cce(5,n))
         ocg1(n) = 1./ccg(1,n)
         ocg2(n) = 1./ccg(2,n)
      enddo

      cre(1) = bmr + 1.
      cre(2) = mur + 1.
      cre(3) = bmr + mur + 1.
      cre(4) = bmr*2. + mur + 1.
      cre(5) = mur + bvr + 1.
      cre(6) = bmr + mur + bvr + 1.
      cre(7) = bmr*0.5 + mur + bvr + 1.
      cre(8) = bmr + mur + bvr + 3.
      cre(9) = mur + bvr + 3.
      cre(10) = mur + 2.
      cre(11) = 0.5*(bvr + 5. + 2.*mur)
      cre(12) = bmr*0.5 + mur + 1.
      cre(13) = bmr*2. + mur + bvr + 1.
      do n = 1, 13
         call WGAMMA(crg(n), cre(n))
      enddo
      obmr = 1./bmr
      ore1 = 1./cre(1)
      org1 = 1./crg(1)
      org2 = 1./crg(2)
      org3 = 1./crg(3)

      cse(1) = bms + 1.
      cse(2) = bms + 2.
      cse(3) = bms*2.
      cse(4) = bms + bvs + 1.
      cse(5) = bms*2. + bvs + 1.
      cse(6) = bms*2. + 1.
      cse(7) = bms + mus + 1.
      cse(8) = bms + mus + 2.
      cse(9) = bms + mus + 3.
      cse(10) = bms + mus + bvs + 1.
      cse(11) = bms*2. + mus + bvs + 1.
      cse(12) = bms*2. + mus + 1.
      cse(13) = bvs + 2.
      cse(14) = bms + bvs
      cse(15) = mus + 1.
      cse(16) = 1.0 + (1.0 + bvs)/2.
      cse(17) = cse(16) + mus + 1.
      cse(18) = bvs + mus + 3.
      do n = 1, 18
         call WGAMMA(csg(n), cse(n))
      enddo
      oams = 1./ams
      obms = 1./bms
      ocms = oams**obms

      cge(1) = bmg + 1.
      cge(2) = mug + 1.
      cge(3) = bmg + mug + 1.
      cge(4) = bmg*2. + mug + 1.
      cge(5) = bmg*2. + mug + bvg + 1.
      cge(6) = bmg + mug + bvg + 1.
      cge(7) = bmg + mug + bvg + 2.
      cge(8) = bmg + mug + bvg + 3.
      cge(9) = mug + bvg + 3.
      cge(10) = mug + 2.
      cge(11) = 0.5*(bvg + 5. + 2.*mug)
      cge(12) = 0.5*(bvg + 5.) + mug
      do n = 1, 12
         call WGAMMA(cgg(n), cge(n))
      enddo
      oamg = 1./amg
      obmg = 1./bmg
      ocmg = oamg**obmg
      oge1 = 1./cge(1)
      ogg1 = 1./cgg(1)
      ogg2 = 1./cgg(2)
      ogg3 = 1./cgg(3)

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         rhof(k) = SQRT(RHONOT/rho(k))
         rc(k) = MAX(R1, qc1d(k)*rho(k))
         if (qr1d(k) .gt. R1 ) then
            rr(k) = qr1d(k)*rho(k)
            call make_rainnumber(nr1d(k),rr(k),temp(k))
            nr(k) = MAX(R2, nr1d(k))
            lamr = (amr*crg(3)*org2*nr(k)/rr(k))**obmr
            ilamr(k) = 1./lamr
            N0r(k) = nr(k)*org2*lamr**cre(2)
            mvdr(k) = (3.0 + mur + 0.672) * ilamr(k)
            Lqr(k) = .true.
         else
            rr(k) = R1
            nr(k) = R1
            mvdr(k) = 50.E-6
            Lqr(k) = .false.
         endif
         if (qs1d(k) .gt. R2) then
            rs(k) = qs1d(k)*rho(k)
            Lqs(k) = .true.
         else
            rs(k) = R1
            Lqs(k) = .false.
         endif
         if (qg1d(k) .gt. R2) then
            rg(k) = qg1d(k)*rho(k)
            Lqg(k) = .true.
         else
            rg(k) = R1
            Lqg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         smo2(k) = 0.
         smob(k) = 0.
         smoc(k) = 0.
         smoz(k) = 0.
      enddo
      if (ANY(Lqs .eqv. .true.)) then
      do k = kts, kte
         if (.not. Lqs(k)) CYCLE
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams

!..All other moments based on reference, 2nd moment.  If bms.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bms.gt.(2.0-1.e-3) .and. bms.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga = sa(1) + sa(2)*tc0 + sa(3)*bms &
     &         + sa(4)*tc0*bms + sa(5)*tc0*tc0 &
     &         + sa(6)*bms*bms + sa(7)*tc0*tc0*bms &
     &         + sa(8)*tc0*bms*bms + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*bms*bms*bms
            a = 10.0**loga
            b = sb(1) + sb(2)*tc0 + sb(3)*bms &
     &         + sb(4)*tc0*bms + sb(5)*tc0*tc0 &
     &         + sb(6)*bms*bms + sb(7)*tc0*tc0*bms &
     &         + sb(8)*tc0*bms*bms + sb(9)*tc0*tc0*tc0 &
     &         + sb(10)*bms*bms*bms
            smo2(k) = (smob(k)/a)**(1./b)
         endif

!..Calculate bms+1 (th) moment.  Useful for diameter calcs.
         loga = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
     &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
     &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(1)*cse(1)*cse(1)
         a = 10.0**loga
         b = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
     &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a * smo2(k)**b

!..Calculate bms*2 (th) moment.  Useful for reflectivity.
         loga = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
     &         + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
     &         + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(3)*cse(3)*cse(3)
         a = 10.0**loga
         b = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
     &        + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
         smoz(k) = a * smo2(k)**b
      enddo
      endif

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+

      if (ANY(Lqg .eqv. .true.)) then
      do k = kte, kts, -1
         ygra1 = alog10(max(1.E-9, rg(k)))
         zans1 = 3.0 + 2./7.*(ygra1+8.) + rand1
         N0exp = 10.**(zans1)
         N0exp = MAX(DBLE(gonvmin), MIN(N0exp, DBLE(gonvmax)))
         lamexp = (N0exp*amg*cgg(1)/rg(k))**oge1
         lamg = lamexp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0g(k) = N0exp/(cgg(2)*lamexp) * lamg**cge(2)
      enddo
      endif

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k0 is level above).
!+---+-----------------------------------------------------------------+
      k0 = kts
      K_LOOP:do k = kte-1, kts, -1 
        if ((temp(k).gt.273.15) .and. Lqr(k)                         &
     &                            .and. (Lqs(k+1).or.Lqg(k+1)) ) then
           k0 = MAX(k+1, k0)
           EXIT K_LOOP
        endif
      enddo K_LOOP
!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         zerain(k) = 1.e-22
         zesnow(k) = 1.e-22
         zegraupel(k) = 1.e-22
         if (Lqr(k)) zerain(k) = N0r(k)*crg(4)*ilamr(k)**cre(4)
         if (Lqs(k)) zesnow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
     &                           * (ams/900.0)*(ams/900.0)*smoz(k)
         if (Lqg(k)) zegraupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
     &                              * (amg/900.0)*(amg/900.0)         &
     &                              * N0g(k)*cgg(4)*ilamg(k)**cge(4)
      enddo

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleighsoakwetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (.not. iiwarm .and. k0.ge.2) then
       do k = k0-1, kts, -1

!..Reflectivity contributed by melting snow
          if (allowwetsnow .and. Lqs(k) .and. Lqs(k0) ) then
           fmelts = MAX(0.05d0, MIN(1.0d0-rs(k)/rs(k0), 0.99d0))
           eta = 0.d0
           oM3 = 1./smoc(k)
           M0 = (smob(k)*oM3)
           Mrat = smob(k)*M0*M0*M0
           slam1 = M0 * Lam0
           slam2 = M0 * Lam1
           do n = 1, nrbins
              x = ams * xxDs(n)**bms
              call rayleighsoakwetgraupel (x, DBLE(ocms), DBLE(obms), &
     &              fmelts, meltoutsides, mw0, mi0, lamdaradar, &
     &              CBACK, mixingrulestrings, matrixstrings,          &
     &              inclusionstrings, hoststrings,                    &
     &              hostmatrixstrings, hostinclusionstrings)
              fd = Mrat*(Kap0*DEXP(-slam1*xxDs(n))                     &
     &              + Kap1*(M0*xxDs(n))**mus * DEXP(-slam2*xxDs(n)))
              eta = eta + fd * CBACK * simpson(n) * xdts(n)
           enddo
           zesnow(k) = SNGL(lamda4 / (pi5 * Kw) * eta)
          endif

!..Reflectivity contributed by melting graupel
          if (allowwetgraupel .and. Lqg(k) .and. Lqg(k0) ) then
           fmeltg = MAX(0.05d0, MIN(1.0d0-rg(k)/rg(k0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = amg * xxDg(n)**bmg
              call rayleighsoakwetgraupel (x, DBLE(ocmg), DBLE(obmg), &
     &              fmeltg, meltoutsideg, mw0, mi0, lamdaradar, &
     &              CBACK, mixingrulestringg, matrixstringg,          &
     &              inclusionstringg, hoststringg,                    &
     &              hostmatrixstringg, hostinclusionstringg)
              fd = N0g(k)*xxDg(n)**mug * DEXP(-lamg*xxDg(n))
              eta = eta + fd * CBACK * simpson(n) * xdtg(n)
           enddo
           zegraupel(k) = SNGL(lamda4 / (pi5 * Kw) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((zerain(k)+zesnow(k)+zegraupel(k))*1.d18)
      enddo

end subroutine calcrefl10cm1d



     subroutine rayleighsoakwetgraupel (xg, ageo, bgeo, fmelt,    &
                     meltratiooutside, mw, mi, lambda, Cback,       &
                     mixingrule,matrix,inclusion,                       &
                     host,hostmatrix,hostinclusion)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in):: xg, ageo, bgeo, fmelt, lambda,  &
                                     meltratiooutside
      DOUBLE PRECISION, INTENT(out):: Cback
      COMPLEX*16, INTENT(in):: mw, mi
      CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion,     &
                                     host, hostmatrix, hostinclusion

      COMPLEX*16:: mcore, mair
      DOUBLE PRECISION:: Dlarge, Dg, rhog, xw, xwa, fm, fmgrenz,    &
                         volg, vg, volair, volice, volwater,            &
                         meltratiooutsidegrenz, mra
      INTEGER:: error
      DOUBLE PRECISION:: PIx
      
      PIx = 3.1415926535897932384626434d0

!     refractive index of air:
      mair = (1.0d0,0.0d0)

!     Limiting the degree of melting --- for safety:
      fm = DMAX1(DMIN1(fmelt, 1.0d0), 0.0d0)
!     Limiting the ratio of (melting on outside)/(melting on inside):
      mra = DMAX1(DMIN1(meltratiooutside, 1.0d0), 0.0d0)

!    ! The relative portion of meltwater melting at outside should increase
!    ! from the given input value (between 0 and 1)
!    ! to 1 as the degree of melting approaches 1,
!    ! so that the melting particle "converges" to a water drop.
!    ! Simplest assumption is linear:
      mra = mra + (1.0d0-mra)*fm

      xw = xg * fm

      Dg = ageo * xg**bgeo

      if (Dg .ge. 1d-12) then

       vg = PIx/6. * Dg**3
       rhog = DMAX1(DMIN1(xg / vg, 900.0d0), 10.0d0)
             vg = xg / rhog

       meltratiooutsidegrenz = 1.0d0 - rhog / 1000.

       if (mra .le. meltratiooutsidegrenz) then
        !..In this case, it cannot happen that, during melting, all the
        !.. air inclusions within the ice particle get filled with
        !.. meltwater. This only happens at the end of all melting.
        volg = vg * (1.0d0 - mra * fm)

       else
        !..In this case, at some melting degree fm, all the air
        !.. inclusions get filled with meltwater.
        fmgrenz=(900.0-rhog)/(mra*900.0-rhog+900.0*rhog/1000.)

        if (fm .le. fmgrenz) then
         !.. not all air pockets are filled:
         volg = (1.0 - mra * fm) * vg
        else
         !..all air pockets are filled with meltwater, now the
         !.. entire ice sceleton melts homogeneously:
         volg = (xg - xw) / 900.0 + xw / 1000.
        endif

       endif

       Dlarge  = (6.0 / PIx * volg) ** (1./3.)
       volice = (xg - xw) / (volg * 900.0)
       volwater = xw / (1000. * volg)
       volair = 1.0 - volice - volwater

       !..complex index of refraction for the ice-air-water mixture
       !.. of the particle:
       call get_m_mix_nested (mcore, mair, mi, mw, volair, volice,      &
                         volwater, mixingrule,host, matrix, inclusion, &
                         hostmatrix,hostinclusion,error)
       if (error .ne. 0) then
        Cback = 0.0d0
        return
       endif

       !..Rayleigh-backscattering coefficient of melting particle:
       Cback = (ABS((mcore**2-1.0d0)/(mcore**2+2.0d0)))**2           &
                * PI5 * Dlarge**6 / lamda4

      else
       Cback = 0.0d0
      endif

      end subroutine rayleighsoakwetgraupel

      subroutine get_m_mix_nested (mmixn,ma, mi, mw, volair,      &
                     volice, volwater, mixingrule, host, matrix,        &
                     inclusion, hostmatrix, hostinclusion, cumulerror)

      IMPLICIT NONE

      COMPLEX*16, INTENT(OUT) :: mmixn
      DOUBLE PRECISION, INTENT(in):: volice, volair, volwater
      COMPLEX*16, INTENT(in):: ma, mi, mw
      CHARACTER(len=*), INTENT(in):: mixingrule, host, matrix,          &
                     inclusion, hostmatrix, hostinclusion
      INTEGER, INTENT(out):: cumulerror

      DOUBLE PRECISION:: vol1, vol2
      COMPLEX*16:: mtmp
      INTEGER:: error

      !..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be
      !.. air, ice, or water

      cumulerror = 0
      mmixn = CMPLX(1.0d0,0.0d0)

      if (host .eq. 'air') then

       if (matrix .eq. 'air') then
        write(*,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        cumulerror = cumulerror + 1
       else
        vol1 = volice / MAX(volice+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        call get_m_mix (mtmp, ma, mi, mw, 0.0d0, vol1, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'air') then
         call get_m_mix (mmixn,ma, mtmp, 2.0*ma,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'icewater') then
         call get_m_mix (mmixn,ma, mtmp, 2.0*ma,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error
        else
         write(*,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
         cumulerror = cumulerror + 1
          endif
       endif

      elseif (host .eq. 'ice') then

       if (matrix .eq. 'ice') then
        write(*,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volair+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        call get_m_mix (mtmp,ma, mi, mw, vol1, 0.0d0, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'ice') then
         call get_m_mix (mmixn,mtmp, mi, 2.0*ma,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airwater') then
         call get_m_mix (mmixn, mtmp, mi, 2.0*ma,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         'air', hostinclusion, error)
         cumulerror = cumulerror + error
        else
         write(*,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'water') then

       if (matrix .eq. 'water') then
        write(*,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volice+volair,1d-10)
        vol2 = 1.0d0 - vol1
        call get_m_mix (mtmp,ma, mi, mw, vol1, vol2, 0.0d0,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'water') then
         call get_m_mix (mmixn,2*ma, mtmp, mw,                &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airice') then
            call get_m_mix (mmixn,2*ma, mtmp, mw,                &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error
        else
         write(*,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',         &
                           hostmatrix
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'none') then

       call get_m_mix (mmixn,ma, mi, mw,                     &
                       volair, volice, volwater, mixingrule,            &
                       matrix, inclusion, error)
       cumulerror = cumulerror + error

      else
       write(*,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
       cumulerror = cumulerror + 1
      endif

      IF (cumulerror .ne. 0) THEN
       write(*,*) 'GET_M_MIX_NESTED: error encountered'
       mmixn = CMPLX(1.0d0,0.0d0)
      endif

      end subroutine get_m_mix_nested

      SUBROUTINE get_m_mix (mmix, ma, mi, mw, volair, volice,     &
                     volwater, mixingrule, matrix, inclusion, error)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(in):: volice, volair, volwater
      COMPLEX*16, INTENT(in):: ma, mi, mw
      CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion
      INTEGER, INTENT(out):: error
      COMPLEX*16, INTENT(OUT) :: mmix

      error = 0
      mmix = CMPLX(1.0d0,0.0d0)

      if (mixingrule .eq. 'maxwellgarnett') then
       if (matrix .eq. 'ice') then
        call mcomplexmaxwellgarnett(mmix,volice, volair, volwater,  &
                           mi, ma, mw, inclusion, error)
       elseif (matrix .eq. 'water') then
        call mcomplexmaxwellgarnett(mmix,volwater, volair, volice,  &
                           mw, ma, mi, inclusion, error)
       elseif (matrix .eq. 'air') then
        call mcomplexmaxwellgarnett(mmix,volair, volwater, volice,  &
                           ma, mw, mi, inclusion, error)
       else
        write(*,*) 'GET_M_MIX: unknown matrix: ', matrix
        error = 1
       endif

      else
       write(*,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
       error = 2
      endif

      if (error .ne. 0) then
       write(*,*) 'GET_M_MIX: error encountered'
      endif

      END SUBROUTINE get_m_mix

      SUBROUTINE mcomplexmaxwellgarnett(mmg, vol1, vol2, vol3,    &
                     m1, m2, m3, inclusion, error)
      IMPLICIT NONE

      COMPLEX*16, INTENT(OUT) :: mmg
      COMPLEX*16, INTENT(IN) :: m1, m2, m3
      DOUBLE PRECISION, INTENT(IN) :: vol1, vol2, vol3
      CHARACTER(len=*), INTENT(IN) :: inclusion

      COMPLEX*16 :: beta2, beta3, m1t, m2t, m3t
      INTEGER, INTENT(out) :: error

      error = 0

      if (DABS(vol1+vol2+vol3-1.0d0) .gt. 1d-6) then
       write(*,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ',       &
              'partial volume fractions is not 1...ERROR'
       mmg=CMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m1t = m1**2
      m2t = m2**2
      m3t = m3**2

      if (inclusion .eq. 'spherical') then
       beta2 = 3.0d0*m1t/(m2t+2.0d0*m1t)
       beta3 = 3.0d0*m1t/(m3t+2.0d0*m1t)
      elseif (inclusion .eq. 'spheroidal') then
       beta2 = 2.0d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0d0)
       beta3 = 2.0d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0d0)
      else
       write(*,*) 'M_COMPLEX_MAXWELLGARNETT: ',                  &
                         'unknown inclusion: ', inclusion
       mmg=DCMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      mmg = &
       SQRT(((1.0d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
       (1.0d0-vol2-vol3+vol2*beta2+vol3*beta3))

      END SUBROUTINE mcomplexmaxwellgarnett

      subroutine make_rainnumber (rnc,Qrain, temp)

      IMPLICIT NONE

      real, intent(out):: rnc
      real, intent(in):: Qrain, temp
      double precision:: lambda, N0, qnr
      !real, parameter:: PI = 3.1415926536
      real, parameter:: amr = PI*1000./6.

      if (Qrain == 0) then
         rnc = 0
         return
      end if

      !+---+-----------------------------------------------------------------+
      !.. Not thrilled with it, but set Y-intercept parameter to Marshal-Palmer value
      !.. that basically assumes melting snow becomes typical rain. However, for
      !.. -2C < T < 0C, make linear increase in exponent to attempt to keep
      !.. supercooled collision-coalescence (warm-rain) similar to drizzle rather
      !.. than bigger rain drops.  While this could also exist at T>0C, it is
      !.. more difficult to assume it directly from having mass and not number.
      !+---+-----------------------------------------------------------------+

      N0 = 8.E6

      if (temp .le. 271.15) then
         N0 = 8.E8
      elseif (temp .gt. 271.15 .and. temp.lt.273.15) then
         N0 = 8. * 10**(279.15-temp)
      endif

      lambda = SQRT(SQRT(N0*amr*6.0/Qrain))
      qnr = Qrain / 6.0 * lambda*lambda*lambda / amr
      rnc = SNGL(qnr)

      end subroutine make_rainnumber

 end module cmpref_mod
