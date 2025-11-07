! Wall model based on Monin-Obukhov similarity theory, choosing
! correction functions for stable or convective conditions
! based on the computed Richardson number
! Linnea Huusko, 2025-05-22

!> @brief Rough log law in implicit form
      subroutine wmles_set_momentum_flux(i)
      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'FRAMELP'
      include 'WMLES'

      ! sampling height
      real h

      ! obukhov length
      real l_obukhov, l_upper, l_lower, l_backup, l_old

      ! the indices of the gll point
      integer i, ix, iy, iz, ie, count, max_count

      ! sampled velocity
      real magvh

      ! richardson number
      real rib, ri_limit

      real utau, utau_old, g, th, ts, q

      ! newton iteration stuff
      real f, dfdl, fd_h

      ! log law parameters
      real kappa, z0, z1

      ! similarity law for velocity and heat
      real tau_conv, heat_flux_conv
      real tau, heat_flux !similarity_law_u_stable, similarity_law_q_stable
      real l, N
      real fcor, C_f, C_N

      ! dummy variables for retrieving parameters
      integer itmp
      logical ltmp
      character*20 ctmp

      parameter(fcor = 1.39e-4)
      parameter(C_f = 0.185)
      parameter(C_N = 2.0)
!-----------------------------------------------------------------------

      ! assign kappa and B and z0
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,z0,ltmp,ctmp,wmles_z0_id,rpar_real)
      call rprm_rp_get(itmp,z1,ltmp,ctmp,wmles_z1_id,rpar_real)

      ri_limit = 0.001
      g = 9.81

      ix = wmles_indices(i, 1)
      iy = wmles_indices(i, 2)
      iz = wmles_indices(i, 3)
      ie = wmles_indices(i, 4)

      h = wmles_sampling_h(i)

      ! Sample the values at the sampling point
      ! Velocity (only horizontal components)
      magvh = wmles_solh(i, 1)**2 +
     $        wmles_solh(i, 3)**2
      magvh = sqrt(magvh)

      ! Temperature from sampling point
      th = wmles_solh(i, 4)

      ! ts = surface temperature
      ts = wmles_surface_temp

      ! Get uncorrected utau for a first guess other wise from last
      ! timestep (only horizontal components)
      if (ISTEP .eq. 0) then
        utau = wmles_tau(i, 1)**2 +
     $         wmles_tau(i, 3)**2
        utau = sqrt(sqrt(utau))
      else
        utau = magvh*kappa/log(h/z0)
      end if

      if (wmles_forcing_type .eq. surface_temperature) then
        q = kappa*utau*(ts - th)/log(h/z1)
      elseif (wmles_forcing_type .eq. surface_heat_flux) then
        q = wmles_q(i)
      elseif (wmles_forcing_type .eq. no_forcing_type) then
          write(*,*) "No surface forcing type set!"
          stop
      endif

      if (ISTEP .gt. 0) then
        if (wmles_forcing_type .eq. surface_temperature) then
          q = kappa*utau*(ts - th)/log(h/z1)
          rib = g*h/th*(th - ts)/magvh**2
        elseif (wmles_forcing_type .eq. surface_heat_flux) then
          q = wmles_q(i)
          rib = -g*h/th*q/(magvh**3*kappa**2)
        endif

        ! Obukhov l based on the previous-step utau
        l_obukhov = -(wmles_theta0*utau**3)/(kappa*g*q)
        wmles_lobukhov(i) = l_obukhov

        ! In case the iteration diverges we will just use this
        l_backup = l_obukhov

        l_old = 0
        count = 0
        max_count = 20

! ===== Neutral =======================================================
        if (abs(rib).lt.ri_limit) then ! neutral (use log law computed above)
          ! if (i.eq.1) then
          !   write(*,*) "Neutral", rib
          ! endif
          l_obukhov = 0

! ===== Convective ====================================================
        elseif (rib.lt.-ri_limit) then ! convective
          ! if (i.eq.1) then
          !   write(*,*) "Convective", rib
          ! endif

          utau = (tau_conv(wmles_solh(i, 1), wmles_solh(i, 3),
     $            rib, h, z0))**0.5

          q = heat_flux_conv(th, ts, wmles_solh(i, 1),
     $        wmles_solh(i, 3), rib, h, z1)

! ===== Stable ========================================================
        else ! stable
          ! if (i.eq.1) then
          !   write(*,*) "Stable", rib
          ! endif

          N = sqrt(g/wmles_theta0 * (th-ts)/h)
    !       l = 1/(1/(0.4*h)
    !  $     + fcor/(C_f*utau)
    !  $     + N/(C_N*utau))
          l = 0.4 * h
          utau = (tau(wmles_solh(i, 1), wmles_solh(i, 3),
     $            rib, h, z0, l))**0.5
          q = heat_flux(th, ts, rib, h, z1, 1.0, l, utau)

        endif

        ! store the computed Obukhov length and Richardson number
        wmles_lobukhov(i) = l_obukhov
        wmles_ri(i) = rib
        wmles_count(i) = l
        wmles_local_index(i) = N

      endif

      wmles_ustar(i) = utau

      ! Assign tau proportional to the velocity magnitudes at
      ! the sampling point
      wmles_tau(i, 1) = -utau**2*wmles_solh(i, 1)/magvh
      wmles_tau(i, 2) = 0
      wmles_tau(i, 3) = -utau**2*wmles_solh(i, 3)/magvh
      wmles_q(i) = q
      end

!--- Convective --------------------------------------------------------
!    Correction functions based on Louis 1979
      real function tau_conv(u, v, ri, h, z0)
      implicit none

      real a, b, c, F_m
      real u, v, ri, h, z0

      a = 0.4 / log(h/z0)
      b = 2
      c = 7.4 * a**2 * b * (h/z0)**0.5

      F_m = 1 - 2*ri / (1 + c * abs(ri)**0.5)

      tau_conv = a**2 * (u**2 + v**2) * F_m

      end function

!-----------------------------------------------------------------------
      real function heat_flux_conv(theta2, theta1, u, v, ri, h, z1)
      implicit none

      real a, b, c, F_h
      real theta2, theta1, u, v, ri, h, z1

      a = 0.4 / log(h/z1)
      b = 2
      c = 5.3 * a**2 * b * (h/z1)**0.5

      F_h = 1 - 2*ri / (1 + c * abs(ri)**0.5)

      heat_flux_conv = - a**2 / 0.74 * (u**2 + v**2)**0.5 *
     $ (theta2 - theta1) * F_h

      end function

!-----------------------------------------------------------------------
!--- Stable ------------------------------------------------------------
!    Correction functions based on Mauritsen et al. 2007
      real function f_tau(ri)
      implicit none

      real ri

      f_tau = 0.17 * (0.25 + 0.75 / (1 + 4*ri))

      end function
!-----------------------------------------------------------------------
      real function f_theta(ri)
      implicit none

      real ri

      f_theta = -0.145 / (1 + 4 * ri)

      end function
!-----------------------------------------------------------------------
      real function tau(u, v, ri, h, z0, l)
      implicit none

      real u, v, ri, h, z0, l
      real f_tau

      tau = (u**2 + v**2)/(log(h/z0)**2)
     $ * f_tau(ri)/f_tau(0) * (l/h)**2

      end function
!-----------------------------------------------------------------------
      real function heat_flux(theta2, theta1, ri, h, z1, pr, l, utau)
      implicit none

      real theta1, theta2, ri, h, z1, pr, l, utau
      real f_theta

      heat_flux = (theta2 - theta1)/(log(h/z1))
     $ * f_theta(ri)/abs(f_theta(0)) * (l/h)
     $ * utau/pr

      end function
