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
      integer i, ix, iy, iz, ie, count, flag

      ! sampled velocity
      real magvh

      ! richardson number
      real rib

      real utau, utau_old, g, th, ts, q

      ! newton iteration stuff
      real f, dfdl, fd_h

      ! log law parameters
      real kappa, z0, z1

      ! parameters in front of the correction functions for heat and
      ! momentum.
      real a, b


      ! similarity law for velocity and heat
      real similarity_law_u_conv, similarity_law_q_conv
      real similarity_law_u_stable, similarity_law_q_stable

      ! dummy variables for retrieving parameters
      integer itmp
      logical ltmp
      character*20 ctmp
!-----------------------------------------------------------------------
      flag = 1

      ! assign kappa and B and z0
      call rprm_rp_get(itmp,kappa,ltmp,ctmp,wmles_logkappa_id,rpar_real)
      call rprm_rp_get(itmp,z0,ltmp,ctmp,wmles_z0_id,rpar_real)
      call rprm_rp_get(itmp,z1,ltmp,ctmp,wmles_z1_id,rpar_real)

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
      if (ISTEP .lt. 3) then
        utau = wmles_tau(i, 1)**2 +
     $         wmles_tau(i, 3)**2
        utau = sqrt(sqrt(utau))
      else
        utau = magvh*kappa/log(h/z0)
      end if

      ! if (ISTEP .gt. 373) then
      !   if (abs(utau) .gt. 1 .or. flag .eq. 1) then
      !     write(*,*) ISTEP, i, "DEBUG A ", "utau =", utau
      !     flag = 1
      !   endif
      ! endif

      if (wmles_forcing_type .eq. surface_temperature) then
        q = kappa*utau*(ts - th)/log(h/z1)
      elseif (wmles_forcing_type .eq. surface_heat_flux) then
        q = wmles_q(i)
      elseif (wmles_forcing_type .eq. no_forcing_type) then
          write(*,*) "No surface forcing type set!"
          stop
      endif

      ! if (ISTEP .gt. 373) then
      !   if (abs(q) .gt. 1 .or. flag .eq. 1) then
      !     write(*,*) ISTEP, i, "DEBUG B ", "q =", q
      !     flag = 1
      !   endif
      ! endif


      if (ISTEP .gt. 3) then
        if (wmles_forcing_type .eq. surface_temperature) then
          q = kappa*utau*(ts - th)/log(h/z1)
          rib = g*h/th*(th - ts)/magvh**2
        elseif (wmles_forcing_type .eq. surface_heat_flux) then
          q = wmles_q(i)
          rib = -g*h/th*q/(magvh**3*kappa**2)
        endif

        ! if (ISTEP .gt. 373) then
        !   if (abs(q) .gt. 1 .or. abs(rib) .gt. 1 .or. flag .eq. 1) then
        !     write(*,*) ISTEP, i, "DEBUG C ", "q =", q, "rib =", rib
        !     flag = 1
        !   endif
        ! endif

        ! Obukhov l based on the previous-step utau
        l_obukhov = -(wmles_theta0*utau**3)/(kappa*g*q)
        wmles_lobukhov(i) = l_obukhov

        ! if (ISTEP .gt. 373) then
        !   if (abs(l_obukhov) .gt. 1000 .or. flag .eq. 1) then
        !     write(*,*) ISTEP, i, "DEBUG D ", "L =", l_obukhov
        !     flag = 1
        !   endif
        ! endif

        ! In case the iteration diverges we will just use this
        l_backup = l_obukhov

        l_old = 0
        count = 0
        if (abs(rib).lt.0.01) then ! neutral (use log law computed above)
          ! if (i.eq.1) then
          !   write(*,*) "Neutral", rib
          ! endif
          l_obukhov = 0
        elseif (rib.lt.-0.01) then ! convective
          ! if (i.eq.1) then
          !   write(*,*) "Convective", rib
          ! endif
          do while ((abs(l_old - l_obukhov)/abs(l_obukhov) .gt. 1e-3)
     $             .and. (count .lt. 20))

            l_old = l_obukhov
            count = count + 1

            ! for the central diff for evaluating dfdl
            fd_h = 1e-3*l_obukhov
            l_upper = l_obukhov + fd_h
            l_lower = l_obukhov - fd_h
            if (wmles_forcing_type .eq. surface_temperature) then
              f=(rib - h/l_obukhov
     $          *similarity_law_q_conv(l_obukhov, h, z1)
     $          /similarity_law_u_conv(l_obukhov, h, z0)**2)
              dfdl = (-h/l_upper*similarity_law_q_conv(l_upper, h, z1)
     $             /similarity_law_u_conv(l_upper, h, z0)**2)
              dfdl=dfdl + (h/l_lower
     $             *similarity_law_q_conv(l_lower, h, z1)
     $             /similarity_law_u_conv(l_lower, h, z0)**2)
              dfdl = dfdl/(2*fd_h)
            elseif (wmles_forcing_type .eq. surface_heat_flux) then
              f = (rib - h/l_obukhov/
     $        similarity_law_u_conv(l_obukhov, h, z0)**3)

              dfdl = (-h/l_upper
     $               /similarity_law_u_conv(l_upper, h, z0)**3)
              dfdl = dfdl + (h/l_lower/
     $        similarity_law_u_conv(l_lower, h, z0)**3)
              dfdl = dfdl/(2*fd_h)
            endif

            l_obukhov = l_obukhov - f/dfdl

            ! This is an adhoc upper bound for L, at which point we
            ! consider N-R to be diverged
            if (abs(l_obukhov) .gt. 20000 .or.
     $        abs(l_obukhov) .lt. 1e-5) then
              count = 20
            end if
          enddo
        else ! stable
          a = 5.0
          b = 5.0
          ! if (i.eq.1) then
          !   write(*,*) "Stable", rib
          ! endif
          do while ((abs(l_old - l_obukhov)/abs(l_obukhov) .gt. 1e-3)
     $           .and. (count .lt. 20))

            l_old = l_obukhov
            count = count + 1
            if (wmles_forcing_type .eq. surface_temperature) then
    !           if (ISTEP .gt. 373) then
    !             write(*,*) ISTEP, i, "LOOP A, rib =", rib,"L =",
    !  $            l_obukhov
    !  $           ,"sim_q =",similarity_law_q_stable(l_obukhov, h, z1)
    !  $           ,"sim_u =",similarity_law_u_stable(l_obukhov, h, z0)
    !             endif
              f = rib - h/l_obukhov*
     $              similarity_law_q_stable(l_obukhov, h, z1)/
     $              similarity_law_u_stable(l_obukhov, h, z0)**2
              ! if (ISTEP .gt. 373) then
              !   write(*,*) ISTEP, i, "LOOP B, f =", f
              ! endif
              dfdl = ((h*log(h/z0)*(2*a*h - b*h + l_obukhov*log(h/z0)))/
     $            (b*h + l_obukhov*log(h/z0))**3)
              ! if (ISTEP .gt. 373) then
              !   write(*,*) ISTEP, i, "LOOP C, dfdl =", dfdl
              !   write(*,*) ISTEP, i, "LOOP D, count =", count
              ! endif
            elseif (wmles_forcing_type .eq. surface_heat_flux) then
              write(*,*) "Not implemented yet!"
              call exitt
            endif
            l_obukhov = l_obukhov - f/dfdl


            ! This is an adhoc upper bound for L, at which point we
            ! consider N-R to be diverged
            if (abs(l_obukhov) > 20000 ) then
              count = 20
            end if
          enddo

        ! if (ISTEP .gt. 373) then
        !   if (flag .eq. 1) then
        !     write(*,*) ISTEP, i, "DEBUG E ", "f =",f,"dfdl =",dfdl
        !     flag = 1
        !   endif
        ! endif
        ! if (ISTEP .gt. 373) then
        !   if (abs(l_obukhov) .gt. 1000 .or. flag .eq. 1) then
        !     write(*,*) ISTEP, i, "DEBUG F ", "L =", l_obukhov
        !     flag = 1
        !   endif
        ! endif
        endif


        ! if we did not converge
        if (count .eq. 20) then
          ! write(*,*) "Unconverged :("
          l_obukhov = l_backup
        endif

        ! store the computed Obukhov length and Richardson number
        wmles_lobukhov(i) = l_obukhov
        wmles_ri(i) = rib

    !     if (ISTEP .gt. 373) then
    !       if (abs(l_obukhov) .gt. 1000 .or. abs(rib) .gt. 1
    !  $             .or. flag.eq.1) then
    !         write(*,*) ISTEP, i,"DEBUG G ","L =",l_obukhov,"rib =",rib
    !         flag = 1
    !       endif
    !     endif

        ! if the case is neutral the previously calculated
        ! values will be used without correction
        if (rib.lt.-0.01) then ! convective
          ! compute u* with the new obukhov length
          utau = kappa*magvh/similarity_law_u_conv(l_obukhov, h, z0)

          ! compute the surface heat flux if the temperature is prescribed
          if (wmles_forcing_type .eq. surface_temperature) then
            q = kappa*utau*(ts - th)
     $          /similarity_law_q_conv(l_obukhov, h, z1)
          endif
        elseif (rib.gt.0.01) then ! stable
          ! compute u* with the new obukhov length
          utau = kappa*magvh/similarity_law_u_stable(l_obukhov, h, z0)

    !       if (ISTEP .gt. 373) then
    !         if (abs(utau) .gt. 1 .or. flag .eq. 1) then
    !           write(*,*) ISTEP, i,"DEBUG H ","utau =", utau,
    !  $         "magvh =", magvh, "sim_u =",
    !  $          similarity_law_u_stable(l_obukhov, h, z0)
    !           flag = 1
    !         endif
    !       endif

          ! compute the surface heat flux if the temperature is prescribed
          if (wmles_forcing_type .eq. surface_temperature) then
            q = kappa*utau*(ts - th)
     $          /similarity_law_q_stable(l_obukhov, h, z1)

    !         if (ISTEP .gt. 373) then
    !           if (abs(q) .gt. 1 .or. flag .eq. 1) then
    !             write(*,*) ISTEP, i,"DEBUG I ","q =", q, "utau =",utau,
    !  $                "ts =", ts, "th =", th, "sim_q =",
    !  $                similarity_law_q_stable(l_obukhov, h, z1)
    !           endif
    !         endif
          endif
        endif

      endif
      ! if (ISTEP .gt. 373) then
      !   write(*,*) " "
      ! endif
      wmles_ustar(i) = utau
      wmles_local_index(i) = i
      ! if (abs(q) .gt. 5) then
      !   write(*,*) "DEBUG q =", q, "ustar = ", utau,
      ! endif
      ! write(*,*) "ustar =", wmles_ustar(i)

      ! Assign tau proportional to the velocity magnitudes at
      ! the sampling point
      wmles_tau(i, 1) = -utau**2*wmles_solh(i, 1)/magvh
      wmles_tau(i, 2) = 0
      wmles_tau(i, 3) = -utau**2*wmles_solh(i, 3)/magvh
      wmles_q(i) = q
      end

!--- convective --------------------------------------------------------
!> @brief Compute correction for the u log law in the convective case
      real function correction_u_conv(z, l)
      implicit none

      real z, l, xi, pi
!-----------------------------------------------------------------------

      pi = 4*atan(1.0)
      xi = (1.0 - 16.0*z/l)**0.25
      correction_u_conv = 2*log(0.5*(1 + xi)) + log(0.5*(1 + xi**2)) -
     $             2*atan(xi) + pi/2

      end function

!> @brief Compute correction for the q log law in the convective case
      real function correction_q_conv(z, l)
      implicit none

      real z, l, xi, pi
!-----------------------------------------------------------------------

      pi = 4*atan(1.0)
      xi = (1.0 - 16.0*z/l)**0.25
      correction_q_conv = 2*log(0.5*(1 + xi**2))

      end function

!> @brief Compute the similarity law for velocity
      real function similarity_law_u_conv(l_obukhov, h, z0)
      implicit none

      real l_obukhov, h, z0
      real correction_u_conv
!-----------------------------------------------------------------------

      similarity_law_u_conv = log(h/z0)
     $                        - correction_u_conv(h, l_obukhov)
     $                        + correction_u_conv(z0, l_obukhov)

      end function

!> @brief Compute the similarity law for heat
      real function similarity_law_q_conv(l_obukhov, h, z1)
      implicit none

      real l_obukhov, h, z1
      real correction_q_conv
!-----------------------------------------------------------------------

      similarity_law_q_conv = log(h/z1)
     $                        - correction_q_conv(h, l_obukhov)
     $                        + correction_q_conv(z1, l_obukhov)

      end function

!--- Stable ------------------------------------------------------------
!> @brief Compute correction for the u log law in the stable case
      real function correction_u_stable(z, l)
      implicit none

      real z, l
!-----------------------------------------------------------------------

      correction_u_stable = -5*z/l
      end function

!> @brief Compute correction for the u log law in the convective case
      real function correction_q_stable(z, l)
      implicit none

      real z, l
!-----------------------------------------------------------------------

      correction_q_stable = -5*z/l
      end function

!> @brief Compute the similarity law for velocity
      real function similarity_law_u_stable(l_obukhov, h, z0)
      implicit none

      real l_obukhov, h, z0
      real correction_u_stable
!-----------------------------------------------------------------------

      similarity_law_u_stable = log(h/z0)
     $                          - correction_u_stable(h, l_obukhov)
c     $                             + correction_u(z0, l_obukhov)

      end function

!> @brief Compute the similarity law for heat
      real function similarity_law_q_stable(l_obukhov, h, z1)
      implicit none

      real l_obukhov, h, z1
      real correction_q_stable
!-----------------------------------------------------------------------

      similarity_law_q_stable = log(h/z1)
     $                          - correction_q_stable(h, l_obukhov)
c     $                             + correction_q(z1, l_obukhov)


      end function
