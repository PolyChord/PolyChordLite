!> This module contains all the tools to calculate evidence estimators
module evidence_module
    implicit none

    contains

    !> This take the existing estimate of 
    function update_Z_dead(evidence,likelihood,nlive,ndead)
        implicit none
        double precision, intent(in) :: evidence   !> Current evidence value to be updated
        double precision, intent(in) :: likelihood !> Likelihood of the late point
        integer,          intent(in) :: nlive      !> Number of live points
        integer,          intent(in) :: ndead      !> Number of dead points

        double precision :: update_Z_dead

        update_Z_dead = evidence + average_dead_volume(nlive,ndead)
        
    end function update_Z_dead


    function average_dead_volume(i,nlive)
        implicit none
        integer,intent(in) :: nlive      !> Number of live points
        integer,intent(in) :: i          !> index of point
        double precision :: average_dead_volume

        average_dead_volume = (1d0 + 1d0/nlive)**(-i)

    end function average_dead_volume


    !> This function gives the average prior volume contained within the iso-likelihood contour defined by the
    !! \f$\ell^\mathrm{th}\f$ live point:
    !! \f{eqnarray*}{ 
    !!   \left\langle X^{(\ell)}_\mathrm{live}\right\rangle 
    !!   &= \left\langle X^{(N)}_\mathrm{dead} \right\rangle \frac{(n-\ell+1)(n-\ell+2)}{(n+1)(n+2)} \\
    !!   &= {\left(1+\frac{1}{n}\right)}^{-N}                \frac{(n-\ell+1)(n-\ell+2)}{(n+1)(n+2)} \\
    !! \f}
    !! For more details on how this is calculated 
    !!

    function average_live_volume(ell,nlive,ndead)
        implicit none
        integer,intent(in) :: nlive      !> Number of live points
        integer,intent(in) :: ndead      !> Number of dead points
        integer,intent(in) :: ell        !> index of point

        double precision :: average_live_volume

        average_live_volume = &
            average_dead_volume(ndead,nlive) * (nlive-ell+1d0)/(nlive+1d0)

    end function average_live_volume



    function correlation_live_volumes(ell,m,nlive,ndead)
        implicit none
        integer,intent(in) :: nlive      !> Number of live points
        integer,intent(in) :: ndead      !> Number of dead points
        integer,intent(in) :: ell,m      !> indices of points

        double precision :: correlation_live_volumes

        correlation_live_volumes = correlation_dead_volumes(ndead,ndead,nlive,ndead)
        if ( ell >= m ) then
            correlation_live_volumes = correlation_live_volumes &
                * (nlive-ell+1d0) * (nlive-m+2d0) &
                / (nlive+1d0)     / (nlive+2d0)
        else
            correlation_live_volumes = correlation_live_volumes &
                * (nlive-m+1d0) * (nlive-ell+2d0) &
                / (nlive+1d0)   / (nlive+2d0)
        endif



    end function correlation_live_volumes


    function correlation_dead_volumes(i,j,nlive,ndead)
        implicit none
        integer,intent(in) :: nlive      !> Number of live points
        integer,intent(in) :: ndead      !> Number of dead points
        integer,intent(in) :: i,j        !> indices of points

        double precision :: correlation_dead_volumes



    end function correlation_dead_volumes


end module evidence_module
