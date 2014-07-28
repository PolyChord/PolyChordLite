module likelihoods_module
    implicit none


    contains

    function gaussian_loglike(theta,nDims)
        double precision, intent(in), dimension(nDims) :: theta
        integer, intent(in) :: nDims

        double precision :: gaussian_loglike

        double precision, dimension(nDims) :: sigma  
        double precision, dimension(nDims) :: mu

        double precision :: TwoPi 

        integer :: i


        TwoPi = 6.2831853d0
        sigma = 0.01
        mu    = 0.5

        gaussian_loglike = - nDims / 2d0 * log( TwoPi )

        do i = 1, nDims
            gaussian_loglike = gaussian_loglike - log( sigma(i) )
        enddo

        gaussian_loglike = gaussian_loglike - sum( ( ( theta( 1:nDims ) - mu(1:nDims ) ) / sigma( 1:nDims ) ) ** 2d0 ) / 2d0

    end function gaussian_loglike

end module likelihoods_module
