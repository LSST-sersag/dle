!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             
!    The Z-transformed Discrete Correlation Function algorithm (ZDCF)      
!                                                                             
!        Version 2.2 (fortran 95 double precision, 24/2/2013)
!
!                           Tal Alexander
!
!                   Weizmann Institute of Science
!                         Faculty of Physics
!             Dept. of Particle Physics and Astrophysics
!                   (tal.alexander@weizmann.ac.il)
!
!   Reference:                                                                
!                                                                             
!   T. Alexander, 1997, 'Is AGN variability correlated with other AGN
!   properties? - ZDCF analysis of small samples of sparse light
!   curves', in "Astronomical Time Series", eds. D. Maoz, A. Sternberg
!   and E.M. Leibowitz, Kluwer, Dordrecht, p. 163
!                                                                             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             
! SAMPLE RUN (user input marked by "<-USER". See explanations below)          
! ----------                                                                  
!                                                                             
! ZDCF V2.2 begins.                                                           
! Auto-correlation or cross-correlation? (1/2):                               
! 2     <-USER                                                                
! Enter output files prefix:                                                  
! dbg   <-USER                                                                
! Uniform sampling of light curve? (y/n):                                     
! n     <-USER                                                                
! Enter minimal number of points per bin (0 for default):                     
! 0     <-USER                                                                
! Omit zero-lag points? (y/n):                                                
! y     <-USER                                                                
! How many Monte Carlo runs for error estimation?                             
! 100   <-USER                                                                
! Enter name of 1st light curve file:                                         
! inp1  <-USER                                                                
! Enter name of 2nd light curve file:                                         
! inp2  <-USER                                                                
!                                                                             
! dbg.lc1  written (contains   59 points) ...                                 
! dbg.lc2  written (contains   63 points) ...                                 
!                                                                             
! =========================================================================== 
! ZDCF PARAMETERS:                                                            
! Autocorrelation?    F                                                       
! Uniform sampling?   F                                                       
! Omit zero lags?     T                                                       
! Minimal # in bin:   11                                                      
! # of Monte Carlo: 100 Monte Carlo seed: 123 59 points in inp1 63
! points in inp2
! ===========================================================================
!                                                                             
! Binning with minimum of  11 points per bin and resolution of 1.00E-03 .     
!                                                                             
!   71 bins actually used,   307 inter-dependent pairs discarded.             
!                                                                             
!   tau       -sig(tau)  +sig(tau)   dcf        -err(dcf)  +err(dcf) (#bin)   
!                                                                             
!  -5.346E+01  2.233E+01  4.374E-01  4.480E-01  3.394E-01  2.459E-01 (  10)   
!  -5.104E+01  9.915E-01  1.058E-01  6.378E-01  2.524E-01  1.636E-01 (  11)   
!   .                                                      .                  
!   .                                                      .                  
!   .                                                      .                  
!   5.203E+01  1.396E-01  9.683E-01  5.681E-02  3.433E-01  3.304E-01 (  11)   
!   5.498E+01  1.020E+00  2.047E+01  1.553E-01  3.704E-01  3.313E-01 (  10)   
!                                                                             
! dbg.dcf written...                                                          
! Program ended.                                                              
!                                                                             
! EXPLANATION                                                                 
! -----------                                                                 
!                                                                             
! o The input light curve files should be in 3 columns (free, floating        
!   point format):                                                            
!   time (ordered); flux / magnitude; absolute error on flux / magnitude.     
!                                                                             
! o The sign of the time lag is defined as                                    
!   tau = t(2nd light curve) - t(1st light curve)                             
!                                                                             
! o The program creates 3 output files with the same prefix. dbg.lc1, dbg.lc2 
!   are the 2 light curves (in the same format as the input) after points     
!   with identical times are averaged. dbg.dcf is the resulting ZDCF in 7     
!   columns, as in the example above (note: the extreme lags may be based     
!   on less than the minimal 11 points. these are the "left-over" points and  
!   are used just for giving an "impression" of the trends at extreme lags.   
!   However, these ZDCF points are NOT statistically reliable).               
!                                                                             
! o If the light curve sampling IS uniform, it is possible to force the       
!   program not to collect different time lags into a single bin              
!                                                                             
! o The default minimal number of points per bin is 11.                       
!                                                                             
! o It is recommended to omit the zero-lag points.                            
!                                                                             
! o 100 Monte Carlo runs seem to be sufficient in most cases.                 
!                                                                             
! FORTRAN ISSUES                                                              
! --------------                                                              
!                                                                             
! o This implementation is designed for speed at the cost of large
!   memory requirements.
!                                                                             
! o To economize work areas, the code is currently limited to light
!   curves with at most N=10^4 data points each. If your data set is
!   larger, increase the parameter LOG10N=4 in the module ZDCF_VARS
!   (in these comments, FORTRAN names and variables are capitalized).
!
! o The small parameter epsilon (see section 2.2.2 in the extended
!   description of the ZDCF, arXiv:1302.1508), which controls the
!   binning of pairs with close lags, is currently not user-
!   adjustable. It can be changed by modifying the parameter 
!   RSLUTN=0.001 in module ZDCF_VARS:
!   epsilon = RSLUTN * (The maximal time lag in the data).  
!                                                                             
! o The seed for the Monte Carlo random simulations, ISEED=123456, is
!   initialized in the module ZDCF_VARS. Change it to get an
!   independent error estimate.
!                                                                             
! VERSION LOG                                                                 
! -----------                                                                 
!                                                                             
! V1.0 ??/??/95: Original version                                             
! V1.1 09/05/96: Minor corrections for compatibility with Alpha/OSF           
! V1.2 21/10/99: Use open source routines for random numbers and sorting      
! V2.0 05/02/13: Complete rewrite and streamlining in double precision Fortran95 
! V2.1 15/02/13: Fixed index bug in the allocation of w%i in ALLOCATE_WORK; 
!                Economized work areas by using smaller integer type
!                for w%i, controlled by the parameter log10N; Updated
!                the "FORTRAN ISSUES" section in the comments; Fixed
!                error in reading 2nd LC in READ_OBS.
! V2.2 24/02/13: Further reduction of memory requirements in the
!                allocation of work areas wa,wb (maxpts decreased from
!                a*n*b%n to min(a%n,b%n))
!                                                                             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE cons
  IMPLICIT NONE 
  ! Defining double precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 60)
  REAL(dp), PARAMETER :: Pi = 4*ATAN(1.0_dp)
  !
  CHARACTER(len=4), PARAMETER :: zdcf_ver = 'V2.2'
  CHARACTER(len=50), PARAMETER :: maintainer_email = &
       'tal.alexander@weizmann.ac.il'
END MODULE cons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE zdcf_vars
  USE cons
  IMPLICIT NONE
  ! Input parameters
  LOGICAL :: autocf,no0lag,uni_sample
  INTEGER :: nMC
  INTEGER :: minpts,maxpts
  INTEGER, PARAMETER :: enough=11
  REAL(dp), PARAMETER :: rslutn = 0.001
  INTEGER :: nseeds,iseed = 123456
  ! The observed light curves
  TYPE lc_type
     INTEGER :: n
     REAL(dp), ALLOCATABLE :: t(:),x(:),err(:) ! n values
     REAL(dp), ALLOCATABLE :: MC(:)            ! n values
     LOGICAL, ALLOCATABLE :: used(:)           ! n values
  END TYPE lc_type
  ! The binned correlation function
  TYPE dcf_type
     INTEGER :: nbins
     REAL(dp), ALLOCATABLE :: t(:),sigtm(:),sigtp(:) ! nbins values
     REAL(dp), ALLOCATABLE :: r(:),sigrm(:),sigrp(:) ! nbins values
     REAL(dp), ALLOCATABLE :: avz(:),expz(:),sigz(:) ! nbins values
     INTEGER, ALLOCATABLE :: inbin(:)                ! nbins values
     INTEGER(kind=8) :: used,unused
  END TYPE dcf_type
  ! Work areas for clcdcf
  ! Economizing storage by limiting work areas to light curves 
  ! with no more than 10**log10N points each. Currently, log10N=4.
  ! Change as needed.
  INTEGER, PARAMETER :: log10N = 4
  INTEGER, PARAMETER :: iw_kind = SELECTED_INT_KIND(log10N)
  !
  TYPE work_type
     REAL(dp), ALLOCATABLE :: t(:),x(:),err(:)    ! maxpts values
     INTEGER(kind=iw_kind), ALLOCATABLE :: i(:,:) ! maxpts,nbins values
  END TYPE work_type
  TYPE(work_type) :: wa,wb
END MODULE zdcf_vars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Functions for calculating the ZDCF
MODULE zdcf_funcs
  USE cons
  USE zdcf_vars
  IMPLICIT NONE 
CONTAINS 
!//////////////////////////////////////////////////////////////////////////////
! Generating a "true" signal by substracting gaussian noise from observations.
! The error errx is the ABSOLUTE error in x.
    SUBROUTINE simerr(lc)
      USE cons
      USE zdcf_vars
      IMPLICIT NONE
      TYPE(lc_TYPE), INTENT(inout) :: lc
      !
      lc%MC = lc%x-lc%err*rndnrm(lc%MC)
    END SUBROUTINE simerr
    !///////////////////////////////////////////////////////////////////////
    ! Generating an array (size(a)) of standard Gaussian deviates using the 
    ! Box-Muller method.
    FUNCTION rndnrm(a)
      USE cons
      IMPLICIT NONE
      REAL(dp), DIMENSION(:),INTENT(in) :: a
      REAL(dp), DIMENSION(SIZE(a)) :: rndnrm
      !
      REAL(dp), DIMENSION((SIZE(a)+1)/2) :: x1,x2,C
      !
      INTEGER :: n,n2
      REAL(dp), PARAMETER :: pi2 = 2*Pi
      !
      n = size(a)
      n2 = (n+1)/2
      !
      CALL RANDOM_NUMBER(x1)
      CALL RANDOM_NUMBER(x2)
      x2 = x2*pi2
      C =  SQRT(-2.*LOG(x1))
      x1 = C*COS(x2)
      x2 = C*SIN(x2)
      rndnrm(1:n2) = x1
      rndnrm(n2+1:n) = x2(1:n-n2)
    END FUNCTION rndnrm
  !//////////////////////////////////////////////////////////////////////////////
  ! Fisher's small sample approximation for E(z) (Kendall + Stuart Vol. 1 p.391)    
  ELEMENTAL FUNCTION fishe(r,n)
    USE cons
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: r
    INTEGER, INTENT(in) :: n
    REAL(dp) :: fishe
    !
    fishe = LOG((1+r)/(1-r))/2.+ &
         r/2/(n-1)*(1+(5+r**2)/4/(n-1)+ &
         (11+2*r**2+3*r**4)/8/(n-1)**2)
  END FUNCTION fishe
  !//////////////////////////////////////////////////////////////////////////////
  ! Fisher's small sample approximation for s(z) (Kendall + Stuart Vol. 1 p.391)
  ELEMENTAL FUNCTION fishs(r,n)
    USE cons
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: r
    INTEGER, INTENT(in) :: n
    REAL(dp) :: fishs
    !
    fishs = 1.0_dp/(n-1)*(1+(4-r**2)/2./(n-1)+ &
         (22-6*r**2-3*r**4)/6/(n-1)**2)
    fishs = SQRT(MAX(0.0_dp,fishs))
  END FUNCTION fishs
  !//////////////////////////////////////////////////////////////////////////////
  ! Allocating points to bins
  SUBROUTINE alcbin(a,b,dcf)
    USE cons
    USE zdcf_vars
    IMPLICIT NONE
    TYPE(lc_type), INTENT(inout) :: a,b
    TYPE(dcf_type), INTENT(out) :: dcf
    ! Local dynamical work areas
    INTEGER, ALLOCATABLE :: idx(:),waidx(:),wbidx(:)
    REAL(dp), ALLOCATABLE :: wtau(:)
    !
    INTEGER :: nbins,mbins
    INTEGER :: i,np,pfr,pmax,incr,nnegtv,plo,phi,ierr,inb,j,p
    REAL(dp) :: tij,tolrnc
    LOGICAL :: first_time = .TRUE. 
    !
    ! Allocating the dynamic work areas
    !
    ALLOCATE (waidx(a%n*b%n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate waidx(',a%n*b%n,') Error ',ierr
       STOP
    END IF
    ALLOCATE (wbidx(a%n*b%n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate wbidx(',a%n*b%n,') Error ',ierr
       STOP
    END IF
    ALLOCATE (wtau(a%n*b%n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate wtau(',a%n*b%n,') Error ',ierr
       STOP
    END IF
    ALLOCATE (idx(a%n*b%n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate idx(',a%n*b%n,') Error ',ierr
       STOP
    END IF
    !
    ! Calculate all the time lag points
    !
    IF (first_time) THEN
       first_time = .FALSE.
       WRITE (*, &
            '(''Binning with minimum of '',i3,'// &
            ''' points per bin and resolution of '',1p,g8.2,'' .'')') &
            minpts,rslutn
    END IF
    np = 0
    IF (autocf) THEN
       DO i = 1,a%n
          acf_loop: DO j = i,b%n
             tij = b%t(j)-b%t(i)
             IF (no0lag .AND. tij == 0.0_dp) CYCLE acf_loop
             np = np+1
             wtau(np) = tij
             waidx(np) = i
             wbidx(np) = j
          END DO acf_loop
       END DO
    ELSE
       DO i = 1,a%n
          ccf_loop: DO j = 1,b%n
             tij = b%t(j)-a%t(i)
             IF (no0lag .AND. tij == 0.0_dp) CYCLE ccf_loop
             np = np+1
             wtau(np) = tij
             waidx(np) = i
             wbidx(np) = j
          END DO ccf_loop
       END DO
    END IF
    !
    ! Sort index according to increasing time lag
    !
    CALL heap_sort_index(np,wtau,idx,sort=.FALSE.)
    !
    ! calculating the tolerance level for lags to be considered the same
    !
    tij = wtau(idx(np))-wtau(idx(1))
    tolrnc = tij*rslutn
    !
    ! The dcf is allocated at this with the maximal number of bins
    ! possible for the requested minpts. The actual number used is
    ! subsequently determined on the fly and is typically much smaller
    !
    mbins = int((np+1.0_dp)/minpts)+1
    CALL allocate_dcf(dcf,mbins)
    !
    ! maxpts is the maximal number of points in a bin. In principle
    ! maxpts = np ( the number of all pairs: a huge number). In
    ! practice, for non uniform sampling it is not supposed to be much
    ! larger than minpts (sometimes there are a few extra points per
    ! bin if they are closer to the bin's boundary than the resolution
    ! tolerance).  For uniform sampling, it is bound by min(a%n,b%n)
    ! To be on the safe side, min(a%n,b%n) is assumed for both cases.
    !
    !old: maxpts = a%n*b%n
    maxpts = MIN(a%n,b%n)
    ! WARNING! wa, wb can be HUGE! 
    CALL allocate_work(wa,maxpts,mbins)
    CALL allocate_work(wb,maxpts,mbins)
    !
    ! Looping on bins
    !
    nbins = 0
    !
    ! If binned CCF: binning from median time-lag upwards and backwards!
    !
    IF (autocf .OR. uni_sample) THEN
       pfr = 1
       pmax = np
       incr = 1
    ELSE
       pfr = np/2
       pmax = 1
       incr = -1
    END IF
    !
    bin_loop: DO
       inb = 0
       !bug?: tij = incr*huge(0.0_dp)
       tij = wtau(idx(pfr))
       nbins = nbins+1
       ! This shouldn't happen...
       IF (nbins > mbins) THEN
          WRITE (*,*) 'ALCBIN: nbins > mbins (Send bug report to ',&
               maintainer_email,'!'
          STOP          
       END IF
       dcf%t(nbins) = 0.0
       !
       ! Initialize used point flag vectors
       !
       a%used = .FALSE.; b%used = .FALSE. 
       !
       ! Collect points into bins that contain at least "enough" points,
       ! but do not break points with the same lag (up to the tolerance)
       ! into separate bins.
       !
       tau_loop: DO i = pfr,pmax,incr
          p = idx(i)
          ! Check whether bin is full
          IF ((ABS(wtau(p)-tij) > tolrnc .AND. & ! tij = previous lag
               (inb >= minpts .OR. uni_sample)) .OR. &          
               (i == pmax)) THEN
             !
             ! Bin is full: Calculating tau and its std (before
             ! proceeding to the next bin)
             !
             dcf%inbin(nbins) = inb
             dcf%t(nbins) = dcf%t(nbins)/inb
             IF (uni_sample) THEN
                dcf%sigtm(nbins) = 0.0
                dcf%sigtp(nbins) = 0.0
                pfr = i
                ! If not enough points in bin, ignore it
                IF (inb < minpts) nbins = nbins-1
                IF (pfr /= pmax) CYCLE bin_loop
             ELSE
                ! If the last point is alone in its bin, we ignore it (to
                ! avoid subsequent divisions by zero) and get immediately
                ! out of tau loop.
                IF (inb <= 1) THEN
                   nbins = nbins-1
                   EXIT tau_loop
                END IF
                !
                ! Finding the 0.3413 (+/-1sig) time lags above and below
                ! the bin mean time lag
                !
                IF (incr == 1) THEN
                   plo = pfr
                   phi = i-1
                ELSE
                   plo = i+1
                   phi = pfr
                END IF
                midpnt_loop: DO p = plo,phi
                   IF (wtau(idx(p)) >= dcf%t(nbins)) THEN
                      j = p ! The point closest to the bin's mean
                      EXIT midpnt_loop
                   END IF
                END DO midpnt_loop
                !??? Why the factor 2???
                p = MAX(j-NINT(dble(j-plo)*0.3413_dp*2),plo)
                dcf%sigtm(nbins) = MAX(dcf%t(nbins)-wtau(idx(p)),0.0_dp)
                p = MIN(j+NINT(dble(phi-j)*0.3413_dp*2),phi)
                dcf%sigtp(nbins) = MAX(wtau(idx(p))-dcf%t(nbins),0.0_dp)
                pfr = i
                IF (pfr /= pmax) CYCLE bin_loop
             END IF
             ! If no more points - get out of tau loop.
             EXIT tau_loop
          END IF
          !
          ! Adding another point to the bin...
          !
          ! Avoiding correlated pairs
          IF ((.NOT. a%used(waidx(p))) .AND. (.NOT. b%used(wbidx(p)))) THEN
             inb = inb+1
             ! This shouldn't happen...
             IF (inb > maxpts) THEN
                WRITE (*,*) 'ALCBIN: inb > maxpts = ',maxpts
                WRITE (*,*) 'Send bug report to ',maintainer_email,'!)'
                STOP
             END IF
             a%used(waidx(p)) = .TRUE.; b%used(wbidx(p)) = .TRUE.
             tij = wtau(p)
             dcf%t(nbins) = dcf%t(nbins)+tij
             wa%i(inb,nbins) = INT(waidx(p),iw_kind)
             wb%i(inb,nbins) = INT(wbidx(p),iw_kind)
          END IF
       END DO tau_loop
       !
       ! Binning is finished
       !
       IF ( .NOT. (autocf .OR. uni_sample .OR. incr == 1)) THEN
          ! Now, go back and bin the other half lag axis
          pfr = np/2+1
          pmax = np
          incr = 1
          nnegtv = nbins
          CYCLE bin_loop
       END IF
       EXIT bin_loop
    END DO bin_loop
    !
    ! If CCF (and NOT uniform sampling): Sort the bins into increasing
    ! chronological order: The nnegtv negative bins are at the
    ! beginning but at reverse order.
    !
    IF ( .NOT. (autocf .OR. uni_sample)) THEN
       DO i = 1,nnegtv/2
          j = nnegtv+1-i
          CALL swap_i(dcf%inbin(i),dcf%inbin(j))
          CALL swap_r(dcf%t(i),dcf%t(j))
          CALL swap_r(dcf%sigtp(i),dcf%sigtp(j))
          CALL swap_r(dcf%sigtm(i),dcf%sigtm(j))
          DO p = 1,MAX(dcf%inbin(i),dcf%inbin(j))
             CALL swap_iw(wa%i(p,i),wa%i(p,j))
             CALL swap_iw(wb%i(p,i),wb%i(p,j))             
          END DO
       END DO
    END IF
    !
    dcf%nbins = nbins
    !
    ! Deallocating the local work areas
    !
    DEALLOCATE (waidx)
    DEALLOCATE (wbidx)
    DEALLOCATE (wtau)
    DEALLOCATE (idx)
    !
    ! This shouldn't happen...
    !
    IF (dcf%nbins == 0) THEN
       WRITE (*,*) 'ALCBIN: nbins = 0 (Send bug report to ',&
            maintainer_email,'!)'
       STOP
    END IF
  END SUBROUTINE alcbin
  !/////////////////////////////////////////////////////////////////////////////
  ! Swap 2 integer variables
  SUBROUTINE swap_i(x,y)
    USE cons
    IMPLICIT NONE 
    INTEGER, INTENT(inout) :: x,y
    !
    INTEGER :: z
    z = x; x = y; y = z
  END SUBROUTINE swap_i
  !/////////////////////////////////////////////////////////////////////////////
  ! Swap 2 integer(kind=iw_kind) variables
  SUBROUTINE swap_iw(x,y)
    USE cons
    IMPLICIT NONE 
    INTEGER(kind=iw_kind), INTENT(inout) :: x,y
    !
    INTEGER(kind=iw_kind) :: z
    z = x; x = y; y = z
  END SUBROUTINE swap_iw
  !/////////////////////////////////////////////////////////////////////////////
  ! Swap 2 double precision real variables
  SUBROUTINE swap_r(x,y)
    USE cons
    IMPLICIT NONE 
    REAL(dp), INTENT(inout) :: x,y
    !
    REAL(dp) :: z
    z = x; x = y; y = z
  END SUBROUTINE swap_r
  !/////////////////////////////////////////////////////////////////////////////
  ! Calculating the discrete correlation function.
  ! POSITIVE lag values mean b lags after a.
  ! This implementation requires rather big work areas in the interest of a
  ! faster algorithm (and is therefore suitable for Monte Carlo simulations).
  SUBROUTINE clcdcf(a,b,dcf,new_data)
    USE cons
    USE zdcf_vars
    IMPLICIT NONE 
    TYPE(lc_type), INTENT(inout) :: a,b
    TYPE(dcf_type), INTENT(inout) :: dcf
    LOGICAL, INTENT(in) :: new_data
    !
    INTEGER :: n,ibin
    REAL(dp) :: expa,expb,vara,varb,vnorm,expbin,z,expz,sigz
    REAL(dp), PARAMETER :: eps = 1d-7
    !
    ! If new data (i.e. NOT another Monte Carlo run) - allocate pairs to bins
    !
    IF (new_data) THEN
       ! Allocating the lags to the bins
       CALL alcbin(a,b,dcf)
       dcf%used = SUM(dcf%inbin(1:dcf%nbins))
       IF (autocf) THEN
          IF (no0lag) THEN
             dcf%unused = a%n*(a%n-1)/2-dcf%used
          ELSE
             dcf%unused = a%n**2/2-dcf%used
          END IF
       ELSE
          dcf%unused = a%n*b%n-dcf%used
       END IF
    END IF 
    ! After allocating pairs to bins: calculating the dcf
    DO ibin = 1,dcf%nbins
       ! Collecting the points of the bin
       n = dcf%inbin(ibin) ! Assumed >1
       IF (n < 2) CYCLE ! This shouldn't happen...
       wa%t = a%t(wa%i(:n,ibin))
       wb%t = b%t(wb%i(:n,ibin))
       wa%x = a%x(wa%i(:n,ibin))
       wb%x = b%x(wb%i(:n,ibin))
       wa%err = a%err(wa%i(:n,ibin))
       wb%err = b%err(wb%i(:n,ibin))
       !
       expa = SUM(wa%x(:n))/n;        expb = SUM(wb%x(:n))/n
       vara = SUM(wa%x(:n)**2);       varb = SUM(wb%x(:n)**2); 
       ! Dividing by (n-1) for an unbiased estimator of the
       ! correlation coefficient cf Barlow / Statistics, p. 80
       vara = (vara-n*expa**2)/(n-1); varb = (varb-n*expb**2)/(n-1)       
       vnorm = vara*varb
       IF (vnorm <= 0.0) THEN
          ! Pathological case: normalization factor <= 0
          dcf%r(ibin) = 0.0_dp
       ELSE
          expbin = SUM((wa%x(:n)-expa)*(wb%x(:n)-expb))
          expbin = expbin/SQRT(vnorm)/(n-1)
          ! Making sure -1 < r < 1
          expbin = MAX(MIN(expbin,1.0_dp-eps),-1.0_dp+eps) 
          dcf%r(ibin) = expbin
       END IF
       !
       ! Calculating the +/- 1 Sigma limits from Fisher's z
       !
       ! NOTE: This error estimation is by "bootstrapping": fishe & fishs give
       !       the true E(z) and S(z) when the TRUE correlation coefficient is
       !       given. We are using the empirical r itself, similarily to the
       !       common poissonian estimate of n +/- sqrt(n)
       !
       z = LOG((1+expbin)/(1-expbin))/2
       sigz = fishs(expbin,n)
       expz = fishe(expbin,n)
       dcf%sigrm(ibin) = expbin-TANH(expz-sigz)
       dcf%sigrp(ibin) = TANH(expz+sigz)-expbin       
    END DO
  END SUBROUTINE clcdcf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! Create an ascending sort index for an array (and optionally sort
  ! it) using the n*log(n) heap sort algorithm.  
  !
  ! Usage: The first n values in real(dp) array are used to create the
  ! permutation index index (does not need to be initialized) for an
  ! ascending order sort. If optional logical parameter sort is given
  ! as true, array is sorted as well.
  !
  ! Adapted from 
  ! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/sort3_f90.txt
  SUBROUTINE heap_sort_index(n,array,index,sort)
    USE cons
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    REAL(dp), INTENT(inout) :: array(n)
    INTEGER, INTENT(out) :: index(n)
    LOGICAL, OPTIONAL, INTENT(in) :: sort
    !
    REAL(dp), DIMENSION(SIZE(array)) :: a 
    INTEGER :: i,j,l,ia,iaa
    REAL(dp) :: aa
    !
    a = array
    index = [(i,i=1,n)]
    !
    l = n/2+1
    ia = n
    ! The index L will be decremented from its initial value during the
    ! "hiring" (heap creation) phase. Once it reaches 1, the index ia 
    ! will be decremented from its initial value down to 1 during the
    ! "retirement-and-promotion" (heap selection) phase.
    loop: DO
       IF (l > 1) THEN
          l = l-1
          aa = a(l)
          iaa = index(l)
       ELSE
          aa = a(ia)
          iaa = index(ia)
          a(ia) = a(1)
          index(ia) = index(1)
          ia = ia-1
          IF (ia == 1) THEN
             a(1) = aa
             index(1) = iaa
             EXIT loop
          END IF
       END IF
       i = l; j = l+l
       !
       DO WHILE (j <= ia)
          IF (j < ia) THEN
             IF (a(j) < a(j+1)) j = j+1
          END IF
          IF (aa < a(j)) THEN
             a(i) = a(j)
             index(i) = index(j)
             i = j; j = j+j
          ELSE
             j = ia+1
          END IF
       END DO
       a(i) = aa
       index(i) = iaa
    END DO loop
    IF (PRESENT(sort) .AND. sort) array = a
  END SUBROUTINE heap_sort_index
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate the arrays in a type(lc) variable (light curve data)
  SUBROUTINE allocate_lc(lc,n)
    USE cons
    USE zdcf_vars
    IMPLICIT NONE
    TYPE(lc_type), INTENT(out) :: lc
    INTEGER, INTENT(in) :: n
    !
    INTEGER :: ierr
    !
    ALLOCATE(lc%t(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate lc%t(',n,') Error ',ierr
       STOP
    END IF
    ALLOCATE(lc%x(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate lc%x(',n,') Error ',ierr
       STOP
    END IF
    ALLOCATE(lc%err(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate lc%err(',n,') Error ',ierr
       STOP
    END IF
    ALLOCATE(lc%used(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate lc%used(',n,') Error ',ierr
       STOP
    END IF
    ALLOCATE(lc%MC(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate lc%MC(',n,') Error ',ierr
       STOP
    END IF
  END SUBROUTINE allocate_lc
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Deallocate the arrays in a type(lc) variable
  SUBROUTINE deallocate_lc(lc)
    USE zdcf_vars
    IMPLICIT NONE 
    TYPE(lc_type), INTENT(inout) :: lc
    !
    DEALLOCATE(lc%t)
    DEALLOCATE(lc%x)
    DEALLOCATE(lc%err)
    DEALLOCATE(lc%used)
    DEALLOCATE(lc%MC)
  END SUBROUTINE deallocate_lc
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! allocate the arrays in a type(dcf) variable (zdcf results)
  SUBROUTINE allocate_dcf(dcf,n)
    USE cons
    USE zdcf_vars
    IMPLICIT NONE 
    TYPE(dcf_type) :: dcf
    INTEGER :: n
    !
    INTEGER :: ierr
    !
    ALLOCATE(dcf%t(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%t(',n,') Error ',ierr
       STOP
    END IF
    dcf%t = 0.0_dp
    ALLOCATE(dcf%sigtm(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%sigtm(',n,') Error ',ierr
       STOP
    END IF
    dcf%sigtm = 0.0_dp
    ALLOCATE(dcf%sigtp(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%sigtp(',n,') Error ',ierr
       STOP
    END IF
    dcf%sigtp = 0.0_dp
    ALLOCATE(dcf%r(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%r(',n,') Error ',ierr
       STOP
    END IF
    dcf%r = 0.0_dp
    ALLOCATE(dcf%sigrm(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%sigrm(',n,') Error ',ierr
       STOP
    END IF
    dcf%sigrm = 0.0_dp
    ALLOCATE(dcf%sigrp(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%sigrp(',n,') Error ',ierr
       STOP
    END IF
    dcf%sigrp = 0.0_dp
    ALLOCATE(dcf%inbin(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%inbin(',n,') Error ',ierr
       STOP
    END IF
    dcf%inbin = 0
    ALLOCATE(dcf%avz(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%avz(',n,') Error ',ierr
       STOP
    END IF
    dcf%avz = 0.0_dp
    ALLOCATE(dcf%expz(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%expz(',n,') Error ',ierr
       STOP
    END IF
    dcf%expz = 0.0_dp
    ALLOCATE(dcf%sigz(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate dcf%sigz(',n,') Error ',ierr
       STOP
    END IF
    dcf%sigz = 0.0_dp
  END SUBROUTINE allocate_dcf
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Deallocate the arrays in a type(dcf) variable
  SUBROUTINE deallocate_dcf(dcf)
    USE zdcf_vars
    IMPLICIT NONE
    TYPE(dcf_type), INTENT(inout) :: dcf
    !
    DEALLOCATE(dcf%t)
    DEALLOCATE(dcf%sigtm)
    DEALLOCATE(dcf%sigtp)
    DEALLOCATE(dcf%r)
    DEALLOCATE(dcf%sigrm)
    DEALLOCATE(dcf%sigrp)
    DEALLOCATE(dcf%inbin)
    DEALLOCATE(dcf%avz)
    DEALLOCATE(dcf%expz)
    DEALLOCATE(dcf%sigz)
  END SUBROUTINE deallocate_dcf
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate the arrays in a type(work) variable (work areas for binning)
  ! n is the maximal number of data pairs used for the zdcf
  ! m is the maximal number of bins in the zdcf
  SUBROUTINE allocate_work(w,n,m)
    USE cons
    USE zdcf_vars
    IMPLICIT NONE
    TYPE(work_type), INTENT(out) :: w
    INTEGER, INTENT(in) :: n,m
    !
    INTEGER :: ierr
    !
    ALLOCATE(w%t(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate w%t(',n,') Error ',ierr
       STOP
    END IF
    ALLOCATE(w%x(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate w%x(',n,') Error ',ierr
       STOP
    END IF
    ALLOCATE(w%err(n),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate w%err(',n,') Error ',ierr
       STOP
    END IF
    ALLOCATE(w%i(n,m),stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'Cannot allocate w%i(',n,m,') Error ',ierr
       STOP
    END IF
  END SUBROUTINE allocate_work
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Deallocate the arrays in a type(work) variable
  SUBROUTINE deallocate_work(w)
    USE zdcf_vars
    IMPLICIT NONE 
    TYPE(work_type), INTENT(inout) :: w
    !
    DEALLOCATE (w%t)
    DEALLOCATE (w%x)
    DEALLOCATE (w%err)
    DEALLOCATE (w%i)
  END SUBROUTINE deallocate_work
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE zdcf_funcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM zdcf
  USE cons
  USE zdcf_vars
  USE zdcf_funcs
  IMPLICIT NONE
  !
  TYPE(lc_type) :: a,b
  TYPE(dcf_type) :: dcf
  !
  CHARACTER(len=100) :: infil1,infil2,outfil,prefix
  INTEGER :: i,oerr,cerr,lprfix
  REAL(dp), PARAMETER :: eps=1d-7
  LOGICAL :: new_data
  !
  INTEGER :: inpi
  CHARACTER(1) :: inpc
  !
  ! Reading program parameters
  !
  WRITE (*,'("ZDCF ",a4," begins:")') zdcf_ver
  WRITE (*,'("Auto-correlation or cross-correlation? (1/2):")')
  READ *,inpi
  autocf = (inpi == 1)
  
  WRITE (*,'("Enter output files prefix:")')
  READ (*,*) prefix
  lprfix = LEN(TRIM(prefix))
  
  WRITE (*,'("Uniform sampling of light curve? (y/n):")')
  READ (*,'(a1)') inpc
  uni_sample = (inpc .EQ. 'y' .OR. inpc .EQ. 'Y')
  WRITE (*,'("Enter minimal number of points per bin (0 for default):")')
  READ (*,*) minpts
  ! if minpts = 0, use default
  IF (minpts <= 0) minpts = enough
  
  WRITE (*,'("Omit zero-lag points? (y/n):")')
  READ (*,'(a1)') inpc
  no0lag = (inpc .EQ. 'y' .OR. inpc .EQ. 'Y')
  
  WRITE (*,'("How many Monte Carlo runs for error estimation?")')
  READ (*,*) nMC
  IF (nMC <= 1) nMC = 0
  !
  WRITE (*,*)
  WRITE (*,'(75("="))')
  WRITE (*,'("ZDCF PARAMETERS:")')
  WRITE (*,'("Autocorrelation?  ",l)') autocf
  WRITE (*,'("Uniform sampling? ",l)') uni_sample
  WRITE (*,'("Omit zero lags?   ",l)') no0lag
  WRITE (*,'("Minimal # in bin: ",i10)') minpts
  WRITE (*,'("# of Monte Carlo: ",i10)') nMC
  WRITE (*,'("Monte Carlo seed: ",i10)') iseed
 
  !
  ! Reading the observed data
  !
  WRITE (*,'("Enter name of 1st light curve file:")')
  READ (*,*) infil1
  IF (.NOT. autocf) THEN
     WRITE (*,'("Enter name of 2nd light curve file:")')
     READ (*,*) infil2
  END IF
  !
  CALL read_obs(a,infil1,b,infil2) 
  !
  ! Writing the condensed light curve data
  !
  WRITE (*,*)
  outfil = prefix(:lprfix)//'.lc1'
  OPEN (unit=11,file=outfil,status='UNKNOWN',iostat=oerr)
  IF (oerr /= 0) THEN
     WRITE (*,*) 'Error ',oerr,' opening file ',TRIM(outfil)
     STOP
  END IF
  WRITE (11,'(1p,3(1x,e12.5))') (a%t(i),a%x(i),a%err(i),i=1,a%n)
  CLOSE (unit=11,status='KEEP',iostat=cerr)
  IF (cerr /= 0) THEN
     WRITE (*,*) 'Error ',oerr,' closing file ',TRIM(outfil)
     STOP
  END IF
  WRITE  (*,*) outfil(:lprfix+5),' written (contains ',a%n,' points) ...'
  
  outfil = prefix(:lprfix)//'.lc2'
  OPEN (unit=12,file=outfil,status='UNKNOWN',iostat=oerr)
  IF (oerr /= 0) THEN
     WRITE (*,*) 'Error ',oerr,' opening file ',TRIM(outfil)
     STOP
  END IF
  WRITE (12,'(1p,3(1x,e12.5))') (b%t(i),b%x(i),b%err(i),i=1,b%n)
  CLOSE (unit=12,status='KEEP',iostat=cerr)
  IF (cerr /= 0) THEN
     WRITE (*,*) 'Error ',oerr,' closing file ',TRIM(outfil)
     STOP
  END IF
  WRITE (*,*) outfil(:lprfix+5),' written (contains ',b%n,' points) ...'
  
  !old: WRITE (*,'(i5," points in ",a)') a%n,TRIM(infil1)
  !old: IF (.NOT. autocf) WRITE (*,'(i5," points in ",a)') b%n,TRIM(infil2)
  WRITE (*,'(75("="))')
  WRITE (*,*)
  !
  ! Estimating the effects of the measurement errors 
  ! by Monte Carlo simulations
  !
  CALL RANDOM_SEED(size=nseeds)
  CALL RANDOM_SEED(put=[(iseed+i*43,i=1,nseeds)])
  !
  IF (nMC > 1) THEN
     !
     ! Calculate the ZDCF with Monte Carlo errors.
     !      
     new_data = .TRUE.
     MC: DO i = 1,nMC
        IF (autocf) THEN
           CALL simerr(a)
        ELSE
           CALL simerr(a)
           CALL simerr(b)
        END IF
        CALL clcdcf(a,b,dcf,new_data)
        new_data = .FALSE. 
        IF (i == 1) THEN
           ! initialize avz (after dcf is allocated in clcdcf)
           dcf%avz = 0.0_dp            
           WRITE (*,*)
           WRITE (*,'(i10," bins actually used, ",i10,&
                &" inter-dependent pairs discarded.")') &
                dcf%nbins,dcf%unused
        END IF
        !
        ! The summing and averaging is done in z-space.
        !
        ! Making sure -1 < dcf < 1 in all bins
        dcf%r = MAX(MIN(dcf%r,1.0_dp-eps),-1.0_dp+eps) 
        dcf%avz = dcf%avz+LOG((1+dcf%r)/(1-dcf%r))/2
     END DO MC
     !
     dcf%avz = dcf%avz/nMC
     dcf%r = TANH(dcf%avz)
     dcf%r = MAX(MIN(dcf%r,1.0_dp-eps),-1.0_dp+eps) 
     dcf%avz = LOG((1+dcf%r)/(1-dcf%r))/2
     dcf%expz = fishe(dcf%r,dcf%inbin)     
     dcf%sigz = fishs(dcf%r,dcf%inbin)
     dcf%sigrm = dcf%r-TANH(dcf%expz-dcf%sigz)
     dcf%sigrp = TANH(dcf%expz+dcf%sigz)-dcf%r
  ELSE
     !
     ! calculate the ZDCF w/o Monte Carlo errors.
     !
     new_data = .TRUE.
     CALL clcdcf(a,b,dcf,new_data)        
     WRITE (*,*)
     WRITE (*,'(i5," bins actually used, ",i5,&
          &" inter-dependent pairs discarded.")') &
          dcf%nbins,dcf%unused
  END IF
  !
  ! printing the results: tdcf, dcf, sigdcf
  !
  outfil = prefix(:lprfix)//'.dcf'
  OPEN (unit=13,file=outfil,status='UNKNOWN',iostat=oerr)
  IF (oerr /= 0) THEN
     WRITE (*,*) 'Error ',oerr,' opening file ',TRIM(outfil)
     STOP
  END IF
  WRITE (*,*)
  WRITE (*,'(1p,6(1x,a10),'' ('',a4,'')'')') &
       ' tau      ','-sig(tau) ','+sig(tau) ',' dcf      ', &
       ' -err(dcf)',' +err(dcf)','#bin'
  WRITE (*,*)
  DO i = 1,dcf%nbins
     WRITE (*,'(1p,6(1x,e10.3),'' ('',i4,'')'')') &
          dcf%t(i),dcf%sigtm(i),dcf%sigtp(i), &
          dcf%r(i),dcf%sigrm(i),dcf%sigrp(i),dcf%inbin(i)
     WRITE (13,'(1p,6(1x,e10.3),1x,i4)') &
          dcf%t(i),dcf%sigtm(i),dcf%sigtp(i), &
          dcf%r(i),dcf%sigrm(i),dcf%sigrp(i),dcf%inbin(i)
  END DO
  CLOSE (unit=13,status='KEEP',iostat=cerr)
  IF (cerr /= 0) THEN
     WRITE (*,*) 'Error ',oerr,' closing file ',TRIM(outfil)
     STOP
  END IF
  WRITE (*,*)
  WRITE (*,*) outfil(:lprfix+4),' written...'
  !
  CALL deallocate_lc(a);    CALL deallocate_lc(b) 
  CALL deallocate_dcf(dcf)
  CALL deallocate_work(wa); CALL deallocate_work(wb)
  !
  !WRITE (*,'(/,"Please cite the ZDCF as:")') 
  !WRITE (*,&
  !     '("T. Alexander, 1997, Is AGN variability correlated with other AGN",/,&
  !       &"properties? - ZDCF analysis of small samples of sparse light",/,&
  !       &"curves, in Astronomical Time Series, eds. D. Maoz, A. Sternberg",/,&
  !       &"and E.M. Leibowitz, Kluwer, Dordrecht, p. 163",/)')
  !
  WRITE (*,'(/,"ZDCF ended.")')
CONTAINS 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Reading the observed points
  SUBROUTINE read_obs(a,infil1,b,infil2)
    USE cons
    USE zdcf_vars
    IMPLICIT NONE 
    TYPE(lc_type) :: a,b
    CHARACTER(len=*),INTENT(in) :: infil1,infil2
    !
    INTEGER :: n,nsamet,oerr,cerr
    REAL(dp) :: t,t0,x,err
    !
    ! Find the size of the 2 data sets and allocate the areas
    !
    OPEN (unit=1,file=infil1,status='UNKNOWN',iostat=oerr)
     IF (oerr /= 0) THEN
       WRITE (*,*) 'READ_OBS: Error ',oerr,' opening file ',TRIM(infil1)
       STOP
    END IF
    n = 0
    t0 = -HUGE(0.0_dp)
    DO 
       READ (1,*,END=11) t
       ! Counting equal time entries as a single datum
       IF (t > t0) THEN
          n = n+1
       ELSE IF (t < t0) THEN
          WRITE (*,*) 'READ_OBS: Mis-ordered times in file ',&
               &TRIM(infil1),' at line ',n
          STOP
       END IF
       t0 = t
    END DO
11  CONTINUE
    REWIND(1)
    !
    CALL allocate_lc(a,n)
    !
    IF (.NOT. autocf) THEN
       OPEN (unit=2,file=infil2,status='UNKNOWN',iostat=oerr)
       IF (oerr /= 0) THEN
          WRITE (*,*) 'READ_OBS: Error ',oerr,' opening file ',TRIM(infil2)
          STOP
       END IF       
       n = 0
       t0 = -HUGE(0.0_dp)
       DO 
          READ (2,*,END=22) t
          ! Counting equal time entries as a single datum
          IF (t > t0) THEN
             n = n+1
          ELSE IF (t < t0) THEN
             WRITE (*,*) 'READ_OBS: Mis-ordered times in file ',&
                  &TRIM(infil2),' at line ',n
             STOP
          END IF
          t0 = t
       END DO
22     CONTINUE
       REWIND(2)
    END IF
    !
    CALL allocate_lc(b,n)
    !
    ! Reading the 1st light curve
    !
    nsamet = 1
    n = 1
    t0 = -HUGE(0.0_dp)
    DO
       READ (1,*,END=111) t,x,err
       a%t(n) = t; a%x(n) = x; a%err(n) = err
       !
       ! Averaging observations with identical times
       !
       IF (n > 1 .AND. a%t(n) == t0) THEN
          a%x(n-1) = a%x(n-1)+a%x(n)
          a%err(n-1) = a%err(n-1)+a%err(n)
          nsamet = nsamet+1
       ELSE IF (n > 1) THEN
          a%x(n-1) = a%x(n-1)/nsamet
          a%err(n-1) = a%err(n-1)/nsamet
          nsamet = 1
          t0 = a%t(n)
          n = n+1
       ELSE
          n = n+1
       END IF
    END DO
111 CONTINUE
    CLOSE (unit=1,status='KEEP',iostat=cerr)
    IF (cerr /= 0) THEN
       WRITE (*,*) 'READ_OBS: Error ',oerr,' closing file ',TRIM(infil1)
       STOP
    END IF
    n = n-1
    a%n = n
    IF (a%n <= 0) THEN
       WRITE (*,*) 'READ_OBS: No data in ',infil1
       STOP
    END IF
    a%x(n) = a%x(n)/nsamet
    a%err(n) = a%err(n)/nsamet
    !
    ! Reading the 2nd light curve
    !
    IF (autocf) THEN
       b%n = a%n
       b%t = a%t
       b%x = a%x
       b%err = a%err
    ELSE
       nsamet = 1
       n = 1
       t0 = -HUGE(0.0_dp)
       DO
          READ (2,*,END=222) t,x,err
          b%t(n) = t; b%x(n) = x; b%err(n) = err
          !
          ! Averaging observations with identical times
          !
          IF (n > 1 .AND. b%t(n) == t0) THEN
             b%x(n-1) = b%x(n-1)+b%x(n)
             b%err(n-1) = b%err(n-1)+b%err(n)
             nsamet = nsamet+1
          ELSE IF (n > 1) THEN
             b%x(n-1) = b%x(n-1)/nsamet
             b%err(n-1) = b%err(n-1)/nsamet
             nsamet = 1
             t0 = b%t(n)
             n = n+1
          ELSE
             n = n+1
          END IF
       END DO
222    CONTINUE
       CLOSE (unit=2,status='KEEP',iostat=cerr)
       IF (cerr /= 0) THEN
          WRITE (*,*) 'READ_OBS: Error ',oerr,' closing file ',TRIM(infil2)
          STOP
       END IF
       n = n-1
       b%n = n
       IF (n <= 0) THEN
          WRITE (*,*) 'READ_OBS: No data in ',infil2
          STOP
       END IF
       b%x(n) = b%x(n)/nsamet
       b%err(n) = b%err(n)/nsamet
    END IF
  END SUBROUTINE read_obs
END PROGRAM zdcf
