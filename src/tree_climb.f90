! Compiling in Windows from command line (with Rtools installed):
! cd <working directory>
! "C:\Program Files\R\R-3.3.1\bin\R.exe" CMD SHLIB tree_climb.f90

! see R code for an explanation of parameters

recursive subroutine tree_climb(n, n_events, leaves, changes, changer, nodes, times, time, a,&
  edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, trait, sigma, ws, traits)

  implicit none
  integer, intent(IN) :: n, n_events, time, n_samples, n_leaves, ml
  integer, dimension(n_events - 1), intent(IN) :: changes
  integer, dimension(n_events), intent(IN) :: leaves, changer, nodes
  integer, dimension(n_samples), intent(IN) :: se

  integer :: a, place, newtime, tmpws, ws, ows
  integer, dimension(n + n_events - 2, 2) :: edge

  double precision, intent(IN) :: sigma
  double precision, dimension(n_events), intent(IN) :: times
  double precision, dimension(n_samples), intent(IN) :: t_el

  double precision, parameter :: pi = 3.1415926535879
  double precision :: time_in, tmptrait, trait, otrait, sumtime
  double precision, dimension(2) :: rands
  double precision, dimension(n + n_events - 2) :: edge_length, edge_trait
  double precision, dimension(ml, n_samples) :: samples

  logical :: traits

  ! create copies for recursive subroutine call
  ows = ws ! which sample we're currently at (at start of function call)
  otrait = trait ! current value of trait (at start of function call)
  
  time_in = 0.0 ! amount of time past the last event
  ! which branch we're on at the current time 
  ! (note: every instance of the function always starts at a speciation)
  ! to start, we pick the branch with smaller ID but later move to the other branch
  
  place = changer(time)
  call init_random_seed() ! set random seed (note that this currently cannot be user set)

  if (ws <= n_samples .and. traits) then ! if we haven't gone past the last sample yet
    do while (time == se(ws)) ! if we have one or more samples between this event and the next
      call random_number(rands) ! generate some random numbers
      ! add Brownian motion for the amount of time between either the last event
      ! or the last sample and the time of the current sample
      trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
      samples(place, ws) = trait ! place trait in the samples matrix
      time_in = t_el(ws) ! store how far past the last event we are
      ws = ws + 1 ! move to the next sample
      if (ws > n_samples) exit ! skip out of the loop if there are no more samples
    enddo
  endif
  
  newtime = time + 1 ! move to the next event in the tree
  sumtime = times(newtime) - time_in ! amount of time since last sample (or event)
  time_in = 0.0
  
  if (newtime < n_events) then ! if we haven't hit the end of the tree
    do while (place .ne. changer(newtime)) ! keep going till our branch has an event
      ! if an event happens to a branch with a smaller ID than ours, push the ID along
      ! (add one if the event is a speciation, remove if it's an extinction)
      if (changer(newtime) < place) place = place + changes(newtime)
      
      time_in = 0.0
  
      ! look for samples, see comments on lines 41-55
      if (ws <= n_samples .and. traits) then
        do while (newtime == se(ws))
          call random_number(rands)
          ! time since last trait update = sumtime (time from last sample (or starting event if none)
          ! to last event, or zero if a sample has happened since last event) 
          ! + t_el(ws) - time_in (time between last event/sample and current)
          trait = trait + sigma * sqrt(-2.0*(sumtime + t_el(ws) - time_in)*&
            log(rands(1)))*cos(2.0*pi*rands(2))
          sumtime = 0.0 ! reset sumtime as a sample has happened since last event
          samples(place, ws) = trait
          time_in = t_el(ws)
          ws = ws + 1
          if (ws > n_samples) exit
        enddo
      endif
      
      newtime = newtime + 1
      sumtime = sumtime + times(newtime) - time_in ! amount of time since last sample (or event)
      if (newtime == n_events) exit ! if we hit the end of the tree
    enddo
  endif
  
  ! we've now hit either the end of the tree or an event on our branch
  edge_length(a) = sum(times((time + 1):newtime)) ! total time between start and end points for edge
  if (traits) then
    call random_number(rands)
    ! time since last trait update = sumtime
    trait = trait + sigma * sqrt(-2.0*sumtime*log(rands(1)))*cos(2.0*pi*rands(2))
    edge_trait(a) = trait ! store final trait value at same point in vector as the matching edge
  endif
    
  if (newtime < n_events) then ! if we haven't reached the end of the tree
    ! use the pre-defined node IDs at the relevant events
    edge(a,:) = (/ nodes(time), nodes(newtime) /) 
    if (changes(newtime) == 1) then ! if this is a speciation
      a = a + 1 ! move to the next edge
      tmpws = ws ! store current values of ws and trait
      tmptrait = trait
      ! recursively run tree_climb on the split
      call tree_climb(n, n_events, leaves, changes, changer, nodes, times, newtime, a, edge,&
      edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, tmptrait, sigma, tmpws, traits)
    endif ! if changes != 1 then it's an extinction and we can stop
  else ! we've reached the end of the tree, so this is an extant node
    edge(a,:) = (/ nodes(time), place /) ! the end ID is simply which branch ID we are at
    ! now look for any remaining samples at or past the end of the tree
!   if (ws <= n_samples) then ! if we haven't gone past the last sample yet
!     time_in = 0.0 ! amount of time past the last event
!     do while (n_events == se(ws)) ! for all events at/past end of tree
!       call random_number(rands) ! generate some random numbers
!       trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
!       samples(place, ws) = trait ! place trait in the samples matrix
!       time_in = t_el(ws) ! store how far past the last event we are
!       ws = ws + 1 ! move to the next sample
!       if (ws > n_samples) exit ! skip out of the loop if there are no more samples
!     enddo
!   endif
  endif
  
  ! now we repeat the process with the other branch of the split
  
  a = a + 1 ! move to the next edge

  ! return to the original place, in terms of which sample is next
  ! and what the trait value is
  ws = ows 
  trait = otrait
  time_in = 0.0
  place = changer(time) + 1 ! now our branch ID at the event is the adjacent one to before

  ! repeat entire process again, look at comments for lines 41-109

  if (ws <= n_samples .and. traits) then
    do while (time == se(ws))
      call random_number(rands)
      trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
      samples(place, ws) = trait
      time_in = t_el(ws)
      ws = ws + 1
      if (ws > n_samples) exit
    enddo
  endif
  
  newtime = time + 1
  sumtime = times(newtime)
  time_in = 0.0

  if (newtime < n_events) then
    do while (place .ne. changer(newtime))
      if (changer(newtime) < place) place = place + changes(newtime)
      time_in = 0.0

      if (ws <= n_samples .and. traits) then  
        do while (newtime == se(ws))
          call random_number(rands)
          trait = trait + sigma * sqrt(-2.0*(sumtime + t_el(ws) - time_in)*&
            log(rands(1)))*cos(2.0*pi*rands(2))
          sumtime = 0.0
          samples(place, ws) = trait
          time_in = t_el(ws)
          ws = ws + 1
          if (ws > n_samples) exit          
        enddo
      endif

      newtime = newtime + 1
      sumtime = sumtime + times(newtime) - time_in
      if (newtime == n_events) exit
    enddo
  endif

  edge_length(a) = sum(times((time + 1):newtime))
  if (traits) then
    call random_number(rands)
    trait = trait + sigma * sqrt(-2.0*sumtime*log(rands(1)))*cos(2.0*pi*rands(2))
    edge_trait(a) = trait
  endif

  if (newtime < n_events) then
    edge(a,:) = (/ nodes(time), nodes(newtime) /)
    if (changes(newtime) == 1) then
      a = a + 1
      tmpws = ws
      tmptrait = trait
      call tree_climb(n, n_events, leaves, changes, changer, nodes, times, newtime, a, edge,&
      edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, tmptrait, sigma, tmpws, traits)
    endif
  else
    edge(a,:) = (/ nodes(time), place /) 
    ! now look for any remaining samples at or past the end of the tree
!   if (ws <= n_samples) then ! if we haven't gone past the last sample yet
!     time_in = 0.0 ! amount of time past the last event
!     do while (n_events == se(ws)) ! for all events at/past end of tree
!       call random_number(rands) ! generate some random numbers
!       trait = trait + sigma * sqrt(-2.0*(t_el(ws) - time_in)*log(rands(1)))*cos(2.0*pi*rands(2))
!       samples(place, ws) = trait ! place trait in the samples matrix
!       time_in = t_el(ws) ! store how far past the last event we are
!       ws = ws + 1 ! move to the next sample
!       if (ws > n_samples) exit ! skip out of the loop if there are no more samples
!     enddo
!   endif
  endif

end subroutine


! uses GCC GNU example routine for initialising random seed
! note that this is currently not settable by user
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed
