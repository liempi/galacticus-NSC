!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a submodule which provides composite trapezoidal numerical integrators.
!!}

submodule (Numerical_Integration2) Numerical_Integration2_Trapezoidal
  !!{RST
  Provides composite trapezoidal numerical integrators.
  !!}
  implicit none

contains

  module procedure compositeTrapezoidalConstructor1D
    !!{RST
    Constructor for a one-dimensional, composite trapezoidal numerical integrator.
    !!}
    implicit none

    call iterationsMaximumValidate(iterationsMaximum)
    self%iterationsMaximum=iterationsMaximum
    return
  end procedure compositeTrapezoidalConstructor1D

  module procedure compositeTrapezoidalEvaluate1D
    !!{RST
    Evaluate a one-dimension integral using a numerical composite trapezoidal rule.
    !!}
    use :: Error    , only : Error_Report
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    integer :: iteration, j, n
    logical :: converged
    double precision :: cumulant, integrand, integrandPrevious, x, step

    cumulant =0.5d0*(self%integrand(b)+self%integrand(a))
    iteration=0
    converged=.false.
    do while (.not.converged)
       iteration=iteration+1
       if (iteration > self%iterationsMaximum)                    call Error_Report(                                                    'maximum iterations exceeded'//                      {introspection:location}                          )
       n   =2**iteration
       step=(b-a)/dble(n)
       ! Evaluate only every second point, as the others were computed in previous iterations.
       do j=1,n-1,2
          x=a+dble(j)*step
          cumulant=cumulant+self%integrand(x)
       end do
       integrand=step*cumulant
       if (iteration > 2)                                      converged=Values_Agree(                                                 integrand             ,                          integrandPrevious     ,                          self%toleranceAbsolute,                          self%toleranceRelative                          )
      integrandPrevious=integrand
   end do
    compositeTrapezoidalEvaluate1D=integrand
    return
  end procedure compositeTrapezoidalEvaluate1D

  module procedure adaptiveCompositeTrapezoidalConstructor1D
    !!{RST
    Constructor for a one-dimensional, adaptive composite trapezoidal numerical integrator.
    !!}
    implicit none

    self%integratorCompositeTrapezoidal1D=integratorCompositeTrapezoidal1D(iterationsMaximum)
    return
  end procedure adaptiveCompositeTrapezoidalConstructor1D

  module procedure adaptiveCompositeTrapezoidalEvaluate1D
    !!{RST
    Evaluate a one-dimension integral using a numerical composite trapezoidal rule.
    !!}
    use :: Error, only : Error_Report
    implicit none
    logical :: converged
    double precision :: error, width, midpoint, fa, fb, fm
    type            (interval                                ), pointer :: head, newInterval1, newInterval2, current

    ! Create our first estimate of the integral, using a single trapezoid.
    allocate(head)
    head%a        =                  a
    head%b        =                  b
    head%fa       =   self%integrand(a)
    head%fb       =   self%integrand(b)
    head%depth    =   1
    head%next     =>  null()
    head%integral =   0.5d0                         *(head%fb+head%fa)             *(head% b-head% a)
    head%error    =   head%integral
    ! Initialize current integral and error estimates, and flag that we're not converged.
    adaptiveCompositeTrapezoidalEvaluate1D=head%integral
    error                                 =head%integral
    converged                             =.false.
    ! Iterate until convergence is reached.
    do while (.not.converged)
       ! Bisect the head interval. By construction, this will always be the interval with the largest absolute error.
       current  => head      ! Pop the head from the list.
       head     => head%next
       if (current%depth > self%iterationsMaximum)                call Error_Report(                                                    'maximum iterations exceeded'//                      {introspection:location}                          )
       width    =  (current%b-current%a)/2.0d0
       midpoint =  (current%b+current%a)/2.0d0
       ! Evaluate the function at the endpoints and midpoints of the current interval. For endpoints, reuse stored values.
       fa       =  current%fa
       fb       =  current%fb
       fm       =  self%integrand(midpoint)
       ! Create two new subintervals.
       allocate(newInterval1)
       allocate(newInterval2)
       newInterval1%next     => null()
       newInterval2%next     => null()
       newInterval1%a        =  current%a
       newInterval1%b        =  midpoint
       newInterval2%a        =  midpoint
       newInterval2%b        =  current%b
       newInterval1%fa       =  fa
       newInterval1%fb       =  fm
       newInterval2%fa       =  fm
       newInterval2%fb       =  fb
       newInterval1%depth    =  current     %depth+1
       newInterval2%depth    =  newInterval1%depth
       ! Compute the integral and error estimate in each new subinterval.
       newInterval1%integral =  (fm+fa)*0.5d0*width
       newInterval2%integral =  (fm+fb)*0.5d0*width
       newInterval1%error    =  abs(0.5d0*width*fm-0.25d0*current%integral)
       newInterval2%error    =  newInterval1%error
       ! Update the total integral and error estimate for the new subintervals.
       adaptiveCompositeTrapezoidalEvaluate1D           =adaptiveCompositeTrapezoidalEvaluate1D   -current     %integral                    +newInterval1%integral                    +newInterval2%integral
       error                                            =error                                    -current     %error                       +newInterval1%error                       +newInterval2%error
       ! Test for convergence.
       converged=                                                                                  abs(error) <= self%toleranceAbsolute                                               .or.                                                                                 abs(error) <= self%toleranceRelative*abs(adaptiveCompositeTrapezoidalEvaluate1D)
       ! Destroy the old interval.
       deallocate(current)
       ! Insert the new intervals into our stack.
       ! First, link the two subintervals - they will always be consecutive in the list as they have the same error estimate.
       newInterval1%next => newInterval2
       call intervalInsert(head,newInterval1,newInterval2)
    end do
    ! Destroy the stack.
    call intervalListDestroy(head)
    return
  end procedure adaptiveCompositeTrapezoidalEvaluate1D

  module procedure vectorizedCompositeTrapezoidalConstructor1D
    !!{RST
    Constructor for a one-dimensional, vectorized composite trapezoidal numerical integrator.
    !!}
    implicit none
    integer :: workspaceSize, i, j, n
    double precision :: step

    call iterationsMaximumValidate(iterationsMaximum)
    self%iterationsMaximum=iterationsMaximum
    workspaceSize=2**iterationsMaximum
    allocate(self%d(workspaceSize))
    ! Points are stored in order of the level of refinement at which they are first used: those of level n=2ⁱ (the odd
    ! multiples of 1/n) are in elements n/2+1 to n. The end points, used before any refinement, are not stored, so element 1 is
    ! unused.
    self%d(1)=0.0d0
    do i=1,iterationsMaximum
       n=2**i
       step=1.0d0/dble(n)
       forall(j=1:n-1:2)
          self%d(n/2+1+j/2)=dble(j)*step
       end forall
    end do
    return
  end procedure vectorizedCompositeTrapezoidalConstructor1D

  module procedure vectorizedCompositeTrapezoidalEvaluate1D
    !!{RST
    Evaluate a one-dimension integral using a numerical vectorized composite trapezoidal rule.
    !!}
    use :: Error    , only : Error_Report
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    integer :: iteration, n
    logical :: converged
    double precision :: cumulant, integrand, integrandPrevious, step

    cumulant =0.5d0*sum(self%integrand([b,a]))
    iteration=0
    converged=.false.
    do while (.not.converged)
       iteration=iteration+1
       if (iteration > self%iterationsMaximum)                    call Error_Report(                                                    'maximum iterations exceeded'//                      {introspection:location}                          )
       n   =2**iteration
       step=(b-a)/dble(n)
       ! Evaluate only every second point, as the others were computed in previous iterations.
       cumulant =cumulant+sum(self%integrand(self%d(n/2+1:n)*(b-a)+a))
       integrand=step*cumulant
       if (iteration > 2)                                      converged=Values_Agree(                                                 integrand             ,                          integrandPrevious     ,                          self%toleranceAbsolute,                          self%toleranceRelative                          )
       integrandPrevious=integrand
    end do
    vectorizedCompositeTrapezoidalEvaluate1D=integrand
    return
  end procedure vectorizedCompositeTrapezoidalEvaluate1D

  module procedure multiVectorizedCompositeTrapezoidal1DConstructor
    !!{RST
    Constructor for a one-dimensional, multi-integrand, vectorized composite trapezoidal numerical integrator.
    !!}
    use :: Error, only : Error_Report
    implicit none

    if (intervalsMaximum < 1)                                                 call Error_Report(                                                                   'at least 1 intevals is required'//                                 {introspection:location}                                         )
    self%intervalsMaximum=intervalsMaximum
    return
  end procedure multiVectorizedCompositeTrapezoidal1DConstructor

  module procedure multiVectorizedCompositeTrapezoidal1DEvaluate
    !!{RST
    Evaluate a one-dimension integral using a numerical composite trapezoidal rule.
    !!}
    use            :: Error  , only : Error_Report, errorStatusFail       , errorStatusSuccess
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    implicit none
    integer         (c_size_t                                       ), parameter :: intervalListSizeInitial  =1000_c_size_t, bisectionIntegrandsMinimum=8_c_size_t, bisectionIntervalsMaximum=6000_c_size_t
    double precision                                                 , dimension(self%integrandCount  ) :: error, errorScale, errorsMaximum
    logical                                                          , dimension(self%integrandCount  ) :: converged, mustEvaluate, convergedPrevious
    double precision                                                 , dimension(                    2) :: xUnion
    double precision                                                 , dimension(self%integrandCount,2) :: fUnion
    double precision                                                 , dimension(self%integrandCount  ) :: fa, fb, fm
    type            (intervalMultiList                              ), dimension(:)                   , allocatable :: list
    double precision                                                 , dimension(:)                   , allocatable :: listValue
    integer         (c_size_t                                       ), dimension(:)                   , allocatable :: listRank
    type            (intervalMulti                                  ), pointer :: head, newInterval1, newInterval2, current, previous, lastInsertCurrent, lastInsertPrevious, searchStart, searchEnd, searchMidpoint
    double precision :: midpoint, width, errorMaximum
    integer         (c_size_t                                       ) :: intervalCount, i, indexStart, indexEnd, indexMidpoint
    logical :: precisionLost
    !$GLC attributes initialized :: previous

    ! If the interval has zero size, return a zero result.
    if (a == b) then
       integral=0.0d0
       if (present(status)) status=errorStatusSuccess
       return
    end if
    ! Create list of intervals.
    allocate(list     (intervalListSizeInitial))
    allocate(listValue(intervalListSizeInitial))
    allocate(listRank (intervalListSizeInitial))
    ! Create our first estimate of the integral, using a single interval.
    intervalCount=1_c_size_t
    head => intervalMultiNew(a,b,self%integrandCount)
    lastInsertCurrent => null()
    mustEvaluate      =  .true.
    xUnion=[head%a,head%b]
    call self%integrand(self%integrandCount,xUnion,mustEvaluate,fUnion)
    head%fa      =fUnion(:,1)
    head%fb      =fUnion(:,2)
    head%integral=+0.5d0                      *(head%fb+head%fa)          *(head% b-head% a)
    head%error   =abs(head%integral)
    ! Initialize current integral and error estimates, and flag that we're not converged.
    integral =head%integral
    error    =head%error
    converged=                                                    abs(error) <= self%toleranceAbsolute                 .or.                                                   abs(error) <= self%toleranceRelative*abs(integral)
    ! Iterate until convergence is reached.
    precisionLost=.false.
    do while (.not.all(converged) .and. intervalCount < self%intervalsMaximum .and. .not.precisionLost)
       ! Bisect the head interval. By construction, this will always be the interval with the largest absolute error.
       current       => head      ! Pop the head from the list.
       head          => head%next
       intervalCount =  intervalCount-1_c_size_t
       width         =  (current%b-current%a)/2.0d0
       midpoint      =  (current%b+current%a)/2.0d0
       if (associated(lastInsertCurrent,current)) lastInsertCurrent => null()
       ! Evaluate the function at the endpoints and midpoints of the current interval. For endpoints, reuse stored values.
       fa       =  current%fa
       fb       =  current%fb
       ! Create two new subintervals.
       newInterval1 => intervalMultiNew(current%a,midpoint ,self%integrandCount)
       newInterval2 => intervalMultiNew(midpoint ,current%b,self%integrandCount)
       ! Check for loss of precision in interval extent.
       if     (                         midpoint==current%a    .or.                     midpoint==current%b   ) then
          call intervalMultiFailureReport(a,b,current,head,integral,error,converged,present(status))
          deallocate(current     )
          deallocate(newInterval1)
          deallocate(newInterval2)
          ! Force exit.
          precisionLost=.true.
          exit          
       end if
       ! Evaluate the function at the midpoint of the current interval.
       mustEvaluate=.not.converged
       xUnion(1)=midpoint
       call self%integrand(self%integrandCount,xUnion(1:1),mustEvaluate,fUnion(:,1))
       fm=fUnion(:,1)
       ! Compute the integral and error estimate at the midpoint.
       convergedPrevious=converged
       ! In cases where the integrand was not evaluated at the midpoint (i.e. the integral is already converged), estimate the
       ! integrand at the midpoint by linear interpolation - which is consistent with the approximations of the trapezoidal
       ! integration rule.
       where (mustEvaluate)
          ! Set new intervals and errors, for integrands that were evaluated (any not evaluated are already converged so we do not
          ! need to update them further). Error on second integral is identical to that on the first.
          newInterval1%fa      =    fa
          newInterval2%fb      =    fb
          newInterval1%fb      = fm
          newInterval2%fa      = fm
          newInterval1%integral=(fm+fa)*0.5d0*width
          newInterval2%integral=(fm+fb)*0.5d0*width
          newInterval1%error   =abs(0.5d0*width*fm-0.25d0*current%integral)
          newInterval2%error   =newInterval1%error
          ! Updated integrals and error estimates.
          integral =+integral                   -current     %integral      +newInterval1%integral      +newInterval2%integral
          error    =+error                      -current     %error         +newInterval1%error         +newInterval2%error
       end where
       ! Test for convergence.
       converged= converged                                               .or.                                                      abs(error) <= self%toleranceAbsolute                    .or.                                                      abs(error) <= self%toleranceRelative*abs(integral)
       if (.not.all(converged) .and. any(converged .neqv. convergedPrevious) .and. intervalCount > 1_c_size_t) then
          ! One or more integrals have converged. We must resort the linked list of intervals to ensure that it is in order of
          ! descending error measure.
          errorScale=max(                                                self%toleranceAbsolute              ,           self%toleranceRelative*abs(integral)           )
          call intervalMultiResort(head,intervalCount,converged,errorScale,list,listValue,listRank)
       end if
       ! Destroy the old interval.
       deallocate(current)
       ! Insert the new intervals into our stack.
       ! First, link the two subintervals - they will always be consecutive in the list as they have the same error estimate.
       newInterval1%next => newInterval2
       ! Check if the stack is empty.
       if (associated(head)) then
          ! Stack is not empty. Perform an insertion sort to insert our new subinterval into the sorted stack. The sorting is
          ! such that the stack is ordered with the interval for which the scaled error (i.e. error divided by the larger of
          ! the absolute and relative tolerances) of non-converged integrals is descending.
          ! Find the maximum scaled error across all non-converged integrands in our new interval.
          errorMaximum=0.0d0
          do i=1,self%integrandCount
             if (converged(i)) then
                errorScale(i)=1.0d0
             else
                errorScale(i)=max(                                                         self%toleranceAbsolute(i)                 ,              self%toleranceRelative(i)*abs(integral(i))              )
                errorMaximum =max(errorMaximum,newInterval1%error(i)/errorScale(i))
             end if
          end do
          errorsMaximum=+errorScale            *errorMaximum
          ! Search for where to insert this interval in the list. We check if the insertion should happen after the location of
          ! the previous insertion. If it does, we begin searching from that point. This is often efficient as consecutively
          ! evaluated intervals often have similar errors.
          if (associated(lastInsertCurrent)) then
             if (any(.not.converged .and. lastInsertCurrent%error > errorsMaximum)) then
                current   => lastInsertCurrent
                previous  => lastInsertPrevious
             else
                current   => head
                previous  => head
             end if
          else
             current   => head
             previous  => head
          end if
          ! Determine whether to do a simple search through the list (stepping one interval at a time and checking the error
          ! measure), or doing a bisection search. The simple search requires more evaluations of the error measure, but the
          ! bisection search requires construction of an array of pointers to intervals. The optimal strategy depends on the
          ! number of integrands and intervals in use. Heuristically the following condition seems to give reasonably good
          ! performance.
          if (self%integrandCount <= bisectionIntegrandsMinimum .and. intervalCount > bisectionIntervalsMaximum) then
             ! Simple search.
             call intervalMultiSearch(current,previous,converged,errorsMaximum)
          else
             ! Bisection search.
             indexStart  = 1_c_size_t
             searchStart => current
             if (any(.not.converged .and. searchStart%error > errorsMaximum)) then
                call intervalMultiWorkspaceEnsure(list,listValue,listRank,intervalCount)
                indexEnd                 =  1_c_size_t
                searchEnd                => current
                list(indexEnd)%interval_ => head
                do while (associated(searchEnd%next))
                   indexEnd                 =  +indexEnd                              +1_c_size_t
                   previous                 =>  searchEnd
                   searchEnd                =>  searchEnd%next
                   list(indexEnd)%interval_ =>  searchEnd
                end do
                if (any(.not.converged .and. searchEnd%error > errorsMaximum)) then
                   ! Insert at end of list.
                   current => searchEnd
                else
                   ! Bisection search for insertion point.
                   do while (indexEnd > indexStart+1_c_size_t)
                      ! Find the midpoint.
                      indexMidpoint  =  indexStart+(indexEnd-indexStart)/2_c_size_t
                      searchMidPoint => list(indexMidpoint)%interval_
                      if (any(.not.converged .and. searchMidpoint%error > errorsMaximum)) then
                         indexStart  =  indexMidpoint
                         searchStart => searchMidpoint
                      else
                         indexEnd    =  indexMidpoint
                         searchEnd   => searchMidpoint
                      end if
                   end do
                   current  => searchEnd
                   previous => searchStart
                end if
             end if
          end if
          ! Record the insertion point for use next time.
          lastInsertCurrent  => current
          lastInsertPrevious => previous
          call intervalMultiLink(head,current,previous,newInterval1,newInterval2,converged,errorsMaximum)
       else
          ! Stack is empty - simply point the head to our new subinterval.
          head => newInterval1
       end if
       intervalCount=intervalCount+2_c_size_t
    end do
    ! Destroy the stack.
    call intervalMultiListDestroy(head)
    ! Destroy the interval pointer list.
    deallocate(list     )
    deallocate(listValue)
    deallocate(listRank )
    ! Report error if number of intervals was exceeded.
    if (present(status)) status=errorStatusSuccess
    if (precisionLost) then
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report("loss of precision in integration interval"//{introspection:location})
       end if
    else if (intervalCount >= self%intervalsMaximum .and. .not.all(converged)) then
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report('maximum number of intervals exceeded'//{introspection:location})
       end if
    end if
    return
  end procedure multiVectorizedCompositeTrapezoidal1DEvaluate

end submodule Numerical_Integration2_Trapezoidal
