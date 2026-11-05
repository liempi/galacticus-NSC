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
Contains a submodule which provides linked lists of integration intervals, maintained in order of descending error.
!!}

submodule (Numerical_Integration2) Numerical_Integration2_Intervals
  !!{RST
  Provides linked lists of integration intervals, maintained in order of descending error.
  !!}
  implicit none

contains

  module procedure intervalInsert
    !!{RST
    Insert a chain of new intervals, all sharing the error of ``newFirst``, into a list sorted by descending absolute error. The new intervals are placed ahead of any existing interval with an equal error.
    !!}
    implicit none
    type(interval)               , pointer :: current, previous

    if (.not.associated(head)) then
       head => newFirst
       return
    end if
    current  => head
    previous => null()
    do while (associated(current%next).and.abs(current%error) > abs(newFirst%error))
       previous => current
       current  => current%next
    end do
    if (.not.associated(current%next) .and. abs(current%error) > abs(newFirst%error)) then
       ! New intervals go at the end of the list.
       current%next     => newFirst
    else if (associated(previous)) then
       ! New intervals go before the current interval, which is not the head.
       newLast %next    => previous%next
       previous%next    => newFirst
    else
       ! New intervals go before the head.
       newLast %next    => head
       head             => newFirst
    end if
    return
  end procedure intervalInsert

  module procedure intervalListDestroy
    !!{RST
    Destroy a list of intervals.
    !!}
    implicit none
    type(interval)               , pointer :: current, previous

    current => head
    do while (associated(current))
       previous => current
       current  => current%next
       deallocate(previous)
    end do
    head => null()
    return
  end procedure intervalListDestroy

  module procedure intervalMultiSearch
    !!{RST
    Step along a list of multi-integrand intervals from ``current`` until reaching the first interval whose error, for every non-converged integrand, does not exceed ``errorsMaximum`` (or the final interval). On return ``previous`` is the interval preceding ``current``, unless ``current`` did not move.
    !!}
    implicit none

    do while (                                                                        associated(                     current%next                 )      .and.                                                                   any       (.not.converged .and. current%error > errorsMaximum)     )
       previous => current
       current  => current%next
    end do
    return
  end procedure intervalMultiSearch

  module procedure intervalMultiLink
    !!{RST
    Link a chain of new multi-integrand intervals into a list at the insertion point found by a search: before ``current``, unless ``current`` is the final interval and has a larger error than the new intervals.
    !!}
    implicit none

    if (.not.associated(current%next) .and. any(.not.converged .and. current%error > errorsMaximum)) then
       ! New intervals go at the end of the list.
       current %next => newFirst
    else if (associated(current,head)) then
       ! New intervals go before the head.
       newLast %next => head
       head          => newFirst
    else
       ! New intervals go before the current interval.
       newLast %next => previous%next
       previous%next => newFirst
    end if
    return
  end procedure intervalMultiLink

  module procedure intervalMultiNew
    !!{RST
    Allocate a new multi-integrand interval spanning ``a`` to ``b``, not linked into any list.
    !!}
    implicit none

    allocate(interval_                    )
    allocate(interval_%fa      (integrandCount))
    allocate(interval_%fb      (integrandCount))
    allocate(interval_%integral(integrandCount))
    allocate(interval_%error   (integrandCount))
    interval_%a    =  a
    interval_%b    =  b
    interval_%next => null()
    return
  end procedure intervalMultiNew

  module procedure intervalMultiWorkspaceEnsure
    !!{RST
    Ensure that the workspace arrays used for reordering a list of multi-integrand intervals can hold ``intervalCount`` intervals.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none

    if (allocated(list)) then
       if (intervalCount <= size(list)) return
       deallocate(list     )
       deallocate(listValue)
       deallocate(listRank )
    end if
    allocate(list     (2*intervalCount))
    allocate(listValue(2*intervalCount))
    allocate(listRank (2*intervalCount))
    return
  end procedure intervalMultiWorkspaceEnsure

  module procedure intervalMultiResort
    !!{RST
    Reorder a list of multi-integrand intervals by descending maximum scaled error over non-converged integrands. This is needed when an integrand converges, as this changes the measure by which the list is sorted.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Sorting      , only : sortIndex
    implicit none
    type            (intervalMulti    ), pointer :: interval_
    integer         (c_size_t         ) :: iInterval
    integer :: i

    call intervalMultiWorkspaceEnsure(list,listValue,listRank,intervalCount)
    interval_ => head
    do iInterval=1_c_size_t,intervalCount
       list(iInterval)%interval_ => interval_
       ! We avoid use of maxval(,mask=) here (and instead determine the maximum error using a simple loop) because maxval()
       ! evaluates all elements prior to applying the mask, so it wastes time in computing values that we will never need.
       listValue(iInterval)=0.0d0
       do i=1,size(converged)
          if (.not.converged(i)) listValue(iInterval)=max(listValue(iInterval),interval_%error(i)/errorScale(i))
       end do
       interval_ => interval_%next
    end do
    listRank(1_c_size_t:intervalCount)=sortIndex(listValue(1_c_size_t:intervalCount))
    head      => list(listRank(intervalCount))%interval_
    interval_ => head
    do iInterval=2_c_size_t,intervalCount
       interval_%next => list(listRank(intervalCount+1_c_size_t-iInterval))%interval_
       interval_      => interval_%next
    end do
    interval_%next => null()
    return
  end procedure intervalMultiResort

  module procedure intervalMultiListDestroy
    !!{RST
    Destroy a list of multi-integrand intervals.
    !!}
    implicit none
    type(intervalMulti)               , pointer :: current, previous

    current => head
    do while (associated(current))
       previous => current
       current  => current%next
       deallocate(previous)
    end do
    head => null()
    return
  end procedure intervalMultiListDestroy

  module procedure intervalMultiFailureReport
    !!{RST
    Report the state of a multi-integrand integration which has failed through loss of precision in the extent of an interval.
    !!}
    use :: Display           , only : displayIndent  , displayMessage     , displayUnindent       , displayVerbosity,                              displayVerbositySet, verbosityLevelStandard
    use :: ISO_Varying_String, only : assignment(=)  , operator(//)       , varying_string
    implicit none
    type            (intervalMulti ), pointer :: interval_
    type            (varying_string) :: message
    character       (len=32        ) :: label
    integer :: i

    if (displayVerbosity() < verbosityLevelStandard .and. .not.statusPresent) call displayVerbositySet(verbosityLevelStandard)
    call displayIndent('integration failure:')
    call displayIndent('current intervals:')
    message="a/b/(b-a)       ="
    write (label,'(e32.12)') a
    message=message//" "  //label
    write (label,'(e32.12)') b
    message=message//"/"  //label
    write (label,'(e32.12)') b-a
    message=message//"/(" //label//")"
    call displayMessage(message)
    message="a/b : f(a)/f(b) ="
    write (label,'(e32.12)') current%a
    message=message//" "  //label
    write (label,'(e32.12)') current%b
    message=message//"/"  //label
    call displayMessage(message)
    do i=1,size(current%fa)
       message="  (i="
       write (label,'(i12)') i
       message=message//adjustl(trim(label))//")"
       write (label,'(e32.12)') current%fa(i)
       message=message//" : "//label
       write (label,'(e32.12)') current%fb(i)
       message=message//"/"  //label
       write (label,'(l1)') converged(i)
       message=message//" (converged="//adjustl(trim(label))//")"
       call displayMessage(message)
    end do
    interval_ => head
    do while (associated(interval_))
       message="a/b : f(a)/f(b) ="
       write (label,'(e32.12)') interval_%a
       message=message//" "  //label
       write (label,'(e32.12)') interval_%b
       message=message//"/"  //label
       call displayMessage(message)
       do i=1,size(interval_%fa)
          message="  (i="
          write (label,'(i12)') i
          message=message//adjustl(trim(label))//")"
          write (label,'(e32.12)') interval_%fa(i)
          message=message//" : "//label
          write (label,'(e32.12)') interval_%fb(i)
          message=message//"/"  //label
          write (label,'(l1)') converged(i)
          message=message//" (converged="//adjustl(trim(label))//")"
          call displayMessage(message)
       end do
       interval_ => interval_%next
    end do
    call displayUnindent('done')
    call displayIndent  ('current integrals:')
    do i=1,size(integral)
       message="integral : error  (i="
       write (label,'(i12)') i
       message=message//adjustl(trim(label))//")"
       write (label,'(e32.12)') integral (i)
       message=message//" : "//label
       write (label,'(e32.12)') error    (i)
       message=message//"/"  //label
       write (label,'(l1)'    ) converged(i)
       message=message//" (converged="//adjustl(trim(label))//")"
       call displayMessage(message)
    end do
    call displayUnindent('done')
    call displayUnindent('done')
    return
  end procedure intervalMultiFailureReport

end submodule Numerical_Integration2_Intervals
