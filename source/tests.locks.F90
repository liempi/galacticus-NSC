!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a program to test OpenMP lock functions.

program Test_Locks
  !% Tests of OpenMP locking functions.
  use, intrinsic :: ISO_C_Binding
  use            :: Unit_Tests
  use            :: Galacticus_Display
  use            :: Locks
  use            :: Array_Utilities
  use            :: String_Handling
  use            :: ISO_Varying_String
  !$ use         :: OMP_Lib
  implicit none
  integer         (c_size_t          ), parameter               :: elementCount   =100_c_size_t
  integer         (c_size_t          ), dimension(elementCount) :: orderedCount
  type            (ompIncrementalLock)                          :: incrementalLock
  integer         (c_size_t          )                          :: i                           , unorderedCounter       , &
       &                                                           orderedCounter              , unorderedCounterPrivate
  double precision                                              :: uniformRandom
  integer                                                       :: sleepTime
  type            (varying_string    )                          :: message

  call Unit_Tests_Begin_Group("OpenMP")
  ! Test incremental locks. We generate a counter which is not guaranteed to be ordered in terms of OpenMP threads. Then we use an
  ! incremental lock to force the access back to being ordered by counter value and thereby construct and ordered array.
  call random_seed()
  incrementalLock =ompIncrementalLock()
  orderedCounter  =0_c_size_t
  unorderedCounter=0_c_size_t
  !$omp parallel do private(unorderedCounterPrivate,message,uniformRandom,sleepTime)
  do i=1_c_size_t,elementCount
     !$omp critical(increment)
     unorderedCounter       =unorderedCounter+1_c_size_t
     unorderedCounterPrivate=unorderedCounter
     !$omp end critical(increment)
     ! Sleep for a random duration to force desynchronization of threads.
     call random_number(uniformRandom)
     sleepTime=int(uniformRandom*3.0d0)
     message=var_str("unordered counter retrieved with value ")//unorderedCounterPrivate//var_str("; pausing for ")//sleepTime//var_str(" seconds")
     call Galacticus_Display_Message(message)
     call sleep        (sleepTime)
     call incrementalLock%set  (unorderedCounterPrivate)
     orderedCounter                =orderedCounter         +1_c_size_t
     message=var_str("ordered counter retrieved with value ")//orderedCounter
     call Galacticus_Display_Message(message)
     orderedCount  (orderedCounter)=unorderedCounterPrivate
     call incrementalLock%unset( )     
  end do
  !$omp end parallel do
  call Assert('OpenMP incremental lock',Array_Is_Monotonic(orderedCount,direction=directionIncreasing,allowEqual=.false.),.true.)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Locks