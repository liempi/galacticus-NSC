!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a program to test multidimensional minimizers.

program Test_Multidimensional_Minimizer
  !% Tests of multidimensional minimizers.
  use            :: Galacticus_Display                       , only : Galacticus_Verbosity_Level_Set, verbosityStandard
  use            :: Multidimensional_Minimizer               , only : multiDMinimizer
  use            :: Test_Multidimensional_Minimizer_Functions, only : minimizerFunction_            , minimizeFunctionDerivative_, minimizeFunctionBoth_
  use            :: Unit_Tests                               , only : Assert                        , Unit_Tests_Begin_Group     , Unit_Tests_End_Group, Unit_Tests_Finish, &
       &                                                              compareLessThan
  use, intrinsic :: ISO_C_Binding                            , only : c_size_t
  implicit none
  type            (multiDMinimizer), allocatable  :: minimizer_
  double precision                 , dimension(2) :: x
  integer                                         :: iteration
  logical                                         :: converged

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Multidimensional minimizer")
  !! Build the minimizer.
  allocate(minimizer_)
  minimizer_=multiDMinimizer(2_c_size_t,minimizerFunction_,minimizeFunctionDerivative_,minimizeFunctionBoth_)
  call minimizer_%set(x=[5.0d0,7.0d0],stepSize=0.01d0,tolerance=1.0d-4)
  ! Perform the minimization.
  iteration=0
  converged=.false.
  do while (.not.converged.and.iteration < 100)
     iteration=iteration+1
     call minimizer_%iterate()
     converged=minimizer_%testGradient(toleranceAbsolute=1.0d-3)
  end do
  x=minimizer_%x()
  call Assert('converged',iteration,100,compareLessThan)
  call Assert('minimum',x,[1.0d0,2.0d0],relTol=1.0d-6)
  deallocate(minimizer_)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Multidimensional_Minimizer