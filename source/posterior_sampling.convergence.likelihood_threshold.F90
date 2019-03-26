!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of a posterior sampling convergence class which declares convergence once all likelihoods are above a threshold.

  !# <posteriorSampleConvergence name="posteriorSampleConvergenceLikelihoodThreshold">
  !#  <description>A posterior sampling convergence class which declares convergence once all likelihoods are above a threshold.</description>
  !# </posteriorSampleConvergence>
  type, extends(posteriorSampleConvergenceClass) :: posteriorSampleConvergenceLikelihoodThreshold
     !% Implementation of a posterior sampling convergence class which declares convergence once all likelihoods are above a threshold.
     private
     integer          :: convergedAtStepCount
     logical          :: converged
     double precision :: likelihoodThreshold
   contains
     procedure :: isConverged     => likelihoodThresholdIsConverged
     procedure :: convergedAtStep => likelihoodThresholdConvergedAtStep
     procedure :: reset           => likelihoodThresholdReset
     procedure :: logReport       => likelihoodThresholdLogReport
     procedure :: stateIsOutlier  => likelihoodThresholdStateIsOutlier
  end type posteriorSampleConvergenceLikelihoodThreshold

  interface posteriorSampleConvergenceLikelihoodThreshold
     !% Constructors for the {\normalfont \ttfamily likelihoodThreshold} posterior sampling convergence class.
     module procedure likelihoodThresholdConstructorParameters
     module procedure likelihoodThresholdConstructorInternal
  end interface posteriorSampleConvergenceLikelihoodThreshold

contains

  function likelihoodThresholdConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily likelihoodThreshold} posterior sampling convergence class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (posteriorSampleConvergenceLikelihoodThreshold)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    double precision                                                               :: likelihoodThreshold

    !# <inputParameter>
    !#   <name>likelihoodThreshold</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The threshold log-likelihood above which convergence is declared.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=posteriorSampleConvergenceLikelihoodThreshold(likelihoodThreshold)
    !# <inputParametersValidate source="parameters"/>
    return
  end function likelihoodThresholdConstructorParameters

  function likelihoodThresholdConstructorInternal(likelihoodThreshold) result(self)
    !% Internal constructor for the {\normalfont \ttfamily likelihoodThreshold} posterior sampling convergence class.
    use Input_Parameters
    implicit none
    type            (posteriorSampleConvergenceLikelihoodThreshold)                :: self
    double precision                                               , intent(in   ) :: likelihoodThreshold
    !# <constructorAssign variables="likelihoodThreshold"/>

    self%converged           =.false.
    self%convergedAtStepCount=huge(0)
    return
  end function likelihoodThresholdConstructorInternal

  logical function likelihoodThresholdIsConverged(self,simulationState,logLikelihood)
    !% Returns true if the posterior sampling is converged (which it likelihoodThreshold is).
    use MPI_Utilities
    implicit none
    class           (posteriorSampleConvergenceLikelihoodThreshold), intent(inout)           :: self
    class           (posteriorSampleStateClass                    ), intent(inout), optional :: simulationState
    double precision                                               , intent(in   ), optional :: logLikelihood

    ! If no arguments were provided, return current convergence status without updating.
    if (.not.(present(simulationState).and.present(logLikelihood))) then
       likelihoodThresholdIsConverged=self%converged
       return
    end if
    ! Convergence requires all chains to have a likelihood above the threshold.
    if (.not.self%converged.and.mpiSelf%all(logLikelihood >= self%likelihoodThreshold)) then
       self%converged=.true.
       self%convergedAtStepCount=simulationState%count()
    end if
    likelihoodThresholdIsConverged=self%converged
    return
  end function likelihoodThresholdIsConverged

  integer function likelihoodThresholdConvergedAtStep(self)
    !% Return the step at which the simulation converged.
    implicit none
    class(posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self

    likelihoodThresholdConvergedAtStep=self%convergedAtStepCount
    return
  end function likelihoodThresholdConvergedAtStep

  subroutine likelihoodThresholdReset(self)
    !% Reset the convergence object.
    implicit none
    class(posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self

    self%converged           =.false.
    self%convergedAtStepCount=-1
    return
  end subroutine likelihoodThresholdReset

  subroutine likelihoodThresholdLogReport(self,fileUnit)
    !% Write a convergence report to the given {\normalfont \ttfamily fileUnit}.
    implicit none
    class  (posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self
    integer                                               , intent(in   ) :: fileUnit

    if (self%converged) then
       write (fileUnit,*) 'Convergence: converged'
    else
       write (fileUnit,*) 'Convergence: unconverged'
    end if
    return
  end subroutine likelihoodThresholdLogReport

  logical function likelihoodThresholdStateIsOutlier(self,stateIndex)
    !% Return true if the specified chain is deemed to be an outlier. In this case, chains are likelihoodThreshold outliers.
    implicit none
    class  (posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self
    integer                                               , intent(in   ) :: stateIndex
    !GCC$ attributes unused :: self, stateIndex

    likelihoodThresholdStateIsOutlier=.false.
    return
  end function likelihoodThresholdStateIsOutlier