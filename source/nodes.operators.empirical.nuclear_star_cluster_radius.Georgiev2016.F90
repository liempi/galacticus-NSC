!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  Implements a node operator class that implements an an empirical power law relationship between disk stellar mass and
  stellar radius.
  !!}

  !![
  <nodeOperator name="nodeOperatorNuclearStarClusterRadiusGeorgiev2016">
   <description>
    A node operator that implements an an empirical power law relationship between nuclear star cluster radius and host galaxy stellar mass.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorNuclearStarClusterRadiusGeorgiev2016
     !!{
     Implements a power law prescription for the stellar mass--stellar radius relation of disks. Specifically:
     \begin{equation}
       r_\mathrm{s} = \radiusPivot \left( \frac{M_\star+M_\mathrm{gas}}{\massPivot} \right)^\radiusPivot, 
     \end{equation}
     where $r_\mathrm{s}$ is the nuclear star cluster scale radius, $M_\star$ is the stellar mass of the nuclear star cluster, 
     $M_\mathrm{gas}$ is the gaseous mass of the nuclear star cluster, and $\radiusPivot$, and $\massPivot$,
     are free parameters.
     !!}
     private
   contains
     !![
     <methods>
       <method method="update" description="Update the nuclear star cluster radius to be consistent with its stellar mass."/>
     </methods>
     !!]
     procedure :: update                              => nuclearStarClusterRadiusGeorgiev2016Update
     procedure :: nodeInitialize                      => nuclearStarClusterRadiusGeorgiev2016NodeInitialize
     procedure :: differentialEvolutionSolveAnalytics => nuclearStarClusterRadiusGeorgiev2016SolveAnalytics
     procedure :: nodesMerge                          => nuclearStarClusterRadiusGeorgiev2016NodesMerge
  end type nodeOperatorNuclearStarClusterRadiusGeorgiev2016
  
  interface nodeOperatorNuclearStarClusterRadiusGeorgiev2016
     !!{
     Constructors for the \refClass{nodeOperatorNuclearStarClusterRadiusGeorgiev2016} node operator class.
     !!}
     module procedure nuclearStarClusterRadiusGeorgiev2016ConstructorParameters
  end interface nodeOperatorNuclearStarClusterRadiusGeorgiev2016
  
contains

  function nuclearStarClusterRadiusGeorgiev2016ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorNuclearStarClusterRadiusGeorgiev2016} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type         (nodeOperatorNuclearStarClusterRadiusGeorgiev2016) :: self
    type         (inputParameters                                 ), intent(inout) :: parameters

    self=nodeOperatorNuclearStarClusterRadiusGeorgiev2016()
    !![
      <inputParametersValidate source="parameters"/>
    !!]
  end function nuclearStarClusterRadiusGeorgiev2016ConstructorParameters

  subroutine nuclearStarClusterRadiusGeorgiev2016Update(self,node)
    !!{
    Update radius of the nuclear star cluster.
    !!} 
    use :: Galacticus_Nodes, only : nodeComponentNSC, nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSCStandard
    implicit none
    class           (nodeOperatorNuclearStarClusterRadiusGeorgiev2016), intent(inout) :: self
    type            (treeNode                                        ), intent(inout) :: node
    class           (nodeComponentNSC                                ), pointer       :: nuclearStarCluster
    class           (nodeComponentDisk                               ), pointer       :: disk
    class           (nodeComponentSpheroid                           ), pointer       :: spheroid
    double precision                                                                  :: radius               , stellarMass, &
      &                                                                                  stellarMassMorphology
    double precision                                                  , parameter     :: massNormalizationLate    =5.61d9   ! M
    double precision                                                  , parameter     :: massNormalizationEarly   =2.09d9   ! M
    double precision                                                  , parameter     :: radiusNormalizationLate  =+3.44d0  ! Adimensional
    double precision                                                  , parameter     :: radiusNormalizationEarly =+6.11d0  ! Adimensional
    double precision                                                  , parameter     :: slopeLate                =+0.356d0 ! Adimensional
    double precision                                                  , parameter     :: slopeEarly               =+0.326d0 ! Adimensional
    double precision                                                  , parameter     :: interceptCoefficientLate =-0.012d0 ! Adimensional
    double precision                                                  , parameter     :: interceptCoefficientEarly=-0.011d0 ! Adimensional
    
    disk              => node%disk    ()
    spheroid          => node%spheroid()
    nuclearStarCluster=> node%NSC()

    select type(nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nothing to do here.
      return
    class default
      stellarMass= disk    %massStellar() &
        &         +spheroid%massStellar()
    
      if (stellarMass > 0.0d0) then
        stellarMassMorphology= spheroid%massStellar()/stellarMass
        if (stellarMassMorphology>0.2d0) then
          radius=+radiusNormalizationEarly                           &
            &    *1.0d1**(                                           &
            &              slopeEarly                                &
            &             *log10(stellarMass/massNormalizationEarly) &
            &             +interceptCoefficientEarly                 &
            &            )                                           &
            &    *1.0e-6
        else
          radius=+radiusNormalizationLate                           &
            &    *1.0d1**(                                          &
            &              slopeLate                                &
            &             *log10(stellarMass/massNormalizationLate) &
            &             +interceptCoefficientLate                 &
            &            )                                          &
            &    *1.0e-6
        end if 
      else 
        radius=0.0d0
      end if
    call nuclearStarCluster%radiusSet(radius)
    end select
    return  
  end subroutine nuclearStarClusterRadiusGeorgiev2016Update

  subroutine nuclearStarClusterRadiusGeorgiev2016NodeInitialize(self,node)
    !!{
    Initialize the nuclear star cluster.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class(nodeOperatorNuclearStarClusterRadiusGeorgiev2016), intent(inout), target  :: self
    type (treeNode                                        ), intent(inout), target  :: node
    class(nodeComponentNSC                                ),                pointer :: nuclearStarCluster

    nuclearStarCluster => node%NSC()
    call self%update(node)
    return
  end subroutine nuclearStarClusterRadiusGeorgiev2016NodeInitialize
   
  subroutine nuclearStarClusterRadiusGeorgiev2016SolveAnalytics(self,node,time)
    !!{
    Set the radius of the nuclear star cluster.
    !!}
    implicit none
    class           (nodeOperatorNuclearStarClusterRadiusGeorgiev2016), intent(inout) :: self
    type            (treeNode                                        ), intent(inout) :: node
    double precision                                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    call self%update(node)
    return
  end subroutine nuclearStarClusterRadiusGeorgiev2016SolveAnalytics

  subroutine nuclearStarClusterRadiusGeorgiev2016NodesMerge(self,node)
    !!{
    Update the radius of the nuclear star cluster after a merger.
    !!}
    implicit none
    class(nodeOperatorNuclearStarClusterRadiusGeorgiev2016), intent(inout) :: self
    type (treeNode                                        ), intent(inout) :: node

    call self%update(node)
    return
  end subroutine nuclearStarClusterRadiusGeorgiev2016NodesMerge
