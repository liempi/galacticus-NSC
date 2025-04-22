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
  Implements a merger nuclear star cluster movement class which uses a simple calculation.
  !!}

  !![
  <nuclearStarClusterMovements name="nuclearStarClusterMovementsSimple">
   <description>
    A merger mass movements class which implements mass movements according to:
   </description>
  </nuclearStarClusterMovements>
  !!]
  type, extends(nuclearStarClusterMovementsClass) :: nuclearStarClusterMovementsSimple
     !!{
     A merger mass movements class which uses a simple calculation.
     !!}
     private
   contains
     procedure :: isDestroyed => simpleIsDestroyed
  end type nuclearStarClusterMovementsSimple

  interface nuclearStarClusterMovementsSimple
     !!{
     Constructors for the {\normalfont \ttfamily simple} merger mass movements class.
     !!}
     module procedure simpleConstructorParameters
  end interface nuclearStarClusterMovementsSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily simple} merger mass movements class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nuclearStarClusterMovementsSimple)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters

    self=nuclearStarClusterMovementsSimple()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  logical function simpleIsDestroyed(self,node,isDestroyed)
    !!{
    Determine if a nuclear star cluster survives a merger event using a simple criteria.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentNSC, nodeComponentSpheroid, nodeComponentDisk
    use :: Galactic_Structure_Options, only : componentTypeAll, massTypeStellar
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (nuclearStarClusterMovementsSimple), intent(inout)         :: self
    type            (treeNode                         ), intent(inout), target :: node
    logical                                            , intent(  out)         :: isDestroyed
    class           (nodeComponentNSC                 ), pointer               :: nuclearStarCluster
    class           (nodeComponentSpheroid            ), pointer               :: spheroid                  , spheroidHost           
    class           (nodeComponentDisk                ), pointer               :: disk                      , diskHost               
    class           (massDistributionClass            ), pointer               :: massDistribution_
    type            (treeNode                         ), pointer               :: nodeHost
    double precision                                                           :: stellarMassNuclearStarCluster, radiusNuclearStarCluster, &
          &                                                                       halfMassRadiusHost           , halfMassRadiusSatellite , &
          &                                                                       distance                     , stellarMassHost         , &
          &                                                                       stellarMassHost              , tidalRadius

    nodeHost           => node    %mergesWith()
    diskHost           => nodeHost%disk      ()
    spheroidHost       => nodeHost%spheroid  ()
    nuclearStarCluster => node    %NSC       ()
    disk               => node    %disk      ()
    spheroid           => node    %spheroid  ()

    halfMassRadiusHost       = max(                               &  
          &                        diskHost    %halfMassRadius(), &
          &                        spheroidHost%halfMassRadius()  &
          &                       )

    halfMassRadiusSatellite  = max(                           &  
          &                        disk    %halfMassRadius(), &
          &                        spheroid%halfMassRadius()  &
          &                       )

    radiusNuclearStarCluster = nuclearStarCluster%radius()

    ! Assume that distance between galaxies is given by the sum of their half mass radii.
    distance                = halfMassRadiusHost+halfMassRadiusSatellite


    massDistribution_             => nodeHost        %massDistribution(                                &
          &                                                            componentType=componentTypeAll, &
          &                                                            massType     =massTypeStellar   &
          &                                                           )

    stellarMassHost               = massDistribution_ %massTotal  ()
    stellarMassNuclearStarCluster = nuclearStarCluster%massStellar()

    ! Compute the tidal radius as defined in I. King. The structure of star clusters. I. an empirical density law. AJ, 67:471, Oct.1962. doi: 10.1086/108756.
    tidalRadius = + distance                       &
          &       *(                               &
          &         stellarMassNuclearStarCluster/ &
          &         stellarMassHost              / &
          &         2.0d0                          &
          &        )                               &
          &        **(1.0d0/3.0d0)

    if (tidalRadius >= radiusNuclearStarCluster) then 
        isDestroyed=.false.
    else
        isDestroyed=.true.
    end if
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function simpleIsDestroyed
