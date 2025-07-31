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

  !+    Contributions to this file made by: Matías Liempi

  !!{
  Implements a property extractor class the properties of nuclear star cluster when a black hole seed is formed using the model of \cite{vergara_global_2023}.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorBlackHoleSeedingVMS">
   <description>
    A property extractor class for the properties of the nuclear star cluster at the moment of the black hole formation.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorBlackHoleSeedingVMS
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     integer  :: radiusNuclearStarClustersID          , blackHoleSeedMassID              , &
         &       velocityNuclearStarClustersID        , ageNuclearStarClustersID         , &
         &       gasMassNuclearStarClustersID         , stellarMassNuclearStarClustersID         , &
         &       redshiftBlackHoleSeedFormationVMSID  , coreCollapseTimescaleNuclearStarClusterID
   contains
     procedure :: elementCount       => blackHoleSeedingVMSElementCount
     procedure :: extract            => blackHoleSeedingVMSExtract
     procedure :: names              => blackHoleSeedingVMSNames
     procedure :: descriptions       => blackHoleSeedingVMSDescriptions
     procedure :: unitsInSI          => blackHoleSeedingVMSUnitsInSI
  end type nodePropertyExtractorBlackHoleSeedingVMS

  interface nodePropertyExtractorBlackHoleSeedingVMS
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleSeedingVMS} output analysis class.
     !!}
     module procedure blackHoleSeedingVMSConstructorParameters
     module procedure blackHoleSeedingVMSConstructorInternal
  end interface nodePropertyExtractorBlackHoleSeedingVMS

contains

  function blackHoleSeedingVMSConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily blackHoleSeedingVMS} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorBlackHoleSeedingVMS)                :: self
    type (inputParameters                         ), intent(inout) :: parameters

    self=nodePropertyExtractorBlackHoleSeedingVMS()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleSeedingVMSConstructorParameters

  function blackHoleSeedingVMSConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleSeedingVMS} property extractor class.
    !!}
    implicit none
    type          (nodePropertyExtractorBlackHoleSeedingVMS) :: self
    !![
    <addMetaProperty component="NSC" name="blackHoleSeedMassFormed"                 id="self%blackHoleSeedMassID"                       isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="ageNuclearStarClusters"                  id="self%ageNuclearStarClustersID"                  isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="radiusNuclearStarClusters"               id="self%radiusNuclearStarClustersID"               isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="gasMassNuclearStarClusters"              id="self%gasMassNuclearStarClustersID"              isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="stellarMassNuclearStarClusters"          id="self%stellarMassNuclearStarClustersID"          isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="redshiftBlackHoleSeedFormation"          id="self%redshiftBlackHoleSeedFormationVMSID"       isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="NSC" name="coreCollapseTimescaleNuclearStarCluster" id="self%coreCollapseTimescaleNuclearStarClusterID" isEvolvable="no" isCreator="no"/>   
    !!]
    return
  end function blackHoleSeedingVMSConstructorInternal

  integer function blackHoleSeedingVMSElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily blackHoleSeedingVMS} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingVMS), intent(inout) :: self
    double precision                                          , intent(in   ) :: time
    !$GLC attributes unused :: time

    blackHoleSeedingVMSElementCount=7
    return
  end function blackHoleSeedingVMSElementCount

  function blackHoleSeedingVMSExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily blackHoleSeedingVMS} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    double precision                                          , dimension(:) , allocatable :: blackHoleSeedingVMSExtract
    class           (nodePropertyExtractorBlackHoleSeedingVMS), intent(inout), target      :: self
    type            (treeNode                                ), intent(inout), target      :: node
    double precision                                          , intent(in   )              :: time
    type            (multiCounter                            ), intent(inout), optional    :: instance
    class           (nodeComponentNSC                        )               , pointer     :: nuclearStarCluster

    !$GLC attributes unused :: time, instance

    allocate(blackHoleSeedingVMSExtract(7))
    nuclearStarCluster => node%NSC()
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nuclear star cluster does not yet exist.
      blackHoleSeedingVMSExtract=[        & 
        &                          0.0d0, &
        &                          0.0d0, &
        &                          0.0d0, &
        &                          0.0d0, &
        &                          0.0d0, &
        &                          0.0d0, &
        &                          0.0d0  &
        &                        ]
    class default
      blackHoleSeedingVMSExtract=[                                                                                               &
       &                           nuclearStarCluster%floatRank0MetaPropertyGet(self%redshiftBlackHoleSeedFormationVMSID      ), &
       &                           nuclearStarCluster%floatRank0MetaPropertyGet(self%blackHoleSeedMassID                      ), &
       &                           nuclearStarCluster%floatRank0MetaPropertyGet(self%ageNuclearStarClustersID                 ), &
       &                           nuclearStarCluster%floatRank0MetaPropertyGet(self%radiusNuclearStarClustersID              ), &
       &                           nuclearStarCluster%floatRank0MetaPropertyGet(self%gasMassNuclearStarClustersID             ), &
       &                           nuclearStarCluster%floatRank0MetaPropertyGet(self%stellarMassNuclearStarClustersID         ), &
       &                           nuclearStarCluster%floatRank0MetaPropertyGet(self%coreCollapseTimescaleNuclearStarClusterID)  &
       &                          ]
    end select
    return
  end function blackHoleSeedingVMSExtract

  subroutine blackHoleSeedingVMSNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily blackHoleSeedingVMS} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingVMS), intent(inout)                             :: self
    double precision                                          , intent(in   )                             :: time
    type            (varying_string                          ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(7))
    names(1)=var_str('blackHoleFormationRedshift'             )
    names(2)=var_str('blackHoleSeedMass'                      )
    names(3)=var_str('nuclearStarClusterAge'                  )
    names(4)=var_str('nuclearStarClusterRadius'               )
    names(5)=var_str('nuclearStarClusterGasMass'              )
    names(6)=var_str('nuclearStarClusterStellarMass'          )
    names(7)=var_str('nuclearStarClusterCoreCollapseTimescale')
    return
  end subroutine blackHoleSeedingVMSNames

  subroutine blackHoleSeedingVMSDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily blackHoleSeedingVMS} property.
    !!}
    implicit none
    class           (nodePropertyExtractorBlackHoleSeedingVMS), intent(inout)                             :: self
    double precision                                          , intent(in   )                             :: time
    type            (varying_string                          ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(7))
    descriptions(1)=var_str('Redshift at the formation of the black hole seed .'                                     )
    descriptions(2)=var_str('Black hole seed mass formed due to VMS model [M⊙].'                                     )
    descriptions(3)=var_str('Mass-weighted age of the nuclear star cluster at which black hole seed is formed [Gyr].')
    descriptions(4)=var_str('Radius of the nuclear star cluster [Mpc].'                                              )
    descriptions(5)=var_str('Gaseous mass of the nuclear star cluster [M⊙].'                                         )
    descriptions(6)=var_str('Stellar mass of the nuclear star cluster [M⊙].'                                         )
    descriptions(7)=var_str('Core collapse timescale of the nuclear star cluster [Gyr].'                             )
    return
  end subroutine blackHoleSeedingVMSDescriptions

  function blackHoleSeedingVMSUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily blackHoleSeedingVMS} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec, gigayear
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                          , allocatable  , dimension(:) :: blackHoleSeedingVMSUnitsInSI
    class           (nodePropertyExtractorBlackHoleSeedingVMS), intent(inout)               :: self
    double precision                                          , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(blackHoleSeedingVMSUnitsInSI(7))
    blackHoleSeedingVMSUnitsInSI=[             &
     &                             1.0d0     , &
     &                             massSolar , &
     &                             gigayear  , &
     &                             megaParsec, &
     &                             massSolar , &
     &                             massSolar , &
     &                             gigayear    &
     &                           ]
    return
  end function blackHoleSeedingVMSUnitsInSI
