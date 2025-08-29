!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a black hole seed based ...
  !!}
 
  !![
  <blackHoleSeeds name="blackHoleSeedsVeryMassiveStars">
    <description>
      A model of black hole seeds in which seeds are formed due to the...
    </description>
  </blackHoleSeeds>
  !!]

  type, extends(blackHoleSeedsClass) :: blackHoleSeedsVeryMassiveStars
     !!{
     A black hole seeds class in which seeds are formed as a result of very massive stars...
     !!}
     private
     double precision                                                  :: massFraction                       , nuclearStarClusterMaximumAge
     integer                                                           :: ageNuclearStarClustersID           , gasMassNuclearStarClustersID             , &
       &                                                                  stellarMassNuclearStarClustersID   , nuclearStarClusterFormationTimeID        , &
       &                                                                  redshiftBlackHoleSeedFormationVMSID, radiusNuclearStarClustersID              , &
       &                                                                  blackHoleSeedMassID                , stellarMassFormedNSCID                   , &
       &                                                                  timeStellarMassFormedNSCID         , coreCollapseTimescaleNuclearStarClusterID
   contains   
     procedure :: mass             => veryMassiveStarsMass
     procedure :: spin             => veryMassiveStarsSpin
     procedure :: formationChannel => veryMassiveStarsFormationChannel
  end type blackHoleSeedsVeryMassiveStars
  
  interface blackHoleSeedsVeryMassiveStars
     !!{
     Constructors for the {\normalfont \ttfamily veryMassiveStars} black hole seeds class.
     !!}
     module procedure veryMassiveStarsConstructorParameters
     module procedure veryMassiveStarsConstructorInternal
  end interface blackHoleSeedsVeryMassiveStars

contains

  function veryMassiveStarsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily veryMassiveStars} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleSeedsVeryMassiveStars)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: massFraction, nuclearStarClusterMaximumAge 

    !![
    <inputParameter>
      <name>massFraction</name>
      <defaultValue>0.0824d0</defaultValue>
      <description>Specifies the efficiency of very massive stars which form a black hole seed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>nuclearStarClusterMaximumAge</name>
      <defaultValue>5.0d-3</defaultValue>
      <description>Specifies the maximum age (Gyr) of the nuclear star cluster to apply this formation channel.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=blackHoleSeedsVeryMassiveStars(massFraction,nuclearStarClusterMaximumAge)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function veryMassiveStarsConstructorParameters
  
  function veryMassiveStarsConstructorInternal(massFraction,nuclearStarClusterMaximumAge) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily veryMassiveStars} node operator class.
    !!}
    implicit none
    type            (blackHoleSeedsVeryMassiveStars)                :: self
    double precision                                , intent(in   ) :: massFraction
    double precision                                , intent(in   ) :: nuclearStarClusterMaximumAge
    !![
    <constructorAssign variables="massFraction, nuclearStarClusterMaximumAge"/>
    <addMetaProperty component="NSC" name="agesStellarMassFormed"                   id="self%stellarMassFormedNSCID"                    isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="NSC" name="agesTimeStellarMassFormed"               id="self%timeStellarMassFormedNSCID"                isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="NSC" name="blackHoleSeedMassFormed"                 id="self%blackHoleSeedMassID"                       isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="NSC" name="ageNuclearStarClusters"                  id="self%ageNuclearStarClustersID"                  isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="NSC" name="radiusNuclearStarClusters"               id="self%radiusNuclearStarClustersID"               isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="NSC" name="gasMassNuclearStarClusters"              id="self%gasMassNuclearStarClustersID"              isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="NSC" name="redshiftBlackHoleSeedFormation"          id="self%redshiftBlackHoleSeedFormationVMSID"       isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="NSC" name="stellarMassNuclearStarClusters"          id="self%stellarMassNuclearStarClustersID"          isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="NSC" name="coreCollapseTimescaleNuclearStarCluster" id="self%coreCollapseTimescaleNuclearStarClusterID" isEvolvable="no"  isCreator="yes"/>   
    <addMetaProperty component="NSC" name="nuclearStarClusterFormationTime"         id="self%nuclearStarClusterFormationTimeID"         isEvolvable="no"  isCreator="no" />
     !!]
    return
  end function veryMassiveStarsConstructorInternal

  double precision function veryMassiveStarsMass(self,node) result(mass)
      !!{
        Compute the nuclear star cluster collapse condition.
      !!}
    use :: Galacticus_Nodes                , only : nodeComponentNSC                     , nodeComponentBasic                     , nodeComponentNSCStandard, treeNode  
    use :: Abundances_Structure            , only : operator(*)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Galactic_Structure_Options      , only : componentTypenuclearStarCluster      , massTypeStellar
    use :: Numerical_Constants_Prefixes    , only : mega                                 , kilo
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal       , radiusSolar                            , megaParsec              , gigaYear, &
        &                                           parsec
    use :: Nuclear_Star_Clusters_Utilities , only : nuclearStarClusterKroupaNumberOfStars, nuclearStarClusterCoreCollapseTimescale
    implicit none
    class           (blackHoleSeedsVeryMassiveStars), intent(inout)          :: self
    type            (treeNode                      ), intent(inout)          :: node
    class           (nodeComponentNSC              )               , pointer :: nuclearStarCluster
    class           (nodeComponentBasic            )               , pointer :: basic
    double precision                                                         :: radiusNuclearStarCluster         , velocityNuclearStarCluster     , &
        &                                                                       massStellarNuclearStarCluster    , ageNuclearStarCluster          , &
        &                                                                       massTimeStellarNuclearStarCluster, coreCollapseTimescale          , &
        &                                                                       time                             , formationTimeNuclearStarCluster
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    mass=0.0d0

    ! Detect the type of the nuclear star cluster component.
    select type (nuclearStarCluster)
      class default
          ! Generic type, do nothing.
          return
      class is (nodeComponentNSCStandard)
          ! Standard class, get the properties of the nuclear star cluster component.
          ! Unphysical nuclear star cluster, do nothing.
          if   (                                         &
             &   nuclearStarCluster%massStellar()<=0.0d0 &
             &  .or.                                     &
             &   nuclearStarCluster%radius     ()<=0.0d0 &
             &  .or.                                     &
             &   nuclearStarCluster%massGas    ()<=0.0d0 &
             & ) return
          massStellarNuclearStarCluster     =  nuclearStarCluster%floatRank0MetaPropertyGet(self%    stellarMassFormedNSCID       )
          massTimeStellarNuclearStarCluster =  nuclearStarCluster%floatRank0MetaPropertyGet(self%timeStellarMassFormedNSCID       )
          formationTimeNuclearStarCluster   =  nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationTimeID)
          basic                             => node              %basic                    (                                      )
          time                              =  basic             %time                     (                                      )
          if (formationTimeNuclearStarCluster>0.0d0) then
             ageNuclearStarCluster=+time                             &
                  &                -formationTimeNuclearStarCluster
          else 
             ageNuclearStarCluster=+0.0d0
          end if
          ! Do nothing if the nuclear star cluster has an unphysical age or already collapsed and formed a black hole seed.
          if   (                                 &
             &  ageNuclearStarCluster<=0.0d0     &
             &  .or.                             &
             &  nuclearStarCluster%isCollapsed() &
             & ) return

          radiusNuclearStarCluster=nuclearStarCluster%radius()
          ! Get the age of the nuclear star cluster.
          velocityNuclearStarCluster        =  sqrt(                                  &
               &                                    +gravitationalConstant_internal   &
               &                                    *nuclearStarCluster%massStellar() & 
               &                                    /radiusNuclearStarCluster         &
               &                                   ) 
          
          coreCollapseTimescale = nuclearStarClusterCoreCollapseTimescale(                                  &
              &                                                                   radiusNuclearStarCluster, &
              &                                                           nuclearStarCluster%massStellar(), &
              &                                                           nuclearStarCluster%massGas    (), &
              &                                                                                     .true.  &
              &                                                           )

          if (                                                            &
              &  coreCollapseTimescale>0.0d0                              &
              &  .and.                                                    &
              &  ageNuclearStarCluster<=self%nuclearStarClusterMaximumAge &
              &  .and.                                                    &
              &  coreCollapseTimescale<=ageNuclearStarCluster             &
              & ) then
            call nuclearStarCluster%floatRank0MetaPropertySet(self%ageNuclearStarClustersID                 ,                                       ageNuclearStarCluster                                                          )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%gasMassNuclearStarClustersID             , nuclearStarCluster                    %massGas                       (                                              ))
            call nuclearStarCluster%floatRank0MetaPropertySet(self%stellarMassNuclearStarClustersID         , nuclearStarCluster                    %massStellar                   (                                              ))
            !call nuclearStarCluster%floatRank0MetaPropertySet(self%redshiftBlackHoleSeedFormationVMSID      , self              %cosmologyFunctions_%redshiftFromExpansionFactor   (self%cosmologyFunctions_%expansionFactor(time)))
            call nuclearStarCluster%floatRank0MetaPropertySet(self%radiusNuclearStarClustersID              ,                                       radiusNuclearStarCluster                                                      )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%coreCollapseTimescaleNuclearStarClusterID, coreCollapseTimescale)
            ! Here, self%massFraction is computed in the following way: 
            ! For a initial mass function  φ(M), asumming the masses of the star in the range min ≤ M [M☉] ≤ max,
            ! the total mass is given by ∫ ₘᵢₙ ᵐᵃˣ  φ(M) M dM. It is possible to estimate the fraction of the 
            ! total mass of the stars more massive than 100 M☉ as
            ! f = ∫₁₀₀  ᵐᵃˣ  φ(M) M dM / ∫ ₘᵢₙ ᵐᵃˣ  φ(M) M dM
            ! Here we take into account all the massive stars. We should include an efficiency parameter.
            mass   =+self%massFraction  &
                 &  *nuclearStarCluster%massStellar   ()
            call nuclearStarCluster%           isCollapsedSet(                         .true.)
            call nuclearStarCluster%floatRank0MetaPropertySet(self%blackHoleSeedMassID,mass  )
            
            ! Adjust stellar mass of the nuclear star cluster
            call nuclearStarCluster%           massStellarSet(                                         &
                 &                                            +(                                       &
                 &                                              +1.0d0                                 &
                 &                                              -self%massFraction                     &
                 &                                             )                                       &
                 &                                            *nuclearStarCluster%      massStellar()  &
                 &                                           )
            ! Adjust stellar abundances of the nuclear star cluster 
            call nuclearStarCluster%     abundancesStellarSet(                                         &
                 &                                            +(                                       &
                 &                                              +1.0d0                                 &
                 &                                              -self%massFraction                     &
                 &                                             )                                       &
                 &                                            *nuclearStarCluster%abundancesStellar()  &
                 &                                           )            
          else
            mass   =+0.0d0
          end if 
    end select
    return
  end function veryMassiveStarsMass

  double precision function veryMassiveStarsSpin(self,node) result(spin)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsVeryMassiveStars), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Assume zero spin.
    spin=0.0d0
    return
  end function veryMassiveStarsSpin

  function veryMassiveStarsFormationChannel (self,node) result(channel)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    type (enumerationBlackHoleFormationChannelType)                :: channel
    class(blackHoleSeedsVeryMassiveStars          ), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    channel=blackHoleFormationChannelVeryMassiveStars
    return
  end function veryMassiveStarsFormationChannel
