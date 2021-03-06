/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
inline const CloudType& Foam::ParticleTorqueList<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
inline CloudType& Foam::ParticleTorqueList<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
inline const Foam::fvMesh& Foam::ParticleTorqueList<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::dictionary& Foam::ParticleTorqueList<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
inline void Foam::ParticleTorqueList<CloudType>::setCalcCoupled(bool flag)
{
    calcCoupled_ = flag;
}


template<class CloudType>
inline void Foam::ParticleTorqueList<CloudType>::setCalcNonCoupled(bool flag)
{
    calcNonCoupled_ = flag;
}


// ************************************************************************* //
