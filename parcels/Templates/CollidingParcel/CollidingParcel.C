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

#include "CollidingParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

/*template<class ParcelType>
template<class TrackData>
void Foam::CollidingParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::setCellValues(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::CollidingParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
	ParcelType::cellValueSourceCorrection(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::CollidingParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
	ParcelType::calc(td, dt, cellI);
}*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const CollidingParcel<ParcelType>& p
)
:
    ParcelType(p),
    f_(p.f_),
    //angularMomentum_(p.angularMomentum_), //changed
    torque_(p.torque_),
    collisionRecords_(p.collisionRecords_)
{}


template<class ParcelType>
Foam::CollidingParcel<ParcelType>::CollidingParcel
(
    const CollidingParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    f_(p.f_),
    //angularMomentum_(p.angularMomentum_), //changed
    torque_(p.torque_),
    collisionRecords_(p.collisionRecords_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::CollidingParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    switch (td.part())
    {
        case TrackData::tpVelocityHalfStep:
        {
            // First and last leapfrog velocity adjust part, required
            // before and after tracking and force calculation

            p.U() += 0.5*trackTime*p.f()/p.mass();

            //p.angularMomentum() += 0.5*trackTime*p.torque(); //changed
			p.omega() += 0.5*trackTime*p.torque()/p.momentOfInertia();

            td.keepParticle = true;

            break;
        }

        case TrackData::tpLinearTrack:
        {
            ParcelType::move(td, trackTime);

            break;
        }

        case TrackData::tpRotationalTrack: 
        {
			//It is included in "TrackData::tpLinearTrack" via using kinematicRotParcel as superclass
            //notImplemented("TrackData::tpRotationalTrack"); //changed

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "CollidingParcel<ParcelType>::move(TrackData&, const scalar)"
            )   << td.part() << " is an invalid part of the tracking method."
                << abort(FatalError);
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
void Foam::CollidingParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    f_ = transform(T, f_);

    //angularMomentum_ = transform(T, angularMomentum_); //changed

    torque_ = transform(T, torque_);
}


template<class ParcelType>
void Foam::CollidingParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "CollidingParcelIO.C"

// ************************************************************************* //
