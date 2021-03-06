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

#include "FrontTrackingParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::FrontTrackingParcel<ParcelType>::propertyList_ =
    Foam::FrontTrackingParcel<ParcelType>::propertyList();

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::FrontTrackingParcel<ParcelType>::FrontTrackingParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    active_(false),
    typeId_(0),
    currentIndex_(0),//Jafari added
    bubbleIndex_(0),//Jafari added
    nParticle_(0.0),
    d_(0.0),
    dTarget_(0.0),
    U_(vector::zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(vector::zero),
    rhoc_(0.0),
    Uc_(vector::zero),
    muc_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            currentIndex_ = readLabel(is);//Jafari added
            bubbleIndex_ = readLabel(is);//Jafari added
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            dTarget_ = readScalar(is);
            is >> U_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
            tTurb_ = readScalar(is);
            is >> UTurb_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&active_),
                sizeof(active_)
              + sizeof(typeId_)
              + sizeof(currentIndex_)
              + sizeof(bubbleIndex_)
              + sizeof(nParticle_)
              + sizeof(d_)
              + sizeof(dTarget_)
              + sizeof(U_)
              + sizeof(rho_)
              + sizeof(age_)
              + sizeof(tTurb_)
              + sizeof(UTurb_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "FrontTrackingParcel<ParcelType>::FrontTrackingParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::FrontTrackingParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<label> active(c.fieldIOobject("active", IOobject::MUST_READ));
    c.checkFieldIOobject(c, active);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    //Jafari added
    IOField<label> currentIndex(c.fieldIOobject("currentIndex", IOobject::MUST_READ));
    c.checkFieldIOobject(c, currentIndex);

    //Jafari added
    IOField<label> bubbleIndex(c.fieldIOobject("bubbleIndex", IOobject::MUST_READ));
    c.checkFieldIOobject(c, bubbleIndex);

    IOField<scalar>
        nParticle(c.fieldIOobject("nParticle", IOobject::MUST_READ));
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d(c.fieldIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::MUST_READ));
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age(c.fieldIOobject("age", IOobject::MUST_READ));
    c.checkFieldIOobject(c, age);

    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::MUST_READ));
    c.checkFieldIOobject(c, UTurb);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        FrontTrackingParcel<ParcelType>& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.currentIndex_ = currentIndex[i];//Jafari added
        p.bubbleIndex_ = bubbleIndex[i];//Jafari added
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::FrontTrackingParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<label> currentIndex(c.fieldIOobject("currentIndex", IOobject::NO_READ), np);//Jafari added
    IOField<label> bubbleIndex(c.fieldIOobject("bubbleIndex", IOobject::NO_READ), np);//Jafari added
    //IOField<vector> shadowPos(c.fieldIOobject("shadowPos", IOobject::NO_READ), np);//Jafari added
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const FrontTrackingParcel<ParcelType>& p = iter();

        active[i] = p.active();
        typeId[i] = p.typeId();
        currentIndex[i] = p.currentIndex();//Jafari added
        bubbleIndex[i] = p.bubbleIndex();//Jafari added
        //shadowPos[i] = p.shadowPos();//Jafari added
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();

        i++;
    }

    active.write();
    typeId.write();
    currentIndex.write();//Jafari added
    bubbleIndex.write();//Jafari added
    //shadowPos.write();//Jafari added
    nParticle.write();
    d.write();
    dTarget.write();
    U.write();
    rho.write();
    age.write();
    tTurb.write();
    UTurb.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const FrontTrackingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.active()
            << token::SPACE << p.typeId()
            << token::SPACE << p.currentIndex()
            << token::SPACE << p.bubbleIndex()
            //<< token::SPACE << p.shadowPos()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            sizeof(p.active())
          + sizeof(p.typeId())
          + sizeof(p.currentIndex())
          + sizeof(p.bubbleIndex())
          //+ sizeof(p.shadowPos())
          + sizeof(p.nParticle())
          + sizeof(p.d())
          + sizeof(p.dTarget())
          + sizeof(p.U())
          + sizeof(p.rho())
          + sizeof(p.age())
          + sizeof(p.tTurb())
          + sizeof(p.UTurb())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const FrontTrackingParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
