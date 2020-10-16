/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#ifndef foam2ColumnsTemplates_H
#define foam2ColumnsTemplates_H

#include "OFstream.H"
#include "SortableList.H"
#include "IOobjectList.H"
#include "PtrList.H"
#include "Field.H"
#include "hashedWordList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<typename Type>
using VolFieldType = GeometricField<Type, fvPatchField, volMesh>;

template<class Type, template<typename> typename tmpType>
hashedWordList validFields
(
    const List<word>& fieldNames,
    const IOobjectList& objsList
);

template<template<typename> typename tmpType>
HashTable<label, word> validFields
(
    const List<word>& fieldNames,
    const IOobjectList& objsList
);

template<class Type>
void writeValue(OFstream& os, word sep, const scalar& value);

template<class Type>
void writeValue(OFstream& os, const Type& value);

template<class Type, template<typename> typename tmpType>
List<Type> readData
(
	const IOobject& obj,
    const volMesh::Mesh& mesh
);

template<template<typename> typename Type>
void readDatas
(
    const volMesh::Mesh& mesh,
    const IOobjectList& objects,
    const wordList& selectedFields,
    PtrList<List<scalar>>& data
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foam2ColumnsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
