/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "foam2ColumnsTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<typename> typename tmpType>
Foam::hashedWordList
Foam::validFields
(
    const List<word>& fieldNames,
    const IOobjectList& objsList
)
{
    hashedWordList validFieldNames;
    IOobjectList objects
    (
        objsList.lookupClass
        (
            tmpType<Type>::typeName
        )
    );

    forAll(fieldNames, iName)
    {
        const word& fieldName = fieldNames[iName];        
        if (objects.lookup(fieldName) != nullptr)
        {
            validFieldNames.append(fieldName);
        }
    }
    return validFieldNames;
}


template<template<typename> typename tmpType>
HashTable<label, word>
Foam::validFields
(
    const List<word>& fieldNames,
    const IOobjectList& objsList
)
{
    HashTable<label, word> validFieldComponents;
    // Extract the fields name list        
    const hashedWordList scalarFieldNames = 
    validFields<scalar, tmpType>(fieldNames, objsList);
    const hashedWordList vectorFieldNames = 
    validFields<vector, tmpType>(fieldNames, objsList);
    const hashedWordList sphericalTensorFieldNames = 
    validFields<sphericalTensor, tmpType>(fieldNames, objsList);
    const hashedWordList symmTensorFieldNames =
    validFields<symmTensor, tmpType>(fieldNames, objsList);
    const hashedWordList tensorFieldNames = 
    validFields<tensor, tmpType>(fieldNames, objsList);

    forAll(scalarFieldNames, fieldI)
    {
        const label nCpnts = pTraits<scalar>::nComponents;
        const word& fieldName = scalarFieldNames[fieldI];
        validFieldComponents.insert(fieldName, nCpnts);
    }
    forAll(vectorFieldNames, fieldI)
    {
        const label nCpnts = pTraits<vector>::nComponents;;
        const word& fieldName = vectorFieldNames[fieldI];
        validFieldComponents.insert(fieldName, nCpnts);
    }
    forAll(sphericalTensorFieldNames, fieldI)
    {
        const label nCpnts = pTraits<sphericalTensor>::nComponents;
        const word& fieldName = sphericalTensorFieldNames[fieldI];
        validFieldComponents.insert(fieldName, nCpnts);
    }
    forAll(symmTensorFieldNames, fieldI)
    {
        const label nCpnts = pTraits<symmTensor>::nComponents;
        const word& fieldName = symmTensorFieldNames[fieldI];
        validFieldComponents.insert(fieldName, nCpnts);
    }
    forAll(tensorFieldNames, fieldI)
    {
        const label nCpnts = pTraits<tensor>::nComponents;
        const word& fieldName = tensorFieldNames[fieldI];
        validFieldComponents.insert(fieldName, nCpnts);
    }

    return validFieldComponents;
}


template<class Type>
void Foam::writeValue(OFstream& os, word sep, const scalar& value)
{
    os << sep << value;
}

template<class Type>
void Foam::writeValue(OFstream& os, word sep, const Type& value)
{
    for (label i=0; i<pTraits<Type>::nComponents; i++)
    {
        os  << sep << value.component(i);
    }
}


template<class Type, template<typename> typename tmpType>
List<Type> Foam::readData
(
	const IOobject& obj,
    const volMesh::Mesh& mesh
)
{
    if(tmpType<Type>::typeName == VolFieldType<Type>::typeName)
    {
        return VolFieldType<Type>(obj, mesh);
    }
    return IOField<Type>(obj);
}


template<template<typename> typename tmpType>
void Foam::readDatas
(
    const volMesh::Mesh& mesh,
	const IOobjectList& objects,
	const wordList& selectedFields,
	PtrList<List<scalar>>& data
)
{
	const hashedWordList scalarFieldNames = 
	validFields<scalar, tmpType>(selectedFields, objects);
	const hashedWordList vectorFieldNames = 
	validFields<vector, tmpType>(selectedFields, objects);
	const hashedWordList sphericalTensorFieldNames = 
	validFields<sphericalTensor, tmpType>(selectedFields, objects);
	const hashedWordList symmTensorFieldNames =
	validFields<symmTensor, tmpType>(selectedFields, objects);
	const hashedWordList tensorFieldNames = 
	validFields<tensor, tmpType>(selectedFields, objects);
	
	label nColumns(0);
	forAll(selectedFields, fieldI)
	{
		const word& fieldName = selectedFields[fieldI];
		const IOobject& obj = *objects.lookup(fieldName);

		if(scalarFieldNames.found(fieldName))
		{
            List<scalar> tmpField = readData<scalar, tmpType>(obj, mesh);
			data.set(nColumns, new List<scalar>(tmpField.size()));
			forAll(tmpField, cellI)
			{
				data[nColumns][cellI] = tmpField[cellI];
			}
			nColumns++;
		}
		else if(vectorFieldNames.found(fieldName))
		{
            List<vector> tmpField = readData<vector, tmpType>(obj, mesh);
			for (label i=0; i<pTraits<vector>::nComponents; i++)
			{
				data.set(nColumns, new List<scalar>(tmpField.size()));
				forAll(tmpField, cellI)
				{
					data[nColumns][cellI] = tmpField[cellI].component(i);
				}
				nColumns++;
			}
		}
		else if(sphericalTensorFieldNames.found(fieldName))
		{
            List<sphericalTensor> tmpField = readData<sphericalTensor, tmpType>(obj, mesh);
			for (label i=0; i<pTraits<sphericalTensor>::nComponents; i++)
			{
				data.set(nColumns, new List<scalar>(tmpField.size()));
				forAll(tmpField, cellI)
				{
					data[nColumns][cellI] = tmpField[cellI].component(i);
				}
				nColumns++;
			}
		}
		else if(symmTensorFieldNames.found(fieldName))
		{
            List<symmTensor> tmpField = readData<symmTensor, tmpType>(obj, mesh);
			for (label i=0; i<pTraits<symmTensor>::nComponents; i++)
			{
				data.set(nColumns, new List<scalar>(tmpField.size()));
				forAll(tmpField, cellI)
				{
					data[nColumns][cellI] = tmpField[cellI].component(i);
				}
				nColumns++;
			}
		}
		else if(tensorFieldNames.found(fieldName))
		{
            List<tensor> tmpField = readData<tensor, tmpType>(obj, mesh);
			for (label i=0; i<pTraits<tensor>::nComponents; i++)
			{
				data.set(nColumns, new List<scalar>(tmpField.size()));
				forAll(tmpField, cellI)
				{
					data[nColumns][cellI] = tmpField[cellI].component(i);
				}
				nColumns++;
			}
		}
        else
        {
            FatalErrorInFunction
                << "Unable to read field " << fieldName
                << abort(FatalError);
        }
	}
}


// ************************************************************************* //
