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

Application


Description


\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "OFstream.H"
#include "ReadFields.H" // for the define of class::IOobjectList and function::readFields
#include "cloud.H" //cloud::prefix
#include "passiveParticleCloud.H"
#include "foam2ColumnsTemplates.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"
    timeSelector::addOptions();
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of volFields to be processed. Eg, '(p U T)' - "
        "regular expressions not currently supported"
    );
    argList::addOption
    (
        "lagrangianFields",
        "list",
        "specify a list of lagrangian fields to be reconstructed. Eg, '(d T)' -"
        "only scalar field is supported currently, "
        "positions always included."
    );
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    instantList timeDirs = timeSelector::select0(runTime, args);

    List<word> selectedEulerianFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedEulerianFields;
    }
    List<word> selectedLagrangianFields;
    if (args.optionFound("lagrangianFields"))
    {
        args.optionLookup("lagrangianFields")() >> selectedLagrangianFields;
    }

    fileName folderName("foam2Columns");
    fileName postProcessPath = runTime.caseName().path()/"postProcessing"/folderName;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        fileName outputPath = postProcessPath/runTime.timeName();
        mkDir(outputPath);
        #include "createMesh.H"
        Info<< "mesh.nCells() = " << mesh.nCells() << endl;
        const word sep = "\t";
        // For Lagrangian scalarFields foam2Columns
        if(selectedLagrangianFields.size())
        {
            fileName lagrangianDir = runTime.caseName().path()/runTime.timeName()/cloud::prefix;

            fileNameList cloudDirs;
            if (!lagrangianDir.empty())
            {
                cloudDirs = readDir
                (
                    lagrangianDir,
                    fileType::directory
                );
            }

            forAll(cloudDirs, i)
            {
                fileName dir(cloud::prefix/cloudDirs[i]);

                // Do local scan for valid cloud objects
                IOobjectList objects(runTime, runTime.timeName(), dir);

                if (objects.found(word("positions")))
                {
                    const word cloudName = cloudDirs[i];

                    Info<< "Converting lagrangian fields for cloud "
                        << cloudName << endl;

                    fileName outputFile = outputPath/fileName("Lagrangian");
                    outputFile = outputFile + '_' + cloudName;

                    // Extract the fields name list        
                    const hashedWordList scalarFieldNames = 
                    validFields<scalar, IOField>
                    (selectedLagrangianFields, objects);
                    const hashedWordList vectorFieldNames = 
                    validFields<vector, IOField>
                    (selectedLagrangianFields, objects);
                    const hashedWordList sphericalTensorFieldNames = 
                    validFields<sphericalTensor, IOField>
                    (selectedLagrangianFields, objects);
                    const hashedWordList symmTensorFieldNames =
                    validFields<symmTensor, IOField>
                    (selectedLagrangianFields, objects);
                    const hashedWordList tensorFieldNames = 
                    validFields<tensor, IOField>
                    (selectedLagrangianFields, objects);


                    HashTable<label> validFieldComponents = validFields<IOField>(selectedLagrangianFields, objects);
                   
                 
                    // check completeness
                    if ( validFieldComponents.size() < selectedLagrangianFields.size())
                    {
                        Info << "At least one of the Lagrangian fields is missing, "
                             << nl << "break current time loop!" << endl;
                        continue;
                    }



                    wordList sortedHeader(0);
                    word sortedExtensions("");

/*
HashTable<using> aaa;
using typeName = scalar;
using typeName = scalar;//typename typeNameTmp;
*/
                    forAll(selectedLagrangianFields, fieldI)
                    {
                        const word fieldName = selectedLagrangianFields[fieldI];
                        const label nComponents = validFieldComponents[fieldName];
                        for (label i=0; i<nComponents; i++)
                        {
                            sortedHeader.append(fieldName + "_" + name(i));
                        }
                        outputFile = outputFile + "_" + fieldName;
                    }





                    constructFile<scalar>(scalarFieldNames, sortedExtensions, sortedHeader);
                    constructFile<vector>(vectorFieldNames, sortedExtensions, sortedHeader);
                    constructFile<sphericalTensor>(sphericalTensorFieldNames, sortedExtensions, sortedHeader);
                    constructFile<symmTensor>(symmTensorFieldNames, sortedExtensions, sortedHeader);
                    constructFile<tensor>(tensorFieldNames, sortedExtensions, sortedHeader);

PtrList<List<scalar>> data(sortedHeader.size());
readDatas<IOField>(objects, selectedLagrangianFields, data);

                    outputFile = outputFile + sortedExtensions;

                    Info<< "outputFile =" << outputFile << endl;
                    OFstream foam2ColumnsLagrangian(outputFile);

                    foam2ColumnsLagrangian
                        << "x" << sep
                        << "y" << sep
                        << "z" << sep
                        << "origId" << sep
                        << "origProc";

                    forAll(sortedHeader, iHeader)
                    {
                        foam2ColumnsLagrangian << sep << sortedHeader[iHeader];
                    }

                    foam2ColumnsLagrangian << nl;

                    PtrList<Field<scalar>> allScalarFields(scalarFieldNames.size());
                    forAll(allScalarFields, fieldI)
                    {
                        readField<scalar>
                        (
                            allScalarFields[fieldI],
                            scalarFieldNames[fieldI],
                            objects
                        );
                    }

                    PtrList<Field<vector>> allVectorFields(vectorFieldNames.size());
                    forAll(allVectorFields, fieldI)
                    {
                        readField<vector>
                        (
                            allVectorFields[fieldI],
                            vectorFieldNames[fieldI],
                            objects
                        );
                    }
                    PtrList<Field<sphericalTensor>> allSphericalTensorFields(sphericalTensorFieldNames.size());
                    forAll(allSphericalTensorFields, fieldI)
                    {
                        readField<sphericalTensor>
                        (
                            allSphericalTensorFields[fieldI],
                            sphericalTensorFieldNames[fieldI],
                            objects
                        );
                    }
                    PtrList<Field<symmTensor>> allSymmTensorFields(symmTensorFieldNames.size());
                    forAll(allSymmTensorFields, fieldI)
                    {
                        readField<symmTensor>
                        (
                            allSymmTensorFields[fieldI],
                            symmTensorFieldNames[fieldI],
                            objects
                        );
                    }
                    PtrList<Field<tensor>> allTensorFields(tensorFieldNames.size());
                    forAll(allTensorFields, fieldI)
                    {
                        readField<tensor>
                        (
                            allTensorFields[fieldI],
                            tensorFieldNames[fieldI],
                            objects
                        );
                    }

                    passiveParticleCloud parcels(mesh, cloudDirs[i]);

                    forAllConstIter(passiveParticleCloud, parcels, iter)
                    {
                        // Loop for each particle
                        const passiveParticle& ppi = iter();

                        foam2ColumnsLagrangian
                            << ppi.position().component(0) << sep
                            << ppi.position().component(1) << sep
                            << ppi.position().component(2) << sep
                            << ppi.origId() << sep
                            << ppi.origProc();

                        label iParticle(0);

                        forAll(allScalarFields, fieldI)
                        {
                            writeValue<scalar>(foam2ColumnsLagrangian, sep, allScalarFields[fieldI][iParticle]);
                        }
                        forAll(allVectorFields, fieldI)
                        {
                            writeValue<vector>(foam2ColumnsLagrangian, sep, allVectorFields[fieldI][iParticle]);
                        }
                        forAll(allSphericalTensorFields, fieldI)
                        {
                            writeValue<sphericalTensor>(foam2ColumnsLagrangian, sep, allSphericalTensorFields[fieldI][iParticle]);
                        }
                        forAll(allSymmTensorFields, fieldI)
                        {
                            writeValue<symmTensor>(foam2ColumnsLagrangian, sep, allSymmTensorFields[fieldI][iParticle]);
                        }
                        forAll(allTensorFields, fieldI)
                        {
                            writeValue<tensor>(foam2ColumnsLagrangian, sep, allTensorFields[fieldI][iParticle]);
                        }

                        foam2ColumnsLagrangian << nl;
                        iParticle ++;
                    }
                }
            }
        }

        // For Eulerian scalarFields foam2Columns
        if(selectedEulerianFields.size())
        {
            fileName outputFile = outputPath/fileName("Eulerian");

            // Read objects in time directory
            IOobjectList objects(mesh, runTime.timeName());

            // Extract the fields name list        
            const hashedWordList scalarFieldNames = 
            validFields<scalar, VolFieldType>
            (selectedEulerianFields, objects);
            const hashedWordList vectorFieldNames = 
            validFields<vector, VolFieldType>
            (selectedEulerianFields, objects);
            const hashedWordList sphericalTensorFieldNames = 
            validFields<sphericalTensor, VolFieldType>
            (selectedEulerianFields, objects);
            const hashedWordList symmTensorFieldNames =
            validFields<symmTensor, VolFieldType>
            (selectedEulerianFields, objects);
            const hashedWordList tensorFieldNames = 
            validFields<tensor, VolFieldType>
            (selectedEulerianFields, objects);

            // Maintain a stack of the stored objects to clear after executing
            LIFOStack<regIOobject*> storedObjects;

            // Read Fields
            readFields<volScalarField>
            (mesh, objects, scalarFieldNames, storedObjects);
            readFields<volVectorField>
            (mesh, objects, vectorFieldNames, storedObjects);
            readFields<volSphericalTensorField>
            (mesh, objects, sphericalTensorFieldNames, storedObjects);
            readFields<volSymmTensorField>
            (mesh, objects, symmTensorFieldNames, storedObjects);
            readFields<volTensorField>
            (mesh, objects, tensorFieldNames, storedObjects);

            // check completeness
            if ( storedObjects.size() < selectedEulerianFields.size())
            {
                while (!storedObjects.empty())
                {
                    storedObjects.pop()->checkOut();
                }
                Info << "At least one of the Eulerian fields is missing, "
                     << nl << "break current time loop!" << endl;
                continue;
            }

            wordList sortedHeader(0);
            word sortedExtensions("");

            constructFile<scalar>(scalarFieldNames, sortedExtensions, sortedHeader);
            constructFile<vector>(vectorFieldNames, sortedExtensions, sortedHeader);
            constructFile<sphericalTensor>(sphericalTensorFieldNames, sortedExtensions, sortedHeader);
            constructFile<symmTensor>(symmTensorFieldNames, sortedExtensions, sortedHeader);
            constructFile<tensor>(tensorFieldNames, sortedExtensions, sortedHeader);


PtrList<List<scalar>> data(sortedHeader.size());
//readDatas<VolFieldType>(objects, selectedEulerianFields, data);
readDatas<IOField>(objects, selectedEulerianFields, data);




            outputFile = outputFile + sortedExtensions;

            Info<< "outputFile =" << outputFile << endl;
            OFstream foam2ColumnsEulerian(outputFile);

            foam2ColumnsEulerian
                << "x" << sep
                << "y" << sep
                << "z" ;

            forAll(sortedHeader, iHeader)
            {
                foam2ColumnsEulerian << sep << sortedHeader[iHeader];
            }

            foam2ColumnsEulerian << nl;

            PtrList<volScalarField> allScalarFields(scalarFieldNames.size());
            forAll(allScalarFields, fieldI)
            {
                allScalarFields.set(fieldI, mesh.objectRegistry::lookupObject<volScalarField>(scalarFieldNames[fieldI]) );
                allScalarFields[fieldI].writeMinMax(Info);
            }
            PtrList<volVectorField> allVectorFields(vectorFieldNames.size());
            forAll(allVectorFields, fieldI)
            {
                allVectorFields.set(fieldI, mesh.objectRegistry::lookupObject<volVectorField>(vectorFieldNames[fieldI]) );
                allVectorFields[fieldI].writeMinMax(Info);
            }
            PtrList<volSphericalTensorField> allSphericalTensorFields(sphericalTensorFieldNames.size());
            forAll(allSphericalTensorFields, fieldI)
            {
                allSphericalTensorFields.set(fieldI, mesh.objectRegistry::lookupObject<volSphericalTensorField>(sphericalTensorFieldNames[fieldI]) );
                allSphericalTensorFields[fieldI].writeMinMax(Info);
            }
            PtrList<volSymmTensorField> allSymmTensorFields(symmTensorFieldNames.size());
            forAll(allSymmTensorFields, fieldI)
            {
                allSymmTensorFields.set(fieldI, mesh.objectRegistry::lookupObject<volSymmTensorField>(symmTensorFieldNames[fieldI]) );
                allSymmTensorFields[fieldI].writeMinMax(Info);
            }
            PtrList<volTensorField> allTensorFields(tensorFieldNames.size());
            forAll(allTensorFields, fieldI)
            {
                allTensorFields.set(fieldI, mesh.objectRegistry::lookupObject<volTensorField>(tensorFieldNames[fieldI]) );
                allTensorFields[fieldI].writeMinMax(Info);
            }

            forAll(mesh.C(), cellI)
            {
                foam2ColumnsEulerian
                    << mesh.C()[cellI].component(0) << sep
                    << mesh.C()[cellI].component(1) << sep
                    << mesh.C()[cellI].component(2);

                forAll(allScalarFields, fieldI)
                {
                    writeValue<scalar>(foam2ColumnsEulerian, sep, allScalarFields[fieldI][cellI]);
                }
                forAll(allVectorFields, fieldI)
                {
                    writeValue<vector>(foam2ColumnsEulerian, sep, allVectorFields[fieldI][cellI]);
                }
                forAll(allSphericalTensorFields, fieldI)
                {
                    writeValue<sphericalTensor>(foam2ColumnsEulerian, sep, allSphericalTensorFields[fieldI][cellI]);
                }
                forAll(allSymmTensorFields, fieldI)
                {
                    writeValue<symmTensor>(foam2ColumnsEulerian, sep, allSymmTensorFields[fieldI][cellI]);
                }
                forAll(allTensorFields, fieldI)
                {
                    writeValue<tensor>(foam2ColumnsEulerian, sep, allTensorFields[fieldI][cellI]);
                }

                foam2ColumnsEulerian << nl;
            }

            while (!storedObjects.empty())
            {
                storedObjects.pop()->checkOut();
            }
        }
    }
    Info<< "Finished!" << endl;

    return 0;
}


// ************************************************************************* //
