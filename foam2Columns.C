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
    argList::noParallel();
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

                    // Extract the fields name list and the number of Components
                    HashTable<label> validFieldComponents = 
                    validFields<IOField>(selectedLagrangianFields, objects);
                 
                    // check completeness
                    if ( validFieldComponents.size() < selectedLagrangianFields.size())
                    {
                        Info << "At least one of the Lagrangian fields is missing, "
                             << nl << "break current time loop!" << endl;
                        continue;
                    }

                    // Head of writing files, each component has a header
                    wordList sortedHeader(0);

                    forAll(selectedLagrangianFields, fieldI)
                    {
                        const word& fieldName = selectedLagrangianFields[fieldI];
                        const label& nComponents = validFieldComponents[fieldName];
                        if( nComponents > 1 )
                        {
                            for (label i=0; i<nComponents; i++)
                            {
                                sortedHeader.append(fieldName + "_" + name(i));
                            }
                        }
                        else
                        {
                            sortedHeader.append(fieldName);
                        }
                        outputFile = outputFile + "_" + fieldName;
                    }

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

                    PtrList<List<scalar>> dataLagrangian(sortedHeader.size());
                    readDatas<IOField>(mesh, objects, selectedLagrangianFields, dataLagrangian);
                    passiveParticleCloud parcels(mesh, cloudDirs[i]);

                    label iParticle(0);
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

                        forAll(dataLagrangian, columnsI)
                        {
                            writeValue<scalar>
                            (
                                foam2ColumnsLagrangian,
                                sep,
                                dataLagrangian[columnsI][iParticle]
                            );
                        }

                        // Important, must be nl instead of endl, which is much slower
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

            // Extract the fields name list and the number of Components
            HashTable<label> validFieldComponents = 
            validFields<VolFieldType>(selectedEulerianFields, objects);
         
            // check completeness
            if ( validFieldComponents.size() < selectedEulerianFields.size())
            {
                Info << "At least one of the Eulerian fields is missing, "
                     << nl << "break current time loop!" << endl;
                continue;
            }

            // Head of writing files, each component has a header
            wordList sortedHeader(0);

            forAll(selectedEulerianFields, fieldI)
            {
                const word& fieldName = selectedEulerianFields[fieldI];
                const label& nComponents = validFieldComponents[fieldName];
                if( nComponents > 1 )
                {
                    for (label i=0; i<nComponents; i++)
                    {
                        sortedHeader.append(fieldName + "_" + name(i));
                    }
                }
                else
                {
                    sortedHeader.append(fieldName);
                }
                outputFile = outputFile + "_" + fieldName;
            }

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

            PtrList<List<scalar>> dataEulerian(sortedHeader.size());
            readDatas<VolFieldType>(mesh, objects, selectedEulerianFields, dataEulerian);

            forAll(mesh.C(), cellI)
            {
                foam2ColumnsEulerian
                    << mesh.C()[cellI].component(0) << sep
                    << mesh.C()[cellI].component(1) << sep
                    << mesh.C()[cellI].component(2);

                forAll(dataEulerian, columnsI)
                {
                    writeValue<scalar>
                    (
                        foam2ColumnsEulerian,
                        sep,
                        dataEulerian[columnsI][cellI]
                    );
                }
                // Important, must be nl instead of endl, which is much slower
                foam2ColumnsEulerian << nl;
            }
        }
    }
    Info<< "Finished!" << endl;

    return 0;
}


// ************************************************************************* //
