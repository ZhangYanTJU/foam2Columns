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
        "specify a list of volScalarFields to be processed. Eg, '(p T)' - "
        "regular expressions not currently supported"
    );
    argList::addOption
    (
        "vectorFields",
        "list",
        "specify a list of volVectorFields to be processed. Eg, '(U)' - "
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

    List<word> selectedScalarFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedScalarFields;
    }
    List<word> selectedVectorFields;
    if (args.optionFound("vectorFields"))
    {
        args.optionLookup("vectorFields")() >> selectedVectorFields;
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
                IOobjectList sprayObjs(runTime, runTime.timeName(), dir);

                if (sprayObjs.found(word("positions")))
                {
                    const word cloudName = cloudDirs[i];

                    Info<< "Converting lagrangian fields for cloud "
                        << cloudName << endl;

                    fileName outputFile = outputPath/fileName("Lagrangian");
                    outputFile = outputFile + '_' + cloudName;

                    forAll (selectedLagrangianFields, fieldI)
                    {
                        outputFile = outputFile + '_' + selectedLagrangianFields[fieldI];
                    }

                    Info<< "outputFile =" << outputFile << endl;
                    OFstream foam2ColumnsLagrangian(outputFile);

                    foam2ColumnsLagrangian
                        << "x" << "\t"
                        << "y" << "\t"
                        << "z" << "\t"
                        << "origId" << "\t"
                        << "origProc";

                    forAll (selectedLagrangianFields, fieldI)
                    {
                        foam2ColumnsLagrangian << "\t" << selectedLagrangianFields[fieldI];
                    }
                    foam2ColumnsLagrangian << nl;

                    PtrList<IOField<scalar>> allLagrangianFields(selectedLagrangianFields.size());
                    forAll(selectedLagrangianFields, fieldI)
                    {
                        const word fieldName = selectedLagrangianFields[fieldI];

                        // Check object on local mesh
                        IOobject fieldIOobject
                        (
                            fieldName,
                            runTime.timeName(),
                            dir,
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        );

                        if (fieldIOobject.typeHeaderOk<IOField<scalar>>(true))
                        {
                            allLagrangianFields.set
                            (
                                fieldI,
                                new
                                IOField<scalar>
                                (
                                    fieldIOobject
                                )
                            );
                        }
                    }

                    passiveParticleCloud parcels(mesh, cloudDirs[i]);
                    label iParticle(0);
                    forAllConstIter(passiveParticleCloud, parcels, iter)
                    {
                        // Loop for each particle
                        const passiveParticle& ppi = iter();

                        foam2ColumnsLagrangian
                            << ppi.position().component(0) << "\t"
                            << ppi.position().component(1) << "\t"
                            << ppi.position().component(2) << "\t"
                            << ppi.origId() << "\t"
                            << ppi.origProc();

                        forAll(allLagrangianFields, fieldI)
                        {
                            foam2ColumnsLagrangian << "\t" << allLagrangianFields[fieldI][iParticle];
                        }
                        foam2ColumnsLagrangian << nl;
                        iParticle ++;
                    }
                }
            }
        }

        // For Eulerian scalarFields foam2Columns
        if(selectedScalarFields.size()||selectedVectorFields.size())
        {
            // Maintain a stack of the stored objects to clear after executing
            LIFOStack<regIOobject*> storedObjects;
            // Read objects in time directory
            IOobjectList objects(mesh, runTime.timeName());
            // Read Fields
            readFields<volScalarField>(mesh, objects, selectedScalarFields, storedObjects);
            readFields<volVectorField>(mesh, objects, selectedVectorFields, storedObjects);

            Switch contamination = false;
            forAll (selectedScalarFields, fieldI)
            {
                if (!mesh.objectRegistry::foundObject<volScalarField>(selectedScalarFields[fieldI]))
                {
                    contamination = true;
                    Info<< selectedScalarFields[fieldI] << " not found" << endl;
                }
            }
            forAll (selectedVectorFields, fieldI)
            {
                if (!mesh.objectRegistry::foundObject<volVectorField>(selectedVectorFields[fieldI]))
                {
                    contamination = true;
                    Info<< selectedVectorFields[fieldI] << " not found" << endl;
                }
            }
            if (contamination)
            {
                while (!storedObjects.empty())
                {
                    storedObjects.pop()->checkOut();
                }
                Info << "At least one of the Eulerian fields is missing, "
                     << nl << "break current time loop!" << endl;
                continue;
            }

            PtrList<volScalarField> allScalarFields(selectedScalarFields.size());
            forAll(allScalarFields, fieldI)
            {
                allScalarFields.set(fieldI, mesh.objectRegistry::lookupObject<volScalarField>(selectedScalarFields[fieldI]) );
                allScalarFields[fieldI].writeMinMax(Info);
            }
            PtrList<volVectorField> allVectorFields(selectedVectorFields.size());
            forAll(allVectorFields, fieldI)
            {
                allVectorFields.set(fieldI, mesh.objectRegistry::lookupObject<volVectorField>(selectedVectorFields[fieldI]) );
                allVectorFields[fieldI].writeMinMax(Info);
            }

            fileName outputFile = outputPath/fileName("Eulerian");
            forAll (selectedScalarFields, fieldI)
            {
                outputFile = outputFile + '_' + selectedScalarFields[fieldI];
            }
            forAll (selectedVectorFields, fieldI)
            {
                outputFile = outputFile + '_' + selectedVectorFields[fieldI];
            }
            Info<< "outputFile =" << outputFile << endl;
            OFstream foam2Columns(outputFile);

            foam2Columns
                << "x" << "\t"
                << "y" << "\t"
                << "z" ;
            forAll (selectedScalarFields, fieldI)
            {
                foam2Columns << "\t" << selectedScalarFields[fieldI];
            }
            forAll (selectedVectorFields, fieldI)
            {
                foam2Columns
                    << "\t" << selectedVectorFields[fieldI] << "_x"
                    << "\t" << selectedVectorFields[fieldI] << "_y"
                    << "\t" << selectedVectorFields[fieldI] << "_z";
            }
            foam2Columns << nl;

            forAll (allScalarFields[0], cellI)
            {
                foam2Columns
                    << mesh.C()[cellI].component(0) << "\t"
                    << mesh.C()[cellI].component(1) << "\t"
                    << mesh.C()[cellI].component(2);

                forAll(allScalarFields, fieldI)
                {
                    foam2Columns << "\t" << allScalarFields[fieldI][cellI];
                }
                forAll(allVectorFields, fieldI)
                {
                    foam2Columns
                        << "\t" << allVectorFields[fieldI][cellI].component(0)
                        << "\t" << allVectorFields[fieldI][cellI].component(1)
                        << "\t" << allVectorFields[fieldI][cellI].component(2);
                }
                foam2Columns << nl;
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
