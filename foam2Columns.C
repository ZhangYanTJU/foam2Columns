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
#include "dynamicFvMesh.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

#include "ReadFields.H" // for the define of class::IOobjectList and function::readFields

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"
    timeSelector::addOptions();
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be processed. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );
    argList::addOption
    (
        "field",
        "word",
        "specify a field to be processed. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    instantList timeDirs = timeSelector::select0(runTime, args);

    List<word> selectedFieldNames;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFieldNames;
    }
    if (args.optionFound("field"))
    {
        word selectedFieldName;
        args.optionLookup("field")() >> selectedFieldName;
        selectedFieldNames.setSize(1, selectedFieldName);
    }

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Maintain a stack of the stored objects to clear after executing
        LIFOStack<regIOobject*> storedObjects;
        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());
        // Read Fields
        readFields<volScalarField>(mesh, objects, selectedFieldNames, storedObjects);

        Switch contamination = false;
        forAll (selectedFieldNames, fieldI)
        {
            if (!mesh.objectRegistry::foundObject<volScalarField>(selectedFieldNames[fieldI]))
            {
                contamination = true;
            }
        }
        if (contamination)
        {
            while (!storedObjects.empty())
            {
                storedObjects.pop()->checkOut();
            }
            continue;
        }

        PtrList<volScalarField> allFields(selectedFieldNames.size());
        forAll(allFields, fieldI)
        {
            allFields.set(fieldI, mesh.objectRegistry::lookupObject<volScalarField>(selectedFieldNames[fieldI]) );
            allFields[fieldI].writeMinMax(Info);
        }

        fileName outputFile = runTime.caseName().path() / runTime.timeName() / fileName("foam2Columns");
        forAll (selectedFieldNames, fieldI)
        {
            outputFile = outputFile + '_' + selectedFieldNames[fieldI];
        }
        Info<<"outputFile======="<<outputFile<<endl;
        OFstream foam2Columns(outputFile);

        //foam2Columns
        //    << "x" << "\t"
        //    << "y" << "\t"
        //    << "z" ;
        forAll (selectedFieldNames, fieldI)
        {
            foam2Columns << selectedFieldNames[fieldI] << "\t";
        }

        foam2Columns << endl;

        forAll (allFields[0], cellI)
        {
            //foam2Columns
            //    << mesh.C()[cellI].component(0) << "\t"
            //    << mesh.C()[cellI].component(1) << "\t"
            //    << mesh.C()[cellI].component(2);

            for (auto Field : allFields)
            {
                foam2Columns << Field[cellI] << "\t";
            }
            foam2Columns << endl;
        }

        while (!storedObjects.empty())
        {
            storedObjects.pop()->checkOut();
        }
    }

    Info<<"Finished!"<<endl;

    return 0;
}


// ************************************************************************* //
