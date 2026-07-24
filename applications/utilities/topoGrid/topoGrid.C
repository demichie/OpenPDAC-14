/*---------------------------------------------------------------------------* \
License
    This file is part of OpenPDAC.

    OpenPDAC is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenPDAC is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenPDAC.  If not, see <http://www.gnu.org/licenses/>.

Application
    topoGrid

Description
    Deforms a polyMesh using an ESRI raster ascii file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "vector.H"
#include "pointFields.H"
#include "IStringStream.H"
#include "volPointInterpolation.H"
#include "UniformTable2.H"
#include "RectangularMatrix.H"
#include <fstream>
#include <sstream>
#include "IOstreams.H"
#include <cstring>
#include "Pstream.H"
#include "tetPointRef.H"
#include "OFstream.H"
#include "globalIndex.H"
#include "polyMeshCheck.H"
#include "syncTools.H"
#include "vectorList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool isProcessorFace(const Foam::polyMesh& mesh, const Foam::label faceI)
{
    // Check if face is internal
    if (faceI < mesh.nInternalFaces())
    {
        return false;
    }

    // If not internal, get its patch ID
    const Foam::label patchID = mesh.boundaryMesh().whichPatch(faceI);

    // Check if the ID is valid
    if (patchID == -1)
    {
        // This should not happen for a valid boundary face
        return false;
    }

    // Get the patch
    const Foam::polyPatch& pp = mesh.boundaryMesh()[patchID];

    // Check if the patch is a processorPolyPatch
    return Foam::isA<Foam::processorPolyPatch>(pp);
}

// =========================================================================

/*
 * Calculate the normalization factors for each edge (1.0 / N),
 * where N is the number of processors sharing the edge.
 * This value is a topological invariant.
 */
scalarList calculateEdgeNormFactors(const fvMesh& mesh)
{
    // 1. Count the processors for each edge using synchronization.
    labelList nProcsPerEdge(mesh.nEdges(), 1);

    syncTools::syncEdgeList(mesh, nProcsPerEdge, plusEqOp<label>(), 0);

    // 2. Convert the count into the normalization factor.
    scalarList normFactors(mesh.nEdges());
    forAll(nProcsPerEdge, edgeI)
    {
        normFactors[edgeI] = 1.0 / scalar(nProcsPerEdge[edgeI]);
    }

    return normFactors;
}

// =========================================================================

/*
 * Computes the physical number of neighbors for each point.
 * This is another topological invariant.
 * Uses pre-calculated edge normalization factors.
 */
scalarList calculateNeighbourCount(const fvMesh& mesh,
                                   const scalarList& normFactors)
{
    // 1. Accumulate fractions locally
    scalarList tempNeighbourCount(mesh.nPoints(), 0.0);
    forAll(mesh.edges(), edgeI)
    {
        const edge& e = mesh.edges()[edgeI];
        tempNeighbourCount[e.start()] += normFactors[edgeI];
        tempNeighbourCount[e.end()] += normFactors[edgeI];
    }

    // 2. Synchronize to add fractions. The result will be the correct
    // count.
    syncTools::syncPointList(mesh, tempNeighbourCount, plusEqOp<scalar>(), 0.0);

    return tempNeighbourCount;
}

// =========================================================================

/*
 * Count how many processors see each MESH POINT ---
 * Returns a labelList indexed by local mesh point ID, containing how many procs
 * share it.
 */
labelList countProcessorsPerMeshPoint(const fvMesh& mesh)
{
    // 1. Initialize a labelList to keep counts.
    // Each local point contributes 1 to its count.
    labelList nProcsPerPoint(mesh.nPoints(), 1);

    // 2. Run sync to add these counts for the shared
    // points.
    syncTools::syncPointList(mesh, nProcsPerPoint, plusEqOp<label>(), 0);

    return nProcsPerPoint;
}

// =========================================================================

scalarList calculateGlobalCellCount(const fvMesh& mesh)
{
    // Accumula localmente il conteggio delle celle
    scalarList tempCellCount(mesh.nPoints(), 0.0);
    forAll(mesh.points(), pointI)
    {
        tempCellCount[pointI] = scalar(mesh.pointCells(pointI).size());
    }

    // Sincronizza per sommare i conteggi. Il risultato sarà il conteggio
    // corretto.
    syncTools::syncPointList(mesh,
                             tempCellCount,
                             plusEqOp<scalar>(),
                             0.0 // NULLVALUE
    );

    return tempCellCount;
}

// =========================================================================
// Funzioni di Smoothing (calcolano solo le posizioni target)
// =========================================================================

// Calcola le posizioni target secondo lo smoothing basato sui punti
tmp<pointField> getLaplacianTargetPoints(const fvMesh& mesh,
                                         const scalarList& normFactors,
                                         const scalarList& neighbourCount)
{
    // --- 1. Accumulo Frazionario Locale ---
    vectorField sumNeighbourCoords(mesh.nPoints(), vector::zero);
    forAll(mesh.edges(), edgeI)
    {
        const edge& e = mesh.edges()[edgeI];
        const label p0 = e.start();
        const label p1 = e.end(); // Corretto: label p1, non scalar p1
        const scalar normFactor = normFactors[edgeI];

        sumNeighbourCoords[p0] += mesh.points()[p1] * normFactor;
        sumNeighbourCoords[p1] += mesh.points()[p0] * normFactor;
    }

    // --- 2. Sincronizzazione Globale dei Contributi ---
    syncTools::syncPointList(mesh,
                             sumNeighbourCoords,
                             plusEqOp<vector>(),
                             vector::zero // NULLVALUE
    );

    // --- 3. Calcolo dei Punti Target Laplaciani ---
    // Pre-allocare la dimensione per il nuovo campo di punti.
    tmp<pointField> tTargetPoints(new pointField(mesh.points().size()));
    // --- CORREZIONE QUI: Usare const_cast per rimuovere il qualificatore const
    // ---
    pointField& targetPoints = const_cast<pointField&>(tTargetPoints());

    forAll(targetPoints, pointI) // Itera su tutti i punti locali
    {
        if (neighbourCount[pointI] > VSMALL) // Verifica che ci siano vicini
        {
            targetPoints[pointI] =
                sumNeighbourCoords[pointI] / neighbourCount[pointI];
        }
        else
        {
            // Se non ci sono vicini, il punto rimane nella sua posizione
            // attuale.
            targetPoints[pointI] = mesh.points()[pointI];
        }
    }

    syncTools::syncPointPositions(
        mesh, targetPoints, minEqOp<point>(), point(great, great, great));

    return tTargetPoints;
}

tmp<pointField> getCentroidalTargetPoints(const fvMesh& mesh)
{
    // --- 1. Accumulo Locale della Somma dei Centri delle Celle ---
    vectorField sumCellCentres(mesh.nPoints(), vector::zero);
    labelField localCellCounts(mesh.nPoints(), 0);

    forAll(mesh.points(), pointI)
    {
        const labelList& pCells = mesh.pointCells(pointI);
        forAll(pCells, i)
        {
            sumCellCentres[pointI] += mesh.cellCentres()[pCells[i]];
            localCellCounts[pointI] += 1; // Incrementa il conteggio locale
        }
    }

    // --- 2. Sincronizzazione Globale dei Contributi ---
    syncTools::syncPointList(
        mesh, sumCellCentres, plusEqOp<vector>(), vector::zero);
    // NUOVO: Sincronizza anche i conteggi delle celle
    syncTools::syncPointList(mesh,
                             localCellCounts,
                             plusEqOp<label>(),
                             0 // NULLVALUE
    );


    // --- 3. Calcolo dei Punti Target Centroidali ---
    tmp<pointField> tTargetPoints(new pointField(mesh.points().size()));
    pointField& targetPoints = const_cast<pointField&>(tTargetPoints());

    forAll(targetPoints, pointI)
    {
        if (localCellCounts[pointI]
            > VSMALL) // Verifica con il conteggio globale sincronizzato
        {
            targetPoints[pointI] =
                sumCellCentres[pointI] / localCellCounts[pointI];
        }
        else
        {
            targetPoints[pointI] = mesh.points()[pointI];
        }
    }
    return tTargetPoints;
}


// =========================================================================
// Lightweight mesh-quality metrics used by the optional quality-aware
// Laplacian smoother.  These checks are intentionally local/simple and are
// not meant to replace checkMesh; they only prevent the smoother from
// accepting a Laplacian step that worsens the most common checkMesh failures.
// =========================================================================

label countConcaveCellsSimple(const fvMesh& mesh, const scalar concavityTol)
{
    const cellList& cells = mesh.cells();
    const vectorField& faceAreas = mesh.faceAreas();
    const vectorField& faceCentres = mesh.faceCentres();
    const labelUList& owner = mesh.faceOwner();

    label localConcave = 0;

    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];
        bool isConcave = false;

        forAll(c, i)
        {
            const label faceI = c[i];
            vector n = faceAreas[faceI];

            // faceAreas are oriented out of the owner cell.  Reverse them
            // when this cell sees the face as a neighbour-side face.
            if (owner[faceI] != cellI)
            {
                n = -n;
            }

            const scalar nMag = mag(n);
            if (nMag <= VSMALL)
            {
                continue;
            }
            n /= nMag;

            const point& fCi = faceCentres[faceI];

            forAll(c, j)
            {
                if (j == i)
                {
                    continue;
                }

                const label otherFaceI = c[j];
                const vector v = faceCentres[otherFaceI] - fCi;

                if ((v & n) > concavityTol)
                {
                    isConcave = true;
                    break;
                }
            }

            if (isConcave)
            {
                break;
            }
        }

        if (isConcave)
        {
            ++localConcave;
        }
    }

    reduce(localConcave, sumOp<label>());
    return localConcave;
}


scalar minFaceTetQualitySimple(const fvMesh& mesh)
{
    const pointField& pts = mesh.points();
    const faceList& faces = mesh.faces();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();
    const labelUList& owner = mesh.faceOwner();
    const labelUList& neighbour = mesh.faceNeighbour();

    scalar localMinQuality = GREAT;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];
        if (f.size() < 3)
        {
            continue;
        }

        const point& fC = faceCentres[faceI];

        forAll(f, fp)
        {
            const point& p0 = pts[f[fp]];
            const point& p1 = pts[f[(fp + 1) % f.size()]];
            const scalar edgeLen = mag(p1 - p0);

            if (edgeLen <= VSMALL)
            {
                localMinQuality = 0.0;
                continue;
            }

            const vector triNormal = 0.5 * ((p0 - fC) ^ (p1 - fC));
            const scalar triArea = mag(triNormal);

            if (triArea <= VSMALL)
            {
                localMinQuality = 0.0;
                continue;
            }

            const label ownI = owner[faceI];
            const scalar ownHeight =
                mag((fC - cellCentres[ownI]) & (triNormal / triArea));
            const scalar ownVol = triArea * ownHeight / 3.0;
            const scalar ownScale = max(pow(edgeLen, 3), VSMALL);
            localMinQuality = min(localMinQuality, ownVol / ownScale);

            if (faceI < mesh.nInternalFaces())
            {
                const label neiI = neighbour[faceI];
                const scalar neiHeight =
                    mag((fC - cellCentres[neiI]) & (triNormal / triArea));
                const scalar neiVol = triArea * neiHeight / 3.0;
                const scalar neiScale = max(pow(edgeLen, 3), VSMALL);
                localMinQuality = min(localMinQuality, neiVol / neiScale);
            }
        }
    }

    reduce(localMinQuality, minOp<scalar>());
    return localMinQuality;
}


struct InterpolationSource
{
    point pCoords;    // Coordinates of the source point (x, y, z)
    scalar dz;        // Z-deformation value
    scalar dx;        // X-deformation value
    scalar dy;        // Y-deformation value
    scalar area;      // Associated area (e.g., face area for top patches)
    bool isTopCentre; // Flag to distinguish z=0 points from top face centers

    // Default constructor
    InterpolationSource()
    : pCoords(point::zero), dz(0.0), dx(0.0), dy(0.0), area(0.0),
      isTopCentre(false)
    {
    }

    // Serialization for Pstream communication
    void write(Foam::Ostream& os) const
    {
        os << pCoords << " " << dz << " " << dx << " " << dy << " " << area
           << " " << isTopCentre;
    }

    // Deserialization for Pstream communication
    void read(Foam::Istream& is)
    {
        is >> pCoords >> dz >> dx >> dy >> area >> isTopCentre;
    }
};

// --- NEW: Global comparison operators for Pstream ---
inline bool operator==(const InterpolationSource& a,
                       const InterpolationSource& b)
{
    return (a.pCoords == b.pCoords && a.dz == b.dz && a.dx == b.dx
            && a.dy == b.dy && a.area == b.area
            && a.isTopCentre == b.isTopCentre);
}

inline bool operator!=(const InterpolationSource& a,
                       const InterpolationSource& b)
{
    return !(a == b);
}
// --- END NEW ---

// Overload of stream operators for easy serialization/deserialization
Foam::Ostream& operator<<(Foam::Ostream& os, const InterpolationSource& s)
{
    s.write(os);
    return os;
}

Foam::Istream& operator>>(Foam::Istream& is, InterpolationSource& s)
{
    s.read(is);
    return is;
}


void generateCroppedDEM(const RectangularMatrix<double>& elevation,
                        scalar xllcorner,
                        scalar yllcorner,
                        scalar cellsize,
                        scalar xVent,
                        scalar yVent, // Translation factors
                        scalar xmin,
                        scalar xmax,
                        scalar ymin,
                        scalar ymax,
                        const Foam::fileName& outputFileName)
{
    // Adjust the domain bounds to match the DEM coordinate system
    xmin += xVent;
    xmax += xVent;
    ymin += yVent;
    ymax += yVent;

    // Compute the new cell size ensuring an integer number of rows/cols
    int ncols_new = round((xmax - xmin) / cellsize);
    int nrows_new = round((ymax - ymin) / cellsize);
    scalar cellsize_new = (xmax - xmin) / ncols_new; // Ensure exact fit

    // Compute the new lower-left corner (adjusting for cell center convention)
    scalar xllcorner_new = xmin + 0.5 * cellsize_new;
    scalar yllcorner_new = ymin + 0.5 * cellsize_new;

    Info << "Generating cropped DEM with cellsize: " << cellsize_new << " ("
         << ncols_new << " x " << nrows_new << " grid)" << endl;

    // Open output file
    std::ofstream file(outputFileName);
    if (!file)
    {
        FatalErrorInFunction << "Cannot open output file " << outputFileName
                             << exit(FatalError);
    }

    // Write ASCII Raster Header
    file << "ncols " << ncols_new << "\n";
    file << "nrows " << nrows_new << "\n";
    file << "xllcorner " << xllcorner_new << "\n";
    file << "yllcorner " << yllcorner_new << "\n";
    file << "cellsize " << cellsize_new << "\n";
    file << "NODATA_value -9999\n";

    // Bilinear interpolation over the new grid
    for (int row = 0; row < nrows_new; ++row)
    {
        for (int col = 0; col < ncols_new; ++col)
        {
            // Compute world coordinates of the new cell center
            scalar x = xllcorner_new + (col + 0.5) * cellsize_new;
            scalar y = yllcorner_new
                     + (nrows_new - row - 0.5) * cellsize_new; // Top to bottom

            // Compute corresponding indices in the original DEM
            int i = (y - (yllcorner + 0.5 * cellsize)) / cellsize;
            int j = (x - (xllcorner + 0.5 * cellsize)) / cellsize;

            if (i >= 0 && i < elevation.m() - 1 && j >= 0
                && j < elevation.n() - 1)
            {
                // Compute interpolation weights
                scalar xLerp = (x - (xllcorner + j * cellsize + 0.5 * cellsize))
                             / cellsize;
                scalar yLerp = (y - (yllcorner + i * cellsize + 0.5 * cellsize))
                             / cellsize;

                // Get the four surrounding elevation values
                scalar v00 = elevation(i, j);
                scalar v01 = elevation(i, j + 1);
                scalar v10 = elevation(i + 1, j);
                scalar v11 = elevation(i + 1, j + 1);

                // Bilinear interpolation
                scalar zInterp =
                    v00 * (1 - xLerp) * (1 - yLerp) + v01 * xLerp * (1 - yLerp)
                    + v10 * (1 - xLerp) * yLerp + v11 * xLerp * yLerp;

                file << zInterp << " ";
            }
            else
            {
                file << "-9999 ";
            }
        }
        file << "\n";
    }

    file.close();
    Info << "Cropped DEM saved as " << outputFileName << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Function to compute the normal vector of a triangle formed by points
// p1, p2, p3
vector computeNormal(const point& p1, const point& p2, const point& p3)
{
    vector v1 = p2 - p1; // Edge vector 1
    vector v2 = p3 - p1; // Edge vector 2

    // Compute the cross product of v1 and v2 to get the normal
    vector normal =
        Foam::vector(v1.y() * v2.z() - v1.z() * v2.y(), // x-component
                     v1.z() * v2.x() - v1.x() * v2.z(), // y-component
                     v1.x() * v2.y() - v1.y() * v2.x()  // z-component
        );

    // Normalize the normal vector
    scalar magnitude = mag(normal);
    if (magnitude > SMALL)
    {
        normal /= magnitude;
    }

    return normal;
}

// Function to write a single triangle in binary format
void writeBinaryTriangle(std::ofstream& stlFile,
                         const vector& normal,
                         const point& p1,
                         const point& p2,
                         const point& p3)
{
    // Write normal vector (12 bytes: 3 floats)
    float nx(normal.x());
    float ny(normal.y());
    float nz(normal.z());
    stlFile.write(reinterpret_cast<const char*>(&nx), 4);
    stlFile.write(reinterpret_cast<const char*>(&ny), 4);
    stlFile.write(reinterpret_cast<const char*>(&nz), 4);

    // Write vertex 1 (12 bytes: 3 floats)
    float P1x(p1.x());
    float P1y(p1.y());
    float P1z(p1.z());
    stlFile.write(reinterpret_cast<const char*>(&P1x), 4);
    stlFile.write(reinterpret_cast<const char*>(&P1y), 4);
    stlFile.write(reinterpret_cast<const char*>(&P1z), 4);

    // Write vertex 2 (12 bytes: 3 floats)
    float P2x(p2.x());
    float P2y(p2.y());
    float P2z(p2.z());
    stlFile.write(reinterpret_cast<const char*>(&P2x), 4);
    stlFile.write(reinterpret_cast<const char*>(&P2y), 4);
    stlFile.write(reinterpret_cast<const char*>(&P2z), 4);

    // Write vertex 3 (12 bytes: 3 floats)
    float P3x(p3.x());
    float P3y(p3.y());
    float P3z(p3.z());
    stlFile.write(reinterpret_cast<const char*>(&P3x), 4);
    stlFile.write(reinterpret_cast<const char*>(&P3y), 4);
    stlFile.write(reinterpret_cast<const char*>(&P3z), 4);

    // Write attribute byte count (2 bytes, set to 0)
    char attribute[2] = "0";
    stlFile.write(attribute, 2);
}

// Function to write STL surface in binary format
void writeBinarySTL(const word& stlFileName,
                    const RectangularMatrix<scalar>& elevation,
                    scalar xOffset,
                    scalar yOffset,
                    scalar cellSize)
{
    std::ofstream stlFile(stlFileName.c_str(), std::ios::binary);
    if (!stlFile)
    {
        FatalErrorInFunction << "Cannot open STL file " << stlFileName
                             << " for writing" << exit(FatalError);
    }

    // Write 80-byte header (just fill with 0 or any message)
    char header[80] = "Generated by OpenFOAM";
    stlFile.write(header, sizeof(header));

    // Get dimensions of the elevation grid
    const label numRows = elevation.m();
    const label numCols = elevation.n();

    // Write triangle count (4 bytes)
    auto triangleCount =
        static_cast<uint32_t>(2 * (numRows - 1) * (numCols - 1));
    stlFile.write(reinterpret_cast<const char*>(&triangleCount),
                  sizeof(triangleCount));

    // Loop over each cell in the grid and create two triangles per cell
    for (label i = 0; i < numRows - 1; ++i)
    {
        for (label j = 0; j < numCols - 1; ++j)
        {
            // Get corner points of the cell

            // Top-left corner
            point p1(xOffset + j * cellSize,
                     yOffset + i * cellSize,
                     elevation(i, j));

            // Top-right corner
            point p2(xOffset + (j + 1) * cellSize,
                     yOffset + i * cellSize,
                     elevation(i, j + 1));

            // Bottom-left corner
            point p3(xOffset + j * cellSize,
                     yOffset + (i + 1) * cellSize,
                     elevation(i + 1, j));

            // Bottom-right corner
            point p4(xOffset + (j + 1) * cellSize,
                     yOffset + (i + 1) * cellSize,
                     elevation(i + 1, j + 1));

            // First triangle (p1, p2, p3)
            vector normal1 = computeNormal(p1, p2, p3);

            // Info << p1 << p2 << p3 << normal1 << endl;
            writeBinaryTriangle(stlFile, normal1, p1, p2, p3);

            // Second triangle (p2, p4, p3)
            vector normal2 = computeNormal(p2, p4, p3);
            writeBinaryTriangle(stlFile, normal2, p2, p4, p3);
        }
    }

    stlFile.close();
    Info << "Binary STL surface written to " << stlFileName << endl;
}

// Function to write STL surface from the elevation grid
void writeSTL(const word& stlFileName,
              const RectangularMatrix<scalar>& elevation,
              scalar xOffset,
              scalar yOffset,
              scalar cellSize)
{
    std::ofstream stlFile(stlFileName.c_str());
    if (!stlFile)
    {
        FatalErrorInFunction << "Cannot open STL file " << stlFileName
                             << " for writing" << exit(FatalError);
    }

    // Write STL file header
    stlFile << "solid topoSurface" << "\n";

    // Get dimensions of the elevation grid
    const label numRows = elevation.m();
    const label numCols = elevation.n();

    // Loop over each cell in the grid and create two triangles per cell
    for (label i = 0; i < numRows - 1; ++i)
    {
        for (label j = 0; j < numCols - 1; ++j)
        {
            // Get corner points of the cell (elevation grid)

            // Top-left corner
            vector p1(xOffset + j * cellSize,
                      yOffset + i * cellSize,
                      elevation(i, j));

            // Top-right corner
            vector p2(xOffset + (j + 1) * cellSize,
                      yOffset + i * cellSize,
                      elevation(i, j + 1));

            // Bottom-left corner
            vector p3(xOffset + j * cellSize,
                      yOffset + (i + 1) * cellSize,
                      elevation(i + 1, j));

            // Bottom-right corner
            vector p4(xOffset + (j + 1) * cellSize,
                      yOffset + (i + 1) * cellSize,
                      elevation(i + 1, j + 1));
            // First triangle (p1, p2, p3) - Top-left, top-right, bottom-left
            vector normal1 = computeNormal(p1, p2, p3);
            stlFile << "  facet normal " << normal1.x() << " " << normal1.y()
                    << " " << normal1.z() << "\n";
            stlFile << "    outer loop" << "\n";
            stlFile << "      vertex " << p1.x() << " " << p1.y() << " "
                    << p1.z() << "\n";
            stlFile << "      vertex " << p2.x() << " " << p2.y() << " "
                    << p2.z() << "\n";
            stlFile << "      vertex " << p3.x() << " " << p3.y() << " "
                    << p3.z() << "\n";
            stlFile << "    endloop" << "\n";
            stlFile << "  endfacet" << "\n";

            // Second triangle (p2, p4, p3) - Top-right, bottom-right,
            // bottom-left
            vector normal2 = computeNormal(p2, p4, p3);
            stlFile << "  facet normal " << normal2.x() << " " << normal2.y()
                    << " " << normal2.z() << "\n";
            stlFile << "    outer loop" << "\n";
            stlFile << "      vertex " << p2.x() << " " << p2.y() << " "
                    << p2.z() << "\n";
            stlFile << "      vertex " << p4.x() << " " << p4.y() << " "
                    << p4.z() << "\n";
            stlFile << "      vertex " << p3.x() << " " << p3.y() << " "
                    << p3.z() << "\n";
            stlFile << "    endloop" << "\n";
            stlFile << "  endfacet" << "\n";
        }
    }

    // Write STL file footer
    stlFile << "endsolid topoSurface" << "\n";

    stlFile.close();
    Info << "STL surface written to " << stlFileName << endl;
}

scalar minQuality(const polyMesh& mesh,
                  const point& cC,
                  const label fI,
                  const bool isOwner,
                  const label faceBasePtI)
{
    // Does fan decomposition of face (starting at faceBasePti) and determines
    // min quality over all resulting tets.

    const pointField& pPts = mesh.points();
    const face& f = mesh.faces()[fI];
    const point& tetBasePt = pPts[f[faceBasePtI]];

    scalar thisBaseMinTetQuality = vGreat;

    for (label tetPtI = 1; tetPtI < f.size() - 1; tetPtI++)
    {
        label facePtI = (tetPtI + faceBasePtI) % f.size();
        label otherFacePtI = f.fcIndex(facePtI);

        label ptAI = -1;
        label ptBI = -1;

        if (isOwner)
        {
            ptAI = f[facePtI];
            ptBI = f[otherFacePtI];
        }
        else
        {
            ptAI = f[otherFacePtI];
            ptBI = f[facePtI];
        }

        const point& pA = pPts[ptAI];
        const point& pB = pPts[ptBI];

        tetPointRef tet(cC, tetBasePt, pA, pB);

        scalar tetQuality = tet.quality();

        if (tetQuality < thisBaseMinTetQuality)
        {
            thisBaseMinTetQuality = tetQuality;
        }
    }
    return thisBaseMinTetQuality;
}

// Function to compute the face flatness
scalar calculateFlatness(const face& f, const pointField& points)
{
    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(f, fp)
    {
        pAvg += points[f[fp]];
    }
    pAvg /= f.size();

    // Compute the face area normal and unit normal
    vector sumA = Zero;
    forAll(f, fp)
    {
        const point& p = points[f[fp]];
        const point& pNext = points[f.nextLabel(fp)];

        const vector a = (pNext - p) ^ (pAvg - p);
        sumA += a;
    }
    const vector sumAHat = normalised(sumA);

    // Compute the area-weighted sum of the triangle centres
    scalar sumAn = 0;
    vector sumAnc = Zero;
    forAll(f, fp)
    {
        const point& p = points[f[fp]];
        const point& pNext = points[f.nextLabel(fp)];

        const vector a = (pNext - p) ^ (pAvg - p);
        const vector c = p + pNext + pAvg;

        const scalar an = a & sumAHat;

        sumAn += an;
        sumAnc += an * c;
    }
    point fc = (1.0 / 3.0) * sumAnc / sumAn;

    // Calculate the sum of the magnitude of areas and compare to magnitude
    // of sum of areas
    scalar summA = 0.0;
    vector sumN = Zero;

    forAll(f, fp)
    {
        const point& thisPoint = points[f[fp]];
        const point& nextPoint = points[f.nextLabel(fp)];

        // Triangle around fc
        const vector n = 0.5 * ((nextPoint - thisPoint) ^ (fc - thisPoint));

        summA += mag(n);
        sumN += n;
    }

    scalar magArea = mag(sumN);
    scalar faceFlatness = magArea / summA;

    return faceFlatness;
}

// Function to compute the interpolation of dz at any mesh point, based on the
// inverse of the distance
point inverseDistanceInterpolationDz(const scalar& Ldef,
                                     const scalar& alpha,
                                     const scalar& coeffVertDeformation,
                                     const point& internalPoint,
                                     const scalarField& boundaryPointsX,
                                     const scalarField& boundaryPointsY,
                                     const scalarField& boundaryPointsZ,
                                     const scalarField& boundaryDz,
                                     const scalarField& boundaryDx,
                                     const scalarField& boundaryDy,
                                     const scalarField& boundaryAreas,
                                     const scalar& distThr)
{
    // Initialize variables
    point DeltaInterp;
    const label n = boundaryDz.size();

    // Precompute alpha^5
    const scalar alpha5 = alpha * alpha * alpha * alpha * alpha;

    // Variables for interpolation
    scalar NumZ(0.0);
    scalar NumX(0.0);
    scalar NumY(0.0);
    scalar Den(0.0);
    scalar Den_z(0.0);
    scalar distance_z, LbyD_z, LbyD3_z, weight_z;
    scalar distance, LbyD, LbyD3, weight;
    scalar dist2_xy;
    scalar dist2_z;

    for (label i = 0; i < n; ++i)
    {

        dist2_xy = sqr(internalPoint.x() - boundaryPointsX[i])
                 + sqr(internalPoint.y() - boundaryPointsY[i]);

        dist2_z = sqr(internalPoint.z() - boundaryPointsZ[i]);

        distance = Foam::sqrt(dist2_xy + dist2_z);

        distance_z = Foam::sqrt(dist2_xy + coeffVertDeformation * dist2_z);

        if (distance < 1.e-3 * distThr)
        {
            DeltaInterp = vector(boundaryDx[i], boundaryDy[i], boundaryDz[i]);
            return DeltaInterp;
        }

        LbyD = Ldef / distance;
        LbyD3 = LbyD * LbyD * LbyD;
        weight = boundaryAreas[i] * (LbyD3 + alpha5 * LbyD3 * LbyD * LbyD);

        NumX += weight * boundaryDx[i];
        NumY += weight * boundaryDy[i];
        Den += weight;

        if (coeffVertDeformation < 1.0)
        {
            LbyD_z = Ldef / distance_z;
            LbyD3_z = LbyD_z * LbyD_z * LbyD_z;
            weight_z = boundaryAreas[i]
                     * (LbyD3_z + alpha5 * LbyD3_z * LbyD_z * LbyD_z);
        }
        else
        {
            weight_z = weight;
        }

        NumZ += weight_z * boundaryDz[i];
        Den_z += weight_z;
    }

    DeltaInterp = vector(NumX / Den, NumY / Den, NumZ / Den_z);

    return DeltaInterp;
}

// Function to compute the interpolation of dz at z=0 points, based on the
// inverse of the distance
Tuple2<scalar, scalar>
inverseDistanceInterpolationDzBottom(const point& internalPoint,
                                     const scalarField& boundaryPointsX,
                                     const scalarField& boundaryPointsY,
                                     const scalarField& boundaryVal1,
                                     const scalarField& boundaryVal2,
                                     const scalar& interpRelRadius,
                                     const scalar& distThr)
{
    scalar interpolatedVal1(0.0);
    scalar interpolatedVal2(0.0);

    const label n = boundaryVal1.size();
    scalar minValue = GREAT;
    label minIndex = -1;
    const scalar eps = 1e-2;

    // Calculate distances and find the minimum in a single loop
    scalarField distances(n);
    for (label i = 0; i < n; ++i)
    {
        distances[i] =
            Foam::sqrt(Foam::sqr(internalPoint.x() - boundaryPointsX[i])
                       + Foam::sqr(internalPoint.y() - boundaryPointsY[i]));

        if (distances[i] < minValue)
        {
            minValue = distances[i];
            minIndex = i;
        }
    }

    // Special case: very close to a boundary point
    if (minValue < 1.e-3 * distThr)
    {
        interpolatedVal1 = boundaryVal1[minIndex];
        interpolatedVal2 = boundaryVal2[minIndex];
    }
    else
    {
        // General case: weighted interpolation
        scalar NumVal2(0.0), NumVal1(0.0), Den(0.0);

        const scalar radiusThreshold = interpRelRadius * minValue;

        for (label i = 0; i < n; ++i)
        {
            scalar distance = distances[i];
            scalar weight = 0.0;

            // Neglect points outside the relative radius
            if ((distance / radiusThreshold - 1.0) < eps)
            {
                weight = minValue / (distance * distance);
            }

            NumVal1 += weight * boundaryVal1[i];
            NumVal2 += weight * boundaryVal2[i];
            Den += weight;
        }
        interpolatedVal1 = NumVal1 / Den;
        interpolatedVal2 = NumVal2 / Den;
    }

    return Tuple2<scalar, scalar>(interpolatedVal1, interpolatedVal2);
}

// Function to calculate the average of a sub-block
double calculateBlockAverage(const RectangularMatrix<double>& elevation,
                             label startRow,
                             label startCol,
                             label blockSize,
                             label maxRows,
                             label maxCols)
{
    double sum = 0.0;
    label count = 0;

    for (label i = startRow; i < startRow + blockSize && i < maxRows; ++i)
    {
        for (label j = startCol; j < startCol + blockSize && j < maxCols; ++j)
        {
            sum += elevation(i, j);
            ++count;
        }
    }
    return sum / count;
}

// Function to subsample a 2D array (used to sub-sample the
// elevation data and save the STL file)
RectangularMatrix<double>
subsampleMatrix(const RectangularMatrix<double>& elevation,
                int ncols,
                int nrows,
                label blockSize)
{
    // New dimensions
    label newRows = nrows / blockSize;
    label newCols = ncols / blockSize;

    // Create a new matrix for the subsampled data
    RectangularMatrix<double> subsampled(newRows, newCols, 0.0);

    // Fill the subsampled matrix
    for (label i = 0; i < newRows; ++i)
    {
        for (label j = 0; j < newCols; ++j)
        {
            subsampled(i, j) = calculateBlockAverage(elevation,
                                                     i * blockSize,
                                                     j * blockSize,
                                                     blockSize,
                                                     nrows,
                                                     ncols);
        }
    }
    return subsampled;
}

Tuple2<scalar, scalar> interpolateNegDeformation(scalar z,
                                                 const bool useNegDeformation,
                                                 const scalarList& zNeg,
                                                 const scalarList& dxNeg,
                                                 const scalarList& dyNeg)
{
    if (!useNegDeformation) return Tuple2<scalar, scalar>(0.0, 0.0);

    if (z >= 0.0)
        return Tuple2<scalar, scalar>(0.0,
                                      0.0); // Above or at z=0 → No deformation

    if (z > zNeg[0]) // Interpolate between (0,0) at z=0 and (dxNeg[0],
                     // dyNeg[0]) at zNeg[0]
    {
        scalar w = z / zNeg[0]; // Weight factor (z=0 → w=0, z=zNeg[0] → w=1)
        scalar interpDx = w * dxNeg[0];
        scalar interpDy = w * dyNeg[0];
        return Tuple2<scalar, scalar>(interpDx, interpDy);
    }

    if (z <= zNeg.last()) // Below lowest zNeg → Constant deformation
    {
        return Tuple2<scalar, scalar>(dxNeg.last(), dyNeg.last());
    }

    // Find the two closest points for linear interpolation
    for (label i = 0; i < zNeg.size() - 1; ++i)
    {
        if (zNeg[i] >= z && z > zNeg[i + 1])
        {
            scalar w = (z - zNeg[i + 1])
                     / (zNeg[i] - zNeg[i + 1]); // Interpolation weight
            scalar interpDx = w * dxNeg[i] + (1 - w) * dxNeg[i + 1];
            scalar interpDy = w * dyNeg[i] + (1 - w) * dyNeg[i + 1];
            return Tuple2<scalar, scalar>(interpDx, interpDy);
        }
    }

    return Tuple2<scalar, scalar>(0.0, 0.0); // Should never reach this point
}

//--------------------------------------------------------------

int main(int argc, char* argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"

    // Read the dictionary file (topoGridDict) from the "system" folder
    IOdictionary topoDict(IOobject("topoGridDict",
                                   runTime.system(),
                                   mesh,
                                   IOobject::MUST_READ,
                                   IOobject::NO_WRITE));

    // Declare variables that need to be accessible outside `if (deform)` block
    // or need to be initialized before `if (deform)`
    scalar Ldef = 0.0;
    scalar maxTopo = 0.0;
    scalar noDeformLevel = 0.0;
    scalar alphaAll = 0.0;
    scalar zMin = 0.0; // Declare zMin here
    scalar zMax = 0.0; // Declare zMax here

    // Final aggregated and sorted source fields for interpolation
    scalarField globalPointsX;
    scalarField globalPointsY;
    scalarField globalPointsZ;
    scalarField globalDz;
    scalarField globalDx;
    scalarField globalDy;
    scalarField globalAreas;


    // Read the switch to deform the mesh
    const Switch deform = topoDict.lookupOrDefault<Switch>("deform", false);

    if (deform)
    {
        // --- Original initial variable declarations ---
        const word rasterFile = topoDict.lookup<word>("rasterFile");
        const scalar xVent = topoDict.lookupOrDefault<scalar>("xVent", 0.0);
        const scalar yVent = topoDict.lookupOrDefault<scalar>("yVent", 0.0);
        const scalar interpRelRadius =
            topoDict.lookupOrDefault<scalar>("interpRelRadius", 4.0);
        const Switch saveSTL =
            topoDict.lookupOrDefault<Switch>("saveSTL", false);
        const Switch saveBinary =
            topoDict.lookupOrDefault<Switch>("saveBinary", false);
        const Switch checkMesh =
            topoDict.lookupOrDefault<Switch>("checkMesh", false);
        const Switch raiseTop =
            topoDict.lookupOrDefault<Switch>("raiseTop", true);
        const Switch orthogonalCorrection =
            topoDict.lookupOrDefault<Switch>("orthogonalCorrection", false);
        const scalar dist_rel1 =
            topoDict.lookupOrDefault<scalar>("dist_rel1", 0.1);
        const scalar dist_rel2 =
            topoDict.lookupOrDefault<scalar>("dist_rel2", 0.2);
        const scalar distC1 = topoDict.lookupOrDefault<scalar>("distC1", 0.0);
        const scalar distC2 = topoDict.lookupOrDefault<scalar>("distC2", 0.0);
        const scalar noDeformCoeff =
            topoDict.lookupOrDefault<scalar>("noDeformCoeff", 0.5);
        const scalar noDeformExp =
            topoDict.lookupOrDefault<scalar>("noDeformExp", 1.0);
        const Switch saveCrop =
            topoDict.lookupOrDefault<Switch>("saveCrop", false);
        const scalar coeffVertDeformation =
            topoDict.lookupOrDefault<scalar>("coeffVertDeformation", 1.0);

        // Optional flattening of the vertical displacement below z=0.
        // Near z=0 the local interpolated displacement is preserved; toward
        // the original global minimum mesh z it tends to the topographic
        // elevation at (flattenBottomX, flattenBottomY).
        const Switch flattenBottom =
            topoDict.lookupOrDefault<Switch>("flattenBottom", false);
        const scalar flattenBottomX =
            topoDict.lookupOrDefault<scalar>("flattenBottomX", 0.0);
        const scalar flattenBottomY =
            topoDict.lookupOrDefault<scalar>("flattenBottomY", 0.0);
        const scalar flattenBottomExp =
            topoDict.lookupOrDefault<scalar>("flattenBottomExp", 2.0);

        if (flattenBottom && flattenBottomExp <= SMALL)
        {
            FatalErrorInFunction << "flattenBottomExp must be greater than zero"
                                 << exit(FatalError);
        }

        scalar dzFlat = 0.0;

        scalarList zNeg, dxNeg, dyNeg;
        bool useNegDeformation = true;
        // --- End original declarations ---

        if (topoDict.found("zNeg") && topoDict.found("dxNeg")
            && topoDict.found("dyNeg"))
        {
            zNeg = topoDict.lookup<scalarList>("zNeg");
            dxNeg = topoDict.lookup<scalarList>("dxNeg");
            dyNeg = topoDict.lookup<scalarList>("dyNeg");

            // Ensure zNeg is sorted in decreasing order
            for (label i = 0; i < zNeg.size() - 1; ++i)
            {
                if (zNeg[i] < zNeg[i + 1])
                {
                    FatalErrorInFunction
                        << "zNeg list must be sorted in "
                           "decreasing order (less negative first)"
                        << exit(FatalError);
                }
            }

            if (zNeg.size() != dxNeg.size() || zNeg.size() != dyNeg.size())
            {
                FatalErrorInFunction
                    << "zNeg, dxNeg, and dyNeg must have the same size"
                    << exit(FatalError);
            }

            Info << "Read " << zNeg.size() << " negative deformation levels."
                 << endl;
        }
        else
        {
            Info << "Missing zNeg, dxNeg, or dyNeg in topoGridDict. Horizontal "
                    "deformation will be set to zero."
                 << endl;
            useNegDeformation = false;
        }

        Foam::fileName pathPrefix = "./constant/DEM/";
        Foam::fileName fullRasterFilePath = pathPrefix / rasterFile;

        Info << "Raster file specified: " << fullRasterFilePath << endl;

        std::ifstream file(fullRasterFilePath);
        if (!file.is_open())
        {
            FatalErrorInFunction
                << "Unable to open the raster file: " << fullRasterFilePath
                << exit(FatalError);
        }

        int ncols = 0, nrows = 0;
        double xllcorner = 0.0, yllcorner = 0.0, cellsize = 0.0;
        double NODATA_value = -9999.0;
        std::string line;

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            std::string key;
            iss >> key;

            if (key == "ncols" || key == "NCOLS") iss >> ncols;
            else if (key == "nrows" || key == "NROWS")
                iss >> nrows;
            else if (key == "xllcorner" || key == "XLLCORNER")
                iss >> xllcorner;
            else if (key == "yllcorner" || key == "YLLCORNER")
                iss >> yllcorner;
            else if (key == "cellsize" || key == "CELLSIZE")
                iss >> cellsize;
            else if (key == "NODATA_value" || key == "NODATA_VALUE")
                iss >> NODATA_value;
            if (key == "NODATA_value" || key == "NODATA_VALUE") break;
        }

        xllcorner -= xVent;
        yllcorner -= yVent;

        RectangularMatrix<double> elevation(nrows, ncols, 0.0);
        for (int i = 0; i < nrows; ++i)
        {
            std::getline(file, line);
            std::istringstream iss(line);
            for (int j = 0; j < ncols; ++j)
            {
                double value;
                iss >> value;
                if (value == NODATA_value) value = 0.0;
                elevation(nrows - 1 - i, j) = value;
            }
        }

        if (saveSTL)
        {
            if (Pstream::master())
            {
                label factor = 2;
                RectangularMatrix<double> elevationSubsampled =
                    subsampleMatrix(elevation, ncols, nrows, factor);
                scalar xllSubsampled(xllcorner + (0.5 * factor) * cellsize);
                scalar yllSubsampled(yllcorner + (0.5 * factor) * cellsize);
                scalar cellsizeSubsampled(factor * cellsize);
                word stlFileName(fullRasterFilePath);
                stlFileName.replace(".asc", ".stl");
                Info << "Saving STL file: " << stlFileName << endl;
                if (saveBinary)
                {
                    writeBinarySTL(stlFileName,
                                   elevationSubsampled,
                                   xllSubsampled,
                                   yllSubsampled,
                                   cellsizeSubsampled);
                }
                else
                {
                    writeSTL(stlFileName,
                             elevationSubsampled,
                             xllSubsampled,
                             yllSubsampled,
                             cellsizeSubsampled);
                }
                Info << "Saving completed" << endl;
            }
            label dummy = 0;
            reduce(dummy, sumOp<label>());
        }
        file.close();

        if (flattenBottom)
        {
            // Reference flat displacement: topographic elevation at the
            // user-selected point. Since the bottom reference surface is z=0,
            // this is also the vertical displacement to use as the deep,
            // laterally uniform limit.
            const scalar x_grid = (flattenBottomX - xllcorner) / cellsize - 0.5;
            const scalar y_grid = (flattenBottomY - yllcorner) / cellsize - 0.5;

            const label colIndex = static_cast<label>(floor(x_grid));
            const label rowIndex = static_cast<label>(floor(y_grid));

            if (colIndex >= 0 && colIndex < ncols - 1 && rowIndex >= 0
                && rowIndex < nrows - 1)
            {
                const scalar xLerp = x_grid - colIndex;
                const scalar yLerp = y_grid - rowIndex;

                const scalar v00 = elevation(rowIndex, colIndex);
                const scalar v01 = elevation(rowIndex, colIndex + 1);
                const scalar v10 = elevation(rowIndex + 1, colIndex);
                const scalar v11 = elevation(rowIndex + 1, colIndex + 1);

                dzFlat = v00 * (1 - xLerp) * (1 - yLerp)
                       + v01 * xLerp * (1 - yLerp) + v10 * (1 - xLerp) * yLerp
                       + v11 * xLerp * yLerp;

                Info << "flattenBottom enabled" << endl;
                Info << "flattenBottomX = " << flattenBottomX << endl;
                Info << "flattenBottomY = " << flattenBottomY << endl;
                Info << "flattenBottomExp = " << flattenBottomExp << endl;
                Info << "Deep flat vertical displacement dzFlat = " << dzFlat
                     << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "flattenBottom point (" << flattenBottomX << ", "
                    << flattenBottomY << ") out of DEM bounds "
                    << "x_grid: " << x_grid << ", y_grid: " << y_grid
                    << ", colIndex: " << colIndex << ", rowIndex: " << rowIndex
                    << ", ncols: " << ncols << ", nrows: " << nrows
                    << exit(FatalError);
            }
        }

        scalar xMin = min(mesh.Cf().component(0)).value();
        scalar xMax = max(mesh.Cf().component(0)).value();
        reduce(xMin, minOp<scalar>());
        reduce(xMax, maxOp<scalar>());
        Info << "xMin = " << xMin << endl;
        Info << "xMax = " << xMax << endl;

        scalar yMin = min(mesh.Cf().component(1)).value();
        scalar yMax = max(mesh.Cf().component(1)).value();
        reduce(yMin, minOp<scalar>());
        reduce(yMax, maxOp<scalar>());
        Info << "yMin = " << yMin << endl;
        Info << "yMax = " << yMax << endl;

        // zMin/zMax from mesh are needed for Ldef calculation below
        zMin = min(mesh.Cf().component(2))
                   .value(); // Access zMin variable declared outside if block
        zMax = max(mesh.Cf().component(2))
                   .value(); // Access zMax variable declared outside if block
        reduce(zMin, minOp<scalar>());
        reduce(zMax, maxOp<scalar>());
        Info << "zMin = " << zMin << endl;
        Info << "zMax = " << zMax << endl;

        label patchID = -1;
        forAll(mesh.boundaryMesh(), patchi)
        {
            if (mesh.boundaryMesh()[patchi].name() == "top")
            {
                patchID = patchi;
                break;
            }
        }
        if (patchID == -1)
        {
            FatalErrorInFunction << "Patch 'top' not found in mesh. "
                                 << "Cannot apply top surface deformation."
                                 << exit(FatalError);
        }

        const fvPatch& patchTop = mesh.boundary()[patchID];
        Info << "Preparing local top face centers for global interpolation "
                "sources..."
             << endl;


        if (saveCrop)
        {
            if (Pstream::master())
            {
                Foam::fileName croppedDEMFile = pathPrefix / "DEMcropped.asc";
                generateCroppedDEM(elevation,
                                   xllcorner + xVent,
                                   yllcorner + yVent,
                                   cellsize,
                                   xVent,
                                   yVent,
                                   xMin,
                                   xMax,
                                   yMin,
                                   yMax,
                                   croppedDEMFile);
            }
            label dummy = 0;
            reduce(dummy, sumOp<label>());
        }

        Ldef = (0.5 * std::sqrt(sqr(xMax - xMin) + sqr(yMax - yMin)));
        Info << "Ldef = " << Ldef << endl;

        noDeformLevel = (noDeformCoeff * Ldef);
        Info << "noDeformCoeff = " << noDeformCoeff << endl;
        Info << "noDeformLevel = " << noDeformLevel << endl << endl;

        const vectorField& faceAreas = mesh.faceAreas();
        const vectorField& faceCentres = mesh.faceCentres();
        const scalarField magFaceAreas(mag(faceAreas));
        const faceList& faces = mesh.faces();
        const pointField& points = mesh.points();

        scalar zFlatten = great;
        forAll(points, pointi)
        {
            zFlatten = min(zFlatten, points[pointi].z());
        }
        reduce(zFlatten, minOp<scalar>());
        Info << "zFlatten = global minimum original mesh point z = " << zFlatten
             << endl;

        if (flattenBottom && zFlatten >= -SMALL)
        {
            FatalErrorInFunction
                << "flattenBottom requires at least one mesh point below z=0. "
                << "The global minimum mesh point z is " << zFlatten
                << exit(FatalError);
        }

        scalar minLenSqr = sqr(great);
        scalar maxLenSqr = -sqr(great);
        labelHashSet smallEdgeSet(mesh.nPoints() / 100);

        forAll(faces, facei)
        {
            const face& f = faces[facei];
            forAll(f, fp)
            {
                label fp1 = f.fcIndex(fp);
                scalar magSqrE = magSqr(points[f[fp]] - points[f[fp1]]);
                minLenSqr = min(minLenSqr, magSqrE);
                maxLenSqr = max(maxLenSqr, magSqrE);
            }
        }
        reduce(minLenSqr, minOp<scalar>());
        reduce(maxLenSqr, maxOp<scalar>());
        Info << "Min/max edge length = " << Foam::sqrt(minLenSqr) << " "
             << Foam::sqrt(maxLenSqr) << endl;

        scalar distThr = Foam::sqrt(minLenSqr);


        // --- SECTION 1: PREPARATION OF Z=0 FACE CENTERS (FOR INITIAL
        // INTERPOLATION TO POINTS) --- This block largely retains the original
        // logic for processing z=0 faces and aggregating their data.

        labelList z0FaceIndices;
        forAll(faces, faceI)
        {
            if (mag(faceCentres[faceI].z()) < 1e-2 * distThr)
            {
                z0FaceIndices.append(faceI);
            }
        }
        Sout << "Proc" << Pstream::myProcNo() << " z=0 faces (local) "
             << z0FaceIndices.size() << endl;

        scalarField bottomCentresX_local(z0FaceIndices.size());
        scalarField bottomCentresY_local(z0FaceIndices.size());
        scalarField bottomCentresZ_local(z0FaceIndices.size());
        scalarField bottomAreas_local(z0FaceIndices.size());
        scalarField bottomCentresDz_local(z0FaceIndices.size());

        forAll(z0FaceIndices, facei)
        {
            point pCentre = faceCentres[z0FaceIndices[facei]];

            // Get x, y coordinates of the point
            scalar x = pCentre.x();
            scalar y = pCentre.y();

            // Transform the point's coordinates into a grid coordinate system
            // where integer values correspond to cell centers.
            scalar x_grid = (x - xllcorner) / cellsize - 0.5;
            scalar y_grid = (y - yllcorner) / cellsize - 0.5;

            // The index of the bottom-left cell for interpolation is the
            // integer part.
            int colIndex = static_cast<int>(floor(x_grid));
            int rowIndex = static_cast<int>(floor(y_grid));


            // Check that the indices are within the matrix bounds for
            // interpolation.
            if (colIndex >= 0 && colIndex < ncols - 1 && rowIndex >= 0
                && rowIndex < nrows - 1)
            {
                // The fractional part is the weight for the interpolation.
                scalar xLerp = x_grid - colIndex;
                scalar yLerp = y_grid - rowIndex;

                // Get the four surrounding elevation values.
                scalar v00 = elevation(rowIndex, colIndex);
                scalar v01 = elevation(rowIndex, colIndex + 1);
                scalar v10 = elevation(rowIndex + 1, colIndex);
                scalar v11 = elevation(rowIndex + 1, colIndex + 1);

                // Bilinear interpolation
                scalar zInterp =
                    v00 * (1 - xLerp) * (1 - yLerp) + v01 * xLerp * (1 - yLerp)
                    + v10 * (1 - xLerp) * yLerp + v11 * xLerp * yLerp;

                bottomCentresX_local[facei] =
                    faceCentres[z0FaceIndices[facei]].x();
                bottomCentresY_local[facei] =
                    faceCentres[z0FaceIndices[facei]].y();
                bottomCentresZ_local[facei] =
                    faceCentres[z0FaceIndices[facei]].z();
                bottomAreas_local[facei] = magFaceAreas[z0FaceIndices[facei]];
                bottomCentresDz_local[facei] =
                    zInterp - faceCentres[z0FaceIndices[facei]].z();
            }
            else
            {
                FatalErrorInFunction
                    << "Coordinates (" << x << ", " << y
                    << ") out of DEM bounds " << "x_grid: " << x_grid
                    << ", y_grid: " << y_grid << ", "
                    << "colIndex: " << colIndex << ", rowIndex: " << rowIndex
                    << ", " << "ncols: " << ncols << ", nrows: " << nrows
                    << exit(FatalError);
            }
        }

        // Aggregate Z=0 face center data from all processors for initial point
        // interpolation

        // Create list with nproc fields for CentresX
        List<scalarField> allProcBottomCentresX(Pstream::nProcs());
        allProcBottomCentresX[Pstream::myProcNo()] = bottomCentresX_local;
        Pstream::gatherList<scalarField>(allProcBottomCentresX);
        Pstream::scatterList<scalarField>(allProcBottomCentresX);

        // Create list with nproc fields for CentresY
        List<scalarField> allProcBottomCentresY(Pstream::nProcs());
        allProcBottomCentresY[Pstream::myProcNo()] = bottomCentresY_local;
        Pstream::gatherList<scalarField>(allProcBottomCentresY);
        Pstream::scatterList<scalarField>(allProcBottomCentresY);

        // Create list with nproc fields for CentresZ
        List<scalarField> allProcBottomCentresZ(Pstream::nProcs());
        allProcBottomCentresZ[Pstream::myProcNo()] = bottomCentresZ_local;
        Pstream::gatherList<scalarField>(allProcBottomCentresZ);
        Pstream::scatterList<scalarField>(allProcBottomCentresZ);

        // Create list with nproc fields for CentresDz
        List<scalarField> allProcBottomCentresDz(Pstream::nProcs());
        allProcBottomCentresDz[Pstream::myProcNo()] = bottomCentresDz_local;
        Pstream::gatherList<scalarField>(allProcBottomCentresDz);
        Pstream::scatterList<scalarField>(allProcBottomCentresDz);

        // Create list with nproc fields for CentresAreas
        List<scalarField> allProcBottomAreas(Pstream::nProcs());
        allProcBottomAreas[Pstream::myProcNo()] = bottomAreas_local;
        Pstream::gatherList<scalarField>(allProcBottomAreas);
        Pstream::scatterList<scalarField>(allProcBottomAreas);

        // Create the global fields
        scalarField globalBottomCentresX_agg;
        scalarField globalBottomCentresY_agg;
        scalarField globalBottomCentresZ_agg;
        scalarField globalBottomCentresDz_agg;
        scalarField globalBottomCentresAreas_agg;

        for (label i = 0; i < Pstream::nProcs(); ++i)
        {
            globalBottomCentresX_agg.append(allProcBottomCentresX[i]);
            globalBottomCentresY_agg.append(allProcBottomCentresY[i]);
            globalBottomCentresZ_agg.append(allProcBottomCentresZ[i]);
            globalBottomCentresDz_agg.append(allProcBottomCentresDz[i]);
            globalBottomCentresAreas_agg.append(allProcBottomAreas[i]);
        }

        Info << "Aggregated z=0 face centers (with duplicates if any): "
             << globalBottomCentresAreas_agg.size() << endl;

        // --- SECTION 2: ORTHOGONAL CORRECTION DX/DY FOR Z=0 FACES (PRESERVED)
        // --- This block is the original logic for calculating Dx/Dy for z=0
        // faces and aggregating it. Its results (globalBottomCentresDx_agg,
        // globalBottomCentresDy_agg) will be used to interpolate dx/dy for z=0
        // mesh points in the next step.


        // Temporary pointField to store the Z-deformed positions of original
        // mesh points at z=0. This is crucial for calculating face normals in
        // the orthogonal correction. Initialize with original mesh points.
        pointField z0Points_deformedZ(mesh.points());
        scalarField z0Points_areas(mesh.points().size());

        // STEP 2a: Deform Z for z=0 mesh points only (preliminary step)
        // This loop calculates only the vertical deformation (Dz) for points at
        // z=0 and updates their Z-coordinate in z0Points_deformedZ.
        Info << "Performing preliminary Z-deformation for local z=0 mesh "
                "points..."
             << endl;

        label localNumZ0Points = 0;

        forAll(mesh.points(), pointi)
        {
            // Original mesh point
            point pEval = mesh.points()[pointi];

            z0Points_deformedZ[pointi].z() = 0.0;
            z0Points_areas[pointi] = 0.0;

            if (mag(pEval.z()) < 1e-2 * distThr)
            {
                Tuple2<scalar, scalar> result =
                    inverseDistanceInterpolationDzBottom(
                        pEval,
                        globalBottomCentresX_agg,
                        globalBottomCentresY_agg,
                        globalBottomCentresDz_agg,
                        globalBottomCentresAreas_agg,
                        interpRelRadius,
                        distThr);

                z0Points_deformedZ[pointi].z() = result.first();
                z0Points_areas[pointi] = result.second();
                localNumZ0Points++;
            }
        }

        Info << "sum(globalDzFace) " << sum(globalBottomCentresDz_agg) << endl;
        Info << "sum(globalAreasFace) " << sum(globalBottomCentresAreas_agg)
             << endl;


        // STEP 2b: Calculate Orthogonal Correction (Dx/Dy) for Z=0 Faces using
        // Z-deformed points
        // We generate the dx/dy for FACES based on the *deformed*
        // z0Points. These will be aggregated and then interpolated to z=0 MESH
        // POINTS.

        scalarField dxFaceNormals_local(z0FaceIndices.size());
        scalarField dyFaceNormals_local(z0FaceIndices.size());

        if (orthogonalCorrection)
        {
            Info << "Calculating orthogonal correction Dx/Dy for z=0 faces "
                    "using Z-deformed points..."
                 << endl;
            // Calculate face normals using the z0Points_deformedZ for
            // connectivity
            forAll(z0FaceIndices, facei)
            {
                const face& f = faces[z0FaceIndices[facei]];

                // Using the calculateFlatness logic with
                // z0Points_deformedZ
                point pAvg = Zero;
                forAll(f, fp)
                {
                    pAvg += z0Points_deformedZ[f[fp]];
                } // Use deformed Z points
                pAvg /= f.size();

                vector sumA = Zero;
                forAll(f, fp)
                {
                    const point& p = z0Points_deformedZ[f[fp]];
                    const point& pNext = z0Points_deformedZ[f.nextLabel(fp)];
                    const vector a = (pNext - p) ^ (pAvg - p);
                    sumA += a;
                }
                const vector sumAHat = normalised(sumA);

                if (sumAHat.z() < 0.0)
                {
                    dxFaceNormals_local[facei] = -sumAHat.x();
                    dyFaceNormals_local[facei] = -sumAHat.y();
                }
                else
                {
                    dxFaceNormals_local[facei] = sumAHat.x();
                    dyFaceNormals_local[facei] = sumAHat.y();
                }
            }
        }
        else
        {
            dxFaceNormals_local = 0.0;
            dyFaceNormals_local = 0.0;
        }

        // Aggregate Dx interpolated from face normals
        List<scalarField> allProcFaceDxNormals(Pstream::nProcs());
        allProcFaceDxNormals[Pstream::myProcNo()] = dxFaceNormals_local;
        Pstream::gatherList<scalarField>(allProcFaceDxNormals);
        Pstream::scatterList<scalarField>(allProcFaceDxNormals);

        // Aggregate Dy interpolated from face normals
        List<scalarField> allProcFaceDyNormals(Pstream::nProcs());
        allProcFaceDyNormals[Pstream::myProcNo()] = dyFaceNormals_local;
        Pstream::gatherList<scalarField>(allProcFaceDyNormals);
        Pstream::scatterList<scalarField>(allProcFaceDyNormals);

        scalarField globalFaceCentresDx_agg;
        scalarField globalFaceCentresDy_agg;
        for (label i = 0; i < Pstream::nProcs(); ++i)
        {
            globalFaceCentresDx_agg.append(allProcFaceDxNormals[i]);
            globalFaceCentresDy_agg.append(allProcFaceDyNormals[i]);
        }
        Info << "Aggregated z=0 face Dx/Dy (from deformed Z-points) count: "
             << globalFaceCentresDx_agg.size() << endl;

        // --- Calcola i fattori di normalizzazione per i mesh points ---
        const labelList nProcsPerMeshPoint = countProcessorsPerMeshPoint(mesh);


        // --- SECTION 3: GLOBAL INTERPOLATION SOURCE PREPARATION (REPRODUCIBLE)
        // --- This is the core new logic to build a reproducible
        // list of sources for the final mesh point deformation.

        // Local list of source points for this processor
        label localNumTopCenters = patchTop.size();
        List<InterpolationSource> rawSources(localNumZ0Points
                                             + localNumTopCenters);

        // Index to fill rawSources sequentially
        label currentRawSourcesIdx = 0;

        // STEP 1: Populate rawSources with MESH POINTS at z=0 (after initial
        // Z-deformation)
        Info << "Populating rawSources with z=0 mesh points (final Dz, Dx, "
                "Dy)..."
             << endl;

        // Iterate over all LOCAL mesh points of this processor
        forAll(mesh.points(), pointi)
        {
            // Original mesh point (X,Y,Z original)
            point pEval_original = mesh.points()[pointi];

            if (mag(pEval_original.z()) < 1e-2 * distThr)
            {
                scalar dzValue = z0Points_deformedZ[pointi].z();
                scalar areaValue = z0Points_areas[pointi];

                scalar interpDx_z0 = 0.0;
                scalar interpDy_z0 = 0.0;
                if (orthogonalCorrection)
                {
                    Tuple2<scalar, scalar> dxdyResult =
                        inverseDistanceInterpolationDzBottom(
                            pEval_original,
                            globalBottomCentresX_agg,
                            globalBottomCentresY_agg,
                            globalFaceCentresDx_agg,
                            globalFaceCentresDy_agg,
                            interpRelRadius,
                            distThr);
                    interpDx_z0 = dxdyResult.first();
                    interpDy_z0 = dxdyResult.second();
                }

                const scalar normFactor =
                    1.0 / scalar(nProcsPerMeshPoint[pointi]);

                InterpolationSource s;
                // X,Y,Z original for the source point
                s.pCoords = pEval_original;

                // The deformation value (delta Z)
                s.dz = dzValue;

                // Derived dx for z=0 points
                s.dx = interpDx_z0;

                // Derived dy for z=0 points
                s.dy = interpDy_z0;

                s.area = areaValue * normFactor;
                s.isTopCentre = false;

                rawSources[currentRawSourcesIdx++] = s;
            }
        }

        // Sout << "Proc" << Pstream::myProcNo()
        //      << "Finished populating rawSources with z=0 mesh points. Count:
        //      "
        //      << rawSources.size() << endl;

        const label totalZ0PointsBeforeGlobalDeduplication =
            returnReduce(localNumZ0Points, sumOp<label>());

        Info << "Total z=0 points before global de-duplication (sum of local "
                "counts): "
             << totalZ0PointsBeforeGlobalDeduplication << endl;


        // PASSO 2: Populate rawSources with TOP FACE CENTERS
        // Retrieve the patchID - ORIGINAL CODE

        // maxTopo must be defined here after globalBottomCentresDz_agg is
        // populated
        if (raiseTop)
        {
            maxTopo = max(globalBottomCentresDz_agg);
        }
        else
        {
            maxTopo = 0.0;
        }

        // Iterate over local top faces of this processor
        forAll(patchTop, facei)
        {
            InterpolationSource s;
            s.pCoords =
                faceCentres[patchTop.start()
                            + facei]; // Original face center coordinates
            s.dz = maxTopo;
            s.dx = 0.0; // Top points have 0 horizontal deformation
            s.dy = 0.0;
            s.area = magFaceAreas[patchTop.start() + facei];
            s.isTopCentre = true;
            rawSources[currentRawSourcesIdx++] = s;
        }
        // Corrected calculation for Info message, manual count needed
        label numTopCenters = 0;
        for (label i = 0; i < rawSources.size(); ++i)
        {
            if (rawSources[i].isTopCentre) numTopCenters++;
        }
        // Sout << "Proc" << Pstream::myProcNo()
        //      << "Finished preparing local top face centers. Count: "
        //      << numTopCenters << endl;

        const label totalTopPointsBeforeGlobalDeduplication =
            returnReduce(numTopCenters, sumOp<label>());

        Info << "Total top points before global de-duplication (sum of local "
                "counts): "
             << totalTopPointsBeforeGlobalDeduplication << endl;


        // --- END OF LOCAL INTERPOLATION SOURCE PREPARATION ---


        // --- START OF GLOBAL AGGREGATION AND CANONICAL SORTING (REPRODUCIBLE)
        // ---

        // PASSO 3: Aggregate all rawSources from all processors
        Info << "Aggregating local interpolation sources from all processors..."
             << endl;
        List<List<InterpolationSource>> allProcRawSources(Pstream::nProcs());
        allProcRawSources[Pstream::myProcNo()] = rawSources;

        Pstream::gatherList<List<InterpolationSource>>(allProcRawSources);
        Pstream::scatterList<List<InterpolationSource>>(allProcRawSources);

        label actualTotalAggregatedSources = 0;
        forAll(allProcRawSources, procI)
        {
            actualTotalAggregatedSources += allProcRawSources[procI].size();
        }
        Info << "Aggregation complete. Total number of sources before GLOBAL "
                "de-duplication: "
             << actualTotalAggregatedSources << endl;


        // --- NUOVO: Semplificazione della combinazione delle sorgenti senza
        // pre-de-duplicazione --- Ora, semplicemente copiamo tutti gli elementi
        // aggregati in una singola lista La de-duplicazione efficace avverrà
        // DOPO l'ordinamento (nel PASSO 4a).
        List<InterpolationSource> finalGlobalSources(
            actualTotalAggregatedSources); // Pre-allocata alla dimensione
                                           // totale

        label currentFinalIdx = 0;
        forAll(allProcRawSources, procI)
        {
            forAll(allProcRawSources[procI], srcI)
            {
                finalGlobalSources[currentFinalIdx++] =
                    allProcRawSources[procI][srcI];
            }
        }

        Info << "Combined all aggregated sources into a single list. Total "
                "count: "
             << finalGlobalSources.size() << endl;


        // PASSO 4: Sort the combined list of sources in a canonical way
        // Criterio: isTopCentre (false prima di true), poi X, poi Y.
        Info << "Sorting global interpolation sources canonically "
                "(isTopCentre, then X, then Y)..."
             << endl;
        std::sort(finalGlobalSources.begin(),
                  finalGlobalSources.end(),
                  [](const InterpolationSource& a, const InterpolationSource& b)
                  {
                      if (a.isTopCentre != b.isTopCentre)
                          return a.isTopCentre < b.isTopCentre;
                      if (a.pCoords.x() != b.pCoords.x())
                          return a.pCoords.x() < b.pCoords.x();
                      return a.pCoords.y() < b.pCoords.y();
                  });
        Info << "Sorting complete." << endl;

        // PASSO 5: Popula le scalarField globali finali dalle sorgenti
        // normalizzate
        globalPointsX.setSize(finalGlobalSources.size());
        globalPointsY.setSize(finalGlobalSources.size());
        globalPointsZ.setSize(finalGlobalSources.size());
        globalDz.setSize(finalGlobalSources.size());
        globalDx.setSize(finalGlobalSources.size());
        globalDy.setSize(finalGlobalSources.size());
        globalAreas.setSize(finalGlobalSources.size());

        forAll(finalGlobalSources, i)
        {
            globalPointsX[i] = finalGlobalSources[i].pCoords.x();
            globalPointsY[i] = finalGlobalSources[i].pCoords.y();
            globalPointsZ[i] = finalGlobalSources[i].pCoords.z();
            globalDz[i] = finalGlobalSources[i].dz;
            globalDx[i] = finalGlobalSources[i].dx;
            globalDy[i] = finalGlobalSources[i].dy;
            globalAreas[i] = finalGlobalSources[i].area;
        }

        Info << "Final global points for deformation (with fractional "
                "normalization): "
             << finalGlobalSources.size() << endl;

        // --- ADAPTATION OF pDeform CALCULATION LOOP ---
        pointField pDeform(0.0 * mesh.points());

        const label totalPoints = mesh.points().size();
        label maxTotalPoints = totalPoints;
        reduce(maxTotalPoints, maxOp<label>());

        label localCount = 0;
        scalar nextPctg(1.0);

        scalar dxMin_rel;
        scalar dxMax_rel;
        scalar dyMin_rel;
        scalar dyMax_rel;
        scalar xCoeff;
        scalar yCoeff;
        scalar distC;
        scalar distCoeff;
        scalar coeffHor;
        Tuple2<scalar, scalar> result;

        // Calculate alphaAll here, after globalDz/globalAreas are ready
        scalar gamma = 5.0;
        Info << "sum(globalAreas) " << sum(globalAreas) << endl;
        Info << "sum(globalDz) " << sum(globalDz) << endl;

        scalar dzMean(sum(globalAreas / sum(globalAreas) * globalDz));
        alphaAll = gamma / Ldef * max(mag(globalDz - dzMean));
        Info << "alpha " << alphaAll << endl;


        // Loop over all points in the mesh to interpolate deformation
        Info
            << "Starting point deformation interpolation for all mesh points..."
            << endl;
        forAll(pDeform, pointi)
        {
            localCount++;

            // The percentage is computed with respect to the maximum number
            // of points among all the processors
            scalar GlobalPercentage =
                100.0 * static_cast<scalar>(localCount) / maxTotalPoints;

            if (GlobalPercentage >= nextPctg)
            {
                Info << "Progress: " << nextPctg << "% completed." << endl;
                nextPctg += 1.0;
            }

            // Use original mesh point for evaluation
            point pEval = mesh.points()[pointi];

            scalar interpDz(0.0);
            scalar interpDx(0.0);
            scalar interpDy(0.0);
            point DeltaInterp;

            if (mag(pEval.z() - zMax) < 1e-2 * distThr)
            {
                interpDz = maxTopo;
                interpDx = 0.0;
                interpDy = 0.0;
            }
            else
            {
                if (pEval.z() < 1e-2 * distThr)
                {
                    Tuple2<scalar, scalar> negDeform =
                        interpolateNegDeformation(
                            pEval.z(), useNegDeformation, zNeg, dxNeg, dyNeg);
                    interpDx = negDeform.first();
                    interpDy = negDeform.second();

                    point pEval_2D = pEval;
                    // For points on or below the topography consider only (x,y)
                    pEval_2D.z() = 0.0;

                    result = inverseDistanceInterpolationDzBottom(
                        pEval_2D,
                        globalBottomCentresX_agg,
                        globalBottomCentresY_agg,
                        globalBottomCentresDz_agg,
                        globalBottomCentresAreas_agg,
                        interpRelRadius,
                        distThr);
                    scalar dzLocal = result.first();

                    if (flattenBottom && pEval.z() < 0.0)
                    {
                        scalar s = pEval.z() / zFlatten;
                        s = min(1.0, max(0.0, s));

                        const scalar flattenWeight =
                            Foam::pow(1.0 - s, flattenBottomExp);

                        interpDz = flattenWeight * dzLocal
                                 + (1.0 - flattenWeight) * dzFlat;
                    }
                    else
                    {
                        interpDz = dzLocal;
                    }

                    if (orthogonalCorrection)
                    {
                        Tuple2<scalar, scalar> dxdyFromSources =
                            inverseDistanceInterpolationDzBottom(
                                pEval_2D,
                                globalPointsX,
                                globalPointsY,
                                globalDx,
                                globalDy,
                                interpRelRadius,
                                distThr);
                        interpDx += dxdyFromSources.first();
                        interpDy += dxdyFromSources.second();
                    }
                }
                else
                {
                    DeltaInterp =
                        inverseDistanceInterpolationDz(Ldef,
                                                       alphaAll,
                                                       coeffVertDeformation,
                                                       pEval,
                                                       globalPointsX,
                                                       globalPointsY,
                                                       globalPointsZ,
                                                       globalDz,
                                                       globalDx,
                                                       globalDy,
                                                       globalAreas,
                                                       distThr);

                    if (pEval.z() > noDeformLevel)
                    {
                        coeffHor = (zMax - pEval.z()) / (zMax - noDeformLevel);
                    }
                    else
                    {
                        coeffHor = 1.0;
                    }

                    coeffHor = Foam::pow(coeffHor, noDeformExp);

                    interpDx = DeltaInterp.x() * coeffHor;
                    interpDy = DeltaInterp.y() * coeffHor;
                    interpDz =
                        coeffHor * DeltaInterp.z() + (1.0 - coeffHor) * maxTopo;
                }
            }

            if (pEval.z() > 1e-2 * distThr)
            {
                dxMin_rel = (pEval.x() - xMin) / (xMax - xMin);
                dxMax_rel = (xMax - pEval.x()) / (xMax - xMin);
                dyMin_rel = (pEval.y() - yMin) / (yMax - yMin);
                dyMax_rel = (yMax - pEval.y()) / (yMax - yMin);

                if (distC2 > distC1)
                {
                    distC = Foam::sqrt(pow(pEval.x(), 2) + pow(pEval.y(), 2));
                    distCoeff = max(
                        0.0, min(1.0, (distC - distC1) / (distC2 - distC1)));
                }
                else
                {
                    distCoeff = 1.0;
                }

                xCoeff = min(distCoeff,
                             max(0.0,
                                 min(1.0,
                                     (min(dxMin_rel, dxMax_rel) - dist_rel1)
                                         / (dist_rel2 - dist_rel1))));
                yCoeff = min(distCoeff,
                             max(0.0,
                                 min(1.0,
                                     (min(dyMin_rel, dyMax_rel) - dist_rel1)
                                         / (dist_rel2 - dist_rel1))));

                pDeform[pointi].x() = xCoeff * interpDx * pEval.z();
                pDeform[pointi].y() = yCoeff * interpDy * pEval.z();
            }
            else
            {
                pDeform[pointi].x() = interpDx;
                pDeform[pointi].y() = interpDy;
            }
            pDeform[pointi].z() = interpDz;
        }
        Info << "Finished point deformation interpolation." << endl;

        // Calculate the new proposed point positions in a temporary field.
        // We are not modifying the mesh itself yet.
        tmp<pointField> tnewPoints = mesh.points() + pDeform;

        // Synchronize the positions of shared points across processors.
        // For each boundary point, this calculates the average of the positions
        // proposed by each processor that shares the point.
        // This "stitches" the mesh back together and eliminates discrepancies.
        syncTools::syncPointPositions(
            mesh, tnewPoints.ref(), maxEqOp<point>(), point::zero);

        // Now that the new positions are consistent, update the mesh.
        // We use movePoints, which is more robust than setPoints for this.
        mesh.movePoints(tnewPoints());

        Sout << "Proc" << Pstream::myProcNo() << " mesh updated" << endl;

        if (checkMesh)
        {
            const faceList& pFaces = mesh.faces();
            const pointField& pPts = mesh.points();
            const vectorField& pC = mesh.cellCentres();
            const labelList& pOwner = mesh.faceOwner();

            scalar minQ(1.0);

            forAll(z0FaceIndices, facei)
            {
                const face& f = pFaces[z0FaceIndices[facei]];
                scalar flatness = calculateFlatness(f, pPts);

                if (flatness < 0.98)
                {
                    Sout << "Proc" << Pstream::myProcNo() << " face " << facei
                         << " flatness " << flatness << endl;
                }

                label oCI = pOwner[z0FaceIndices[facei]];

                point oCc = pC[oCI];

                minQ = 1.0;

                forAll(f, faceBasePtI)
                {
                    minQ = minQuality(
                        mesh, oCc, z0FaceIndices[facei], true, faceBasePtI);
                }

                if (minQ < 1e-15)
                {
                    Sout << "Proc" << Pstream::myProcNo() << " face " << facei
                         << " minQ " << minQ << endl;
                    forAll(f, ip)
                    {
                        Sout << ip << " coord " << pPts[f[ip]] << endl;
                    }
                    Sout << " oCc " << oCc << endl;
                }
            }
        }
    }

    const Switch smoothing =
        topoDict.lookupOrDefault<Switch>("smoothing", false);

    if (smoothing)
    {
        // --- Input Parameters ---
        const int nSmoothIter = topoDict.lookupOrDefault<int>("nIter", 50);
        const scalar maxRotationAngleDeg =
            topoDict.lookupOrDefault<scalar>("maxRotationAngleDeg", 1.0);
        const scalar qualityThresholdDeg =
            topoDict.lookupOrDefault<scalar>("stopOnQualityDeg", 85.0);
        const int laplacianFrequency =
            topoDict.lookupOrDefault<int>("laplacianFrequency", 10);
        const scalar laplacianRelaxFactor =
            topoDict.lookupOrDefault<scalar>("laplacianRelaxFactor", 0.01);
        const scalar blendingFactor =
            topoDict.lookupOrDefault<scalar>("internalBlending", 0.2);
        const Switch normalAwareNearGroundSmoothing =
            topoDict.lookupOrDefault<Switch>("normalAwareNearGroundSmoothing",
                                             false);
        const label normalAwareLayers =
            topoDict.lookupOrDefault<label>("normalAwareLayers", 2);
        const scalar normalAwareNormalRelax =
            topoDict.lookupOrDefault<scalar>("normalAwareNormalRelax", 1.0);

        const Switch qualityAwareLaplacian =
            topoDict.lookupOrDefault<Switch>("qualityAwareLaplacian", false);
        const Switch rejectLaplacianOnConcave =
            topoDict.lookupOrDefault<Switch>("rejectLaplacianOnConcave", true);
        const Switch rejectLaplacianOnLowTetQuality =
            topoDict.lookupOrDefault<Switch>("rejectLaplacianOnLowTetQuality",
                                             true);
        const scalar minAllowedTetQuality =
            topoDict.lookupOrDefault<scalar>("minAllowedTetQuality", 1e-15);
        const label maxAllowedConcaveCells =
            topoDict.lookupOrDefault<label>("maxAllowedConcaveCells", -1);
        const Switch laplacianBacktracking =
            topoDict.lookupOrDefault<Switch>("laplacianBacktracking", false);
        const scalar laplacianBacktrackingFactor =
            topoDict.lookupOrDefault<scalar>("laplacianBacktrackingFactor",
                                             0.5);
        const scalar laplacianMinRelaxFactor =
            topoDict.lookupOrDefault<scalar>("laplacianMinRelaxFactor", 0.003);
        const scalar concavityTolerance = topoDict.lookupOrDefault<scalar>(
            "qualityAwareConcavityTolerance", SMALL);

        const Switch bestStateUseQualityScore =
            topoDict.lookupOrDefault<Switch>("bestStateUseQualityScore", false);
        const scalar bestStateConcaveWeight =
            topoDict.lookupOrDefault<scalar>("bestStateConcaveWeight", 0.005);
        const scalar bestStateLowTetWeight =
            topoDict.lookupOrDefault<scalar>("bestStateLowTetWeight", 1000.0);

        const scalar qualityCosThreshold =
            Foam::cos(degToRad(qualityThresholdDeg));

        // --- User Information ---
        Info << "\nImproving mesh quality using hybrid smoother..." << endl;
        Info << "  - Max smoothing iterations: " << nSmoothIter << endl;
        Info << "  - Max rotation angle (deg): " << maxRotationAngleDeg << endl;
        Info << "  - Stopping when worst non-orthogonality is below "
             << qualityThresholdDeg << " degrees." << endl;
        if (qualityAwareLaplacian)
        {
            Info << "  - Quality-aware Laplacian: enabled" << endl;
            Info << "    reject on concave: " << rejectLaplacianOnConcave
                 << ", reject on low tet quality: "
                 << rejectLaplacianOnLowTetQuality << endl;
            Info << "    minAllowedTetQuality: " << minAllowedTetQuality
                 << ", maxAllowedConcaveCells: " << maxAllowedConcaveCells
                 << endl;
            Info << "    backtracking: " << laplacianBacktracking
                 << ", factor: " << laplacianBacktrackingFactor
                 << ", min relax factor: " << laplacianMinRelaxFactor << endl;
        }
        if (bestStateUseQualityScore)
        {
            Info << "  - Best-state quality score: enabled" << endl;
            Info << "    score = maxNonOrtho[deg] + " << bestStateConcaveWeight
                 << "*nConcave + lowTetPenalty" << endl;
            Info << "    lowTet penalty weight: " << bestStateLowTetWeight
                 << endl;
        }

        // --- 0: Identify all boundary points that must remain fixed ---
        Info << "Identifying fixed boundary points..." << endl;
        boolList isBoundaryPoint(
            mesh.nPoints(),
            false); // Lista di flag per tutti i punti locali della mesh
        forAll(mesh.boundaryMesh(),
               patchI) // Itera su TUTTI i patch di confine locali
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            // --- NUOVA LOGICA: Filtra i patch per tipo ---
            // Un patch è "fisso" se NON è un processorPolyPatch.
            // Quindi, include i patch fisici e i cyclicPolyPatch.
            if (!isA<processorPolyPatch>(pp))
            {
                // Se non è un processorPolyPatch, allora è un patch fisico
                // (wall, inlet, outlet, top, bottom) o un cyclicPolyPatch. I
                // suoi punti devono rimanere fissi.
                Info << "  - Considering patch '" << pp.name()
                     << "' (type: " << pp.type()
                     << ") as a fixed boundary (physical or cyclic)." << endl;
                const labelUList& meshPts =
                    pp.meshPoints(); // Punti della patch attuale (indici
                                     // locali)
                forAll(
                    meshPts,
                    i) // Itera su tutti i punti che appartengono a questa patch
                {
                    isBoundaryPoint[meshPts[i]] =
                        true; // Marca il punto come "boundary point"
                }
            }
            else // Caso else: è un processorPolyPatch
            {
                Info << "  - Patch '" << pp.name() << "' (type: " << pp.type()
                     << ") is a processor boundary. Points are NOT fixed here."
                     << endl;
            }
        }

        // --- 0b. Optional normal-aware mask for near-terrain Laplacian
        // smoothing --- The mask is built topologically from wall patches, e.g.
        // terrain_in/terrain_out. Points in the first normalAwareLayers layers
        // are allowed to smooth only tangentially to the local terrain normal.
        boolList normalAwarePoint(mesh.nPoints(), false);
        vectorField normalAwareNormal(mesh.nPoints(), vector::zero);

        if (normalAwareNearGroundSmoothing && normalAwareLayers > 0)
        {
            Info << "Building normal-aware smoothing mask from wall patches..."
                 << endl;

            boolList frontier(mesh.nPoints(), false);

            const faceList& faces = mesh.faces();
            const vectorField& faceAreas = mesh.faceAreas();

            forAll(mesh.boundaryMesh(), patchI)
            {
                const polyPatch& pp = mesh.boundaryMesh()[patchI];

                if (!isA<processorPolyPatch>(pp) && pp.type() == "wall")
                {
                    Info << "  - Using patch '" << pp.name()
                         << "' for terrain-normal smoothing." << endl;

                    forAll(pp, i)
                    {
                        const label faceI = pp.start() + i;
                        vector n = faceAreas[faceI];

                        if (mag(n) > VSMALL)
                        {
                            n /= mag(n);

                            const face& f = faces[faceI];
                            forAll(f, fp)
                            {
                                const label pI = f[fp];
                                frontier[pI] = true;
                                normalAwareNormal[pI] += n;
                            }
                        }
                    }
                }
            }

            syncTools::syncPointList(mesh, frontier, orEqOp<bool>(), false);
            syncTools::syncPointList(
                mesh, normalAwareNormal, plusEqOp<vector>(), vector::zero);

            forAll(normalAwareNormal, pI)
            {
                if (mag(normalAwareNormal[pI]) > VSMALL)
                {
                    normalAwareNormal[pI] /= mag(normalAwareNormal[pI]);
                }
            }

            const edgeList& edges = mesh.edges();

            for (label layer = 0; layer < normalAwareLayers; ++layer)
            {
                boolList nextFrontier(mesh.nPoints(), false);
                vectorField nextNormal(mesh.nPoints(), vector::zero);

                forAll(edges, edgeI)
                {
                    const edge& e = edges[edgeI];
                    const label p0 = e[0];
                    const label p1 = e[1];

                    if (frontier[p0] && !isBoundaryPoint[p1])
                    {
                        normalAwarePoint[p1] = true;
                        nextFrontier[p1] = true;
                        nextNormal[p1] += normalAwareNormal[p0];
                    }

                    if (frontier[p1] && !isBoundaryPoint[p0])
                    {
                        normalAwarePoint[p0] = true;
                        nextFrontier[p0] = true;
                        nextNormal[p0] += normalAwareNormal[p1];
                    }
                }

                syncTools::syncPointList(
                    mesh, normalAwarePoint, orEqOp<bool>(), false);
                syncTools::syncPointList(
                    mesh, nextFrontier, orEqOp<bool>(), false);
                syncTools::syncPointList(
                    mesh, nextNormal, plusEqOp<vector>(), vector::zero);

                forAll(nextNormal, pI)
                {
                    if (mag(nextNormal[pI]) > VSMALL)
                    {
                        nextNormal[pI] /= mag(nextNormal[pI]);
                        normalAwareNormal[pI] = nextNormal[pI];
                    }
                }

                frontier = nextFrontier;
            }

            label nProtected = 0;
            forAll(normalAwarePoint, pI)
            {
                if (normalAwarePoint[pI])
                {
                    ++nProtected;
                }
            }
            reduce(nProtected, sumOp<label>());

            Info << "  - Normal-aware Laplacian points: " << nProtected << endl;
        }

        // --- 1. Initial Calculation and Setup for Best-State Tracking ---
        mesh.clearGeom();
        tmp<scalarField> t_ortho = meshCheck::faceOrthogonality(
            mesh, mesh.faceAreas(), mesh.cellCentres());
        scalarField& ortho = t_ortho.ref();

        // --- BEST MESH LOGIC: Initialize tracking variables ---
        scalar initialMinOrtho = GREAT;
        for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
        {
            initialMinOrtho = min(initialMinOrtho, ortho[faceI]);
        }
        reduce(initialMinOrtho, minOp<scalar>());

        pointField bestPoints = mesh.points();
        scalar bestMinOrtho = initialMinOrtho;
        label bestConcaveCells = 0;
        scalar bestMinTetQuality = GREAT;
        scalar bestScore =
            Foam::radToDeg(Foam::acos(max(-1.0, min(1.0, bestMinOrtho))));

        if (bestStateUseQualityScore)
        {
            mesh.clearGeom();
            bestConcaveCells =
                countConcaveCellsSimple(mesh, concavityTolerance);
            bestMinTetQuality = minFaceTetQualitySimple(mesh);
            bestScore += bestStateConcaveWeight * bestConcaveCells;
            if (bestMinTetQuality < minAllowedTetQuality)
            {
                bestScore += bestStateLowTetWeight;
            }
        }

        const scalarList normFactors = calculateEdgeNormFactors(mesh);
        const scalarList neighbourCount_points =
            calculateNeighbourCount(mesh, normFactors);
        const scalarList neighbourCount_cells = calculateGlobalCellCount(mesh);

        if (Pstream::master())
        {
            scalar cosVal = max(-1.0, min(1.0, initialMinOrtho));
            scalar initialWorstAngle = Foam::radToDeg(Foam::acos(cosVal));
            Info << "  - Initial worst non-orthogonality: " << initialWorstAngle
                 << " deg." << endl;
        }
        // --- END OF BEST MESH LOGIC ---

        /*
        const point faceA_centroid(
            -112.692518067, -229.315478954, 389.879622638);

        label labelA = -1;
        label procA = -1;

        for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
        {
            const point& currentCentroid = mesh.faceCentres()[faceI];

            // Controlliamo se la faccia attuale è una di quelle che
            // cerchiamo
            if (mag(currentCentroid - faceA_centroid) < 1e-4)
            {
                Sout << "Trovata Face A su proc " << Pstream::myProcNo()
                     << " faccia " << faceI << endl;
                labelA = faceI;
                procA = Pstream::myProcNo();
            }
        }

        reduce(procA, maxOp<label>());

        const point faceB_centroid(
            -201.182612077, -131.209119549, 373.489510406);

        label labelB = -1;
        label procB = -1;

        for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
        {
            const point& currentCentroid = mesh.faceCentres()[faceI];

            // Controlliamo se la faccia attuale è una di quelle che
            // cerchiamo
            if (mag(currentCentroid - faceB_centroid) < 1e-4)
            {
                Sout << "Trovata Face B su proc " << Pstream::myProcNo()
                     << " faccia " << faceI << endl;
                labelB = faceI;
                procB = Pstream::myProcNo();
            }
        }

        label labelA0 = -1;
        label labelB0 = -1;

        reduce(procB, maxOp<label>());
        */
        scalar currentMaxRotDeg = maxRotationAngleDeg;
        scalar lastLaplacianBestMinOrtho = bestMinOrtho;
        Info << "  Adaptive rotation enabled. Starting angle: "
             << currentMaxRotDeg << " deg." << endl;

        // --- Main Hybrid Smoothing Loop ---
        for (int iter = 0; iter < nSmoothIter; ++iter)
        {
            bool doLaplacianSmoothing =
                (laplacianFrequency > 0
                 && (iter + 1) % laplacianFrequency == 0);

            if (doLaplacianSmoothing)
            {


                if (bestMinOrtho > lastLaplacianBestMinOrtho
                                       + 1e-8) // Usiamo una piccola tolleranza
                {
                    currentMaxRotDeg *= 1.01;
                    // CASO POSITIVO: Abbiamo fatto progressi in questo blocco
                    // di iterazioni.
                    if (Pstream::master())
                    {
                        Info << "  Iteration " << iter + 1
                             << " - Progress detected since last Laplacian ("
                             << lastLaplacianBestMinOrtho << " -> "
                             << bestMinOrtho
                             << "). Keeping max rotation: " << currentMaxRotDeg
                             << " deg." << endl;
                    }
                    // Opzionale: Potresti voler ripristinare l'angolo originale
                    // se le cose vanno bene? currentMaxRotDeg =
                    // maxRotationAngleDeg;
                }
                else
                {
                    // CASO NEGATIVO: Nessun record battuto nell'ultimo ciclo.
                    // Stagnazione o Deadlock. Puniamo l'algoritmo riducendo
                    // l'angolo.

                    scalar oldAngle = currentMaxRotDeg;
                    currentMaxRotDeg *= 0.99; // Riduci del 25%

                    if (currentMaxRotDeg < 1e-4)
                        currentMaxRotDeg = 1e-4; // Limite inferiore

                    if (Pstream::master())
                    {
                        Info << "  Iteration " << iter + 1
                             << " - NO Progress since last Laplacian. Reducing "
                                "max rotation: "
                             << oldAngle << " -> " << currentMaxRotDeg
                             << " deg." << endl;
                    }
                }

                // 2. AGGIORNA IL CHECKPOINT PER IL PROSSIMO GIRO
                lastLaplacianBestMinOrtho = bestMinOrtho;

                // ======================================================= //
                //             LAPLACIAN SMOOTHING ITERATION               //
                // ======================================================= //

                // --- BEST MESH LOGIC: Revert to best state before Laplacian
                // step --- Questo blocco rimane invariato e funziona come da
                // tua logica originale.
                if (Pstream::master())
                {
                    scalar bestAngleSoFar = Foam::radToDeg(
                        Foam::acos(max(-1.0, min(1.0, bestMinOrtho))));
                    Info << "  Iteration " << iter + 1
                         << " - Reverting to best state (quality: "
                         << bestAngleSoFar << " deg)"
                         << " before applying Laplacian smoothing." << endl;
                }

                /*
                if (Pstream::myProcNo() == procB)
                {
                    const label worstFaceI = labelB;

                    Sout << "faceB center " << mesh.faceCentres()[worstFaceI]
                         << endl;

                    const face& faceB = mesh.faces()[worstFaceI];
                    forAll(faceB, fp)
                    {
                        label pI = faceB[fp];
                        labelB0 = pI;
                        const point& p = mesh.points()[pI];
                        Sout << "pointB " << pI << " coords " << p << endl;
                    }
                }


                if (Pstream::myProcNo() == procA)
                {
                    const label worstFaceI = labelA;

                    Sout << "faceA center " << mesh.faceCentres()[worstFaceI]
                         << endl;

                    const face& faceA = mesh.faces()[worstFaceI];
                    forAll(faceA, fp)
                    {
                        label pI = faceA[fp];
                        labelA0 = pI;
                        const point& p = mesh.points()[pI];
                        Sout << "pointA " << pI << " coords " << p << endl;
                    }
                }
                */


                const_cast<pointField&>(mesh.points()) = bestPoints;
                syncTools::syncPointPositions(
                    mesh,
                    const_cast<pointField&>(mesh.points()),
                    minEqOp<point>(),
                    point(great, great, great));
                mesh.clearGeom();
                // --- END OF BEST MESH LOGIC ---

                if (Pstream::master())
                {
                    Info
                        << "  Iteration " << iter + 1
                        << " - Performing global blended Laplacian smoothing..."
                        << endl;
                }

                // --- 1. Compute Target Laplacian Points (Point-based) ---
                tmp<pointField> tLaplacianTarget = getLaplacianTargetPoints(
                    mesh, normFactors, neighbourCount_points);
                const pointField& P_laplacian_target = tLaplacianTarget();

                /*
                if (Pstream::myProcNo() == procA)
                {
                    Sout << "targetA " << labelA0 << " "
                         << P_laplacian_target[labelA0] << endl;
                }

                if (Pstream::myProcNo() == procB)
                {
                    Sout << "targetB " << labelB0 << " "
                         << P_laplacian_target[labelB0] << endl;
                }
                */


                // --- 2. Compute Target Centroidal Points (Cell-based) ---
                tmp<pointField> tCentroidalTarget = getCentroidalTargetPoints(
                    mesh); // Passa il global cell count
                const pointField& P_internal_target = tCentroidalTarget();


                // --- 3. Compute blended displacement ---
                boolList pointsToMove(mesh.nPoints(), false);
                forAll(pointsToMove, pI)
                {
                    if (!isBoundaryPoint[pI])
                    {
                        pointsToMove[pI] = true;
                    }
                }

                pointField proposedDisplacement(mesh.nPoints(), vector::zero);

                forAll(pointsToMove, pointI)
                {
                    if (pointsToMove[pointI])
                    {
                        const point& currentPos = mesh.points()[pointI];

                        point P_ideal_blended =
                            (1.0 - blendingFactor) * P_laplacian_target[pointI]
                            + blendingFactor * P_internal_target[pointI];

                        proposedDisplacement[pointI] =
                            laplacianRelaxFactor
                            * (P_ideal_blended - currentPos);

                        if (normalAwareNearGroundSmoothing
                            && normalAwarePoint[pointI]
                            && mag(normalAwareNormal[pointI]) > VSMALL)
                        {
                            const vector nHat = normalAwareNormal[pointI]
                                              / mag(normalAwareNormal[pointI]);

                            proposedDisplacement[pointI] -=
                                normalAwareNormalRelax
                                * (proposedDisplacement[pointI] & nHat) * nHat;
                        }
                    }
                }

                // --- 4. Sync the displacements among the processors
                syncTools::syncPointList(mesh,
                                         proposedDisplacement,
                                         minEqOp<point>(),
                                         point(great, great, great));


                /*
                if (Pstream::myProcNo() == procA)
                {
                    Sout << "prop displA " << labelA0 << " "
                         << proposedDisplacement[labelA0] << endl;
                }

                if (Pstream::myProcNo() == procB)
                {
                    Sout << "prop displB " << labelB0 << " "
                         << proposedDisplacement[labelB0] << endl;
                }
                */


                // --- 5. Apply the Laplacian step.  Optionally reject or
                // backtrack if the step worsens checkMesh-like quality metrics.
                const pointField pointsBeforeLaplacian = mesh.points();

                bool acceptLaplacianStep = true;
                scalar acceptedScale = 1.0;

                label oldConcaveCells = 0;
                scalar oldMinTetQuality = GREAT;

                if (qualityAwareLaplacian)
                {
                    mesh.clearGeom();
                    oldConcaveCells =
                        countConcaveCellsSimple(mesh, concavityTolerance);
                    oldMinTetQuality = minFaceTetQualitySimple(mesh);
                }

                scalar trialScale = 1.0;

                while (true)
                {
                    const_cast<pointField&>(mesh.points()) =
                        pointsBeforeLaplacian
                        + trialScale * proposedDisplacement;

                    syncTools::syncPointPositions(
                        mesh,
                        const_cast<pointField&>(mesh.points()),
                        minEqOp<point>(),
                        point(great, great, great));

                    if (!qualityAwareLaplacian)
                    {
                        acceptLaplacianStep = true;
                        acceptedScale = trialScale;
                        break;
                    }

                    mesh.clearGeom();
                    const label newConcaveCells =
                        countConcaveCellsSimple(mesh, concavityTolerance);
                    const scalar newMinTetQuality =
                        minFaceTetQualitySimple(mesh);

                    bool badStep = false;

                    if (rejectLaplacianOnConcave
                        && newConcaveCells > oldConcaveCells)
                    {
                        badStep = true;
                    }

                    if (rejectLaplacianOnConcave && maxAllowedConcaveCells >= 0
                        && newConcaveCells > maxAllowedConcaveCells)
                    {
                        badStep = true;
                    }

                    if (rejectLaplacianOnLowTetQuality
                        && newMinTetQuality < minAllowedTetQuality)
                    {
                        badStep = true;
                    }

                    if (!badStep)
                    {
                        acceptLaplacianStep = true;
                        acceptedScale = trialScale;

                        if (Pstream::master())
                        {
                            Info << "    - Accepted quality-aware Laplacian"
                                 << " scale " << acceptedScale << ": concave "
                                 << oldConcaveCells << " -> " << newConcaveCells
                                 << ", minTetQuality " << oldMinTetQuality
                                 << " -> " << newMinTetQuality << endl;
                        }
                        break;
                    }

                    if (laplacianBacktracking
                        && trialScale * laplacianRelaxFactor
                               > laplacianMinRelaxFactor + SMALL)
                    {
                        if (Pstream::master())
                        {
                            Info << "    - Rejecting Laplacian trial scale "
                                 << trialScale << ": concave "
                                 << oldConcaveCells << " -> " << newConcaveCells
                                 << ", minTetQuality " << oldMinTetQuality
                                 << " -> " << newMinTetQuality
                                 << ". Backtracking." << endl;
                        }

                        trialScale *= laplacianBacktrackingFactor;
                        continue;
                    }

                    acceptLaplacianStep = false;

                    if (Pstream::master())
                    {
                        Info << "    - Rejected quality-aware Laplacian"
                             << ": concave " << oldConcaveCells << " -> "
                             << newConcaveCells << ", minTetQuality "
                             << oldMinTetQuality << " -> " << newMinTetQuality
                             << endl;
                    }
                    break;
                }

                if (!acceptLaplacianStep)
                {
                    const_cast<pointField&>(mesh.points()) =
                        pointsBeforeLaplacian;

                    syncTools::syncPointPositions(
                        mesh,
                        const_cast<pointField&>(mesh.points()),
                        minEqOp<point>(),
                        point(great, great, great));
                }
            }
            else
            {
                // ======================================================= //
                //              ROTATIONAL SMOOTHING ITERATION             //
                // ======================================================= //

                scalar localMinOrtho = GREAT;
                label localWorstFaceI = -1;
                point localWorstFaceCentroid = point::max;

                for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
                {
                    if (ortho[faceI] < localMinOrtho)
                    {
                        localMinOrtho = ortho[faceI];
                        localWorstFaceI = faceI;
                        localWorstFaceCentroid =
                            mesh.faceCentres()[faceI]; // Save the centroid
                    }
                }

                // --- PHASE 1: Find the global min non-orthogonality ---
                scalar globalMinOrtho = localMinOrtho;
                reduce(globalMinOrtho, minOp<scalar>());


                // --- PHASE 2: Find the min Z coord among the faces with
                // globalMinOrtho ---
                scalar localMinZCoord = GREAT; // Initialize with large value
                if (localWorstFaceI != -1
                    && mag(localMinOrtho - globalMinOrtho) < SMALL)
                {
                    localMinZCoord = localWorstFaceCentroid.z();
                }
                scalar globalMinZCoord = localMinZCoord;
                reduce(globalMinZCoord, minOp<scalar>());


                // --- PHASE 3: Find the master processor among those
                // satisfying both the criteria ---
                label masterProcID = -1; // Default: no processor is master
                if (localWorstFaceI != -1
                    && mag(localMinOrtho - globalMinOrtho) < 1e-6
                    && mag(localMinZCoord - globalMinZCoord) < 1e-6)
                {
                    masterProcID =
                        Pstream::myProcNo(); // This processor is a candidate
                    /*
                    Sout << "  Iteration " << iter + 1
                         << "myProc " << masterProcID << " localWorstFaceI "
                         << localWorstFaceI << " localMinOrtho "
                         << localMinOrtho << " localMinZCoord "
                         << localMinZCoord << endl;
                    */
                }
                // Select the processor with the lowest ID among the candidates
                // for deterministic tie-breaking. If your MPI/OpenFOAM
                // implementation prefers maxOp for tie-breaking, use that
                // instead. We use minOp here for the procID for a deterministic
                // tie-break (e.g., proc 0, then proc 1, etc.).
                reduce(masterProcID, maxOp<label>());

                if (globalMinOrtho > qualityCosThreshold)
                {
                    Info << "  Mesh quality is acceptable. Stopping smoothing."
                         << endl;
                    break;
                }

                List<point> neighbourCellCentres;
                syncTools::swapBoundaryCellPositions(
                    mesh, mesh.cellCentres(), neighbourCellCentres);

                pointField displacement(mesh.nPoints(), vector::zero);

                if (Pstream::myProcNo() == masterProcID)
                {
                    const label worstFaceI = localWorstFaceI;
                    const label ownI = mesh.faceOwner()[worstFaceI];
                    const point& ownC = mesh.cellCentres()[ownI];
                    point neiC;

                    Sout << "  Iteration " << iter + 1 << " Worst angle "
                         << localMinOrtho << " Worst face center "
                         << mesh.faceCentres()[worstFaceI] << endl;

                    if (worstFaceI >= mesh.nInternalFaces())
                    {
                        // The face is on the boundary of the local mesh
                        if (isProcessorFace(mesh, worstFaceI))
                        {
                            Sout << "Face Type: processor boundary face"
                                 << endl;
                        }
                        else
                        {
                            const label patchID =
                                mesh.boundaryMesh().whichPatch(worstFaceI);
                            Sout << "Face Type: physical boundary face on "
                                    "patch '"
                                 << mesh.boundaryMesh()[patchID].name() << "'"
                                 << endl;
                        }
                        Sout << "not internal face" << endl;
                        // Get the neig.centroid from other proc.
                        neiC = neighbourCellCentres[worstFaceI
                                                    - mesh.nInternalFaces()];
                    }
                    else
                    {
                        // The face is internal for the local mesh
                        Sout << "internal face" << endl;
                        neiC = mesh.cellCentres()
                                   [mesh.faceNeighbour()[worstFaceI]];
                    }

                    const vector& S = mesh.faceAreas()[worstFaceI];
                    const point& fC = mesh.faceCentres()[worstFaceI];
                    vector d = neiC - ownC;

                    if (mag(d) > SMALL && mag(S) > SMALL)
                    {
                        vector S_hat = normalised(S);
                        vector d_hat = normalised(d);
                        scalar cosAngle = S_hat & d_hat;

                        if (cosAngle < 0.9999)
                        {
                            vector axis = S_hat ^ d_hat;
                            if (mag(axis) > SMALL)
                            {
                                scalar totalCorrectionAngle =
                                    Foam::acos(max(-1.0, min(1.0, cosAngle)));
                                scalar currentMaxRotRad =
                                    Foam::degToRad(currentMaxRotDeg);
                                scalar rotationStepAngle =
                                    min(totalCorrectionAngle, currentMaxRotRad);

                                // scalar rotationStepAngle = min(
                                //     totalCorrectionAngle,
                                //     maxRotationAngleRad);
                                quaternion R(normalised(axis),
                                             rotationStepAngle);

                                const face& worstFace =
                                    mesh.faces()[worstFaceI];
                                forAll(worstFace, fp)
                                {
                                    label pI = worstFace[fp];
                                    if (!isBoundaryPoint[pI])
                                    {
                                        const point& p = mesh.points()[pI];
                                        vector rotated_vec =
                                            R.transform(p - fC);
                                        point p_new = fC + rotated_vec;
                                        displacement[pI] = p_new - p;
                                    }
                                    else
                                    {
                                        // Sout << "boundary point " << endl;
                                    }
                                }
                            }
                        }
                    }
                }

                syncTools::syncPointList(mesh,
                                         displacement,
                                         plusEqOp<vector>(),
                                         vector::zero // NULLVALUE
                );

                const_cast<pointField&>(mesh.points()) += displacement;
                syncTools::syncPointPositions(
                    mesh,
                    const_cast<pointField&>(mesh.points()),
                    minEqOp<point>(),
                    point(great, great, great));
            }

            // --- Recalculate orthogonality and check for new best state ---
            mesh.clearGeom();
            tmp<scalarField> t_ortho_new = meshCheck::faceOrthogonality(
                mesh, mesh.faceAreas(), mesh.cellCentres());
            ortho = t_ortho_new.ref();

            // --- BEST MESH LOGIC: Check and save if current state is better
            scalar localMinOrthoForBestMesh = GREAT;
            for (label faceI = 0; faceI < mesh.nFaces(); ++faceI)
            {
                localMinOrthoForBestMesh =
                    min(localMinOrthoForBestMesh, ortho[faceI]);
            }

            // --- Global aggregation
            scalar currentMinOrtho = localMinOrthoForBestMesh;
            reduce(currentMinOrtho, minOp<scalar>());

            const scalar currentMaxNonOrthoDeg = Foam::radToDeg(
                Foam::acos(max(-1.0, min(1.0, currentMinOrtho))));

            bool currentIsBest = false;

            if (bestStateUseQualityScore)
            {
                mesh.clearGeom();
                const label currentConcaveCells =
                    countConcaveCellsSimple(mesh, concavityTolerance);
                const scalar currentMinTetQuality =
                    minFaceTetQualitySimple(mesh);

                scalar currentScore =
                    currentMaxNonOrthoDeg
                    + bestStateConcaveWeight * currentConcaveCells;

                if (currentMinTetQuality < minAllowedTetQuality)
                {
                    currentScore += bestStateLowTetWeight;
                }

                if (currentScore < bestScore)
                {
                    currentIsBest = true;
                    bestScore = currentScore;
                    bestConcaveCells = currentConcaveCells;
                    bestMinTetQuality = currentMinTetQuality;
                }
            }
            else if (currentMinOrtho > bestMinOrtho)
            {
                currentIsBest = true;
            }

            if (currentIsBest)
            {
                bestMinOrtho = currentMinOrtho;
                bestPoints = mesh.points();
                if (Pstream::master())
                {
                    if (bestStateUseQualityScore)
                    {
                        Info << "    - New best mesh quality score found at "
                                "iteration "
                             << iter + 1 << ": score " << bestScore
                             << ", nonOrtho " << currentMaxNonOrthoDeg
                             << " deg, concave " << bestConcaveCells
                             << ", minTetQuality " << bestMinTetQuality << endl;
                    }
                    else
                    {
                        Info
                            << "    - New best mesh quality found at iteration "
                            << iter + 1 << ": " << currentMaxNonOrthoDeg
                            << " deg." << endl;
                    }
                }
            }
            else
            {
                Info << "curretnMaxAngle " << currentMaxNonOrthoDeg << " deg."
                     << endl;
            }
            // --- END OF BEST MESH LOGIC ---
        }

        // --- BEST MESH LOGIC: Final restoration to the best found
        // configuration ---
        Info << "\nSmoothing iterations completed." << endl;
        const_cast<pointField&>(mesh.points()) = bestPoints;
        syncTools::syncPointPositions(mesh,
                                      const_cast<pointField&>(mesh.points()),
                                      minEqOp<point>(),
                                      point(great, great, great));
        mesh.clearGeom();

        if (Pstream::master())
        {
            scalar finalBestAngle =
                Foam::radToDeg(Foam::acos(max(-1.e-6, min(1.0, bestMinOrtho))));
            Info << "Restored mesh to best configuration with worst "
                    "non-orthogonality: "
                 << finalBestAngle << " deg." << endl;
            if (bestStateUseQualityScore)
            {
                Info << "Best-state score: " << bestScore
                     << ", concave cells: " << bestConcaveCells
                     << ", minTetQuality: " << bestMinTetQuality << endl;
            }
        }
        // --- END OF BEST MESH LOGIC ---
    }

    const Switch shiftMinZToZero =
        topoDict.lookupOrDefault<Switch>("shiftMinZToZero", false);

    if (shiftMinZToZero)
    {
        scalar zMinFinal = great;

        const pointField& finalPoints = mesh.points();
        forAll(finalPoints, pointi)
        {
            zMinFinal = min(zMinFinal, finalPoints[pointi].z());
        }
        reduce(zMinFinal, minOp<scalar>());

        Info << "Shifting final mesh vertically so that global min(z) is 0"
             << endl;
        Info << "  - Global final min(z) before shift: " << zMinFinal << endl;
        Info << "  - Applied vertical shift: " << -zMinFinal << endl;

        pointField shiftedPoints(mesh.points());
        forAll(shiftedPoints, pointi)
        {
            shiftedPoints[pointi].z() -= zMinFinal;
        }

        syncTools::syncPointPositions(
            mesh, shiftedPoints, maxEqOp<point>(), point::zero);

        mesh.movePoints(shiftedPoints);

        scalar zMinCheck = great;
        const pointField& shiftedCheckPoints = mesh.points();
        forAll(shiftedCheckPoints, pointi)
        {
            zMinCheck = min(zMinCheck, shiftedCheckPoints[pointi].z());
        }
        reduce(zMinCheck, minOp<scalar>());

        Info << "  - Global final min(z) after shift: " << zMinCheck << endl;
    }

    Info << "Writing new mesh" << endl;
    mesh.setInstance("constant");
    mesh.write();

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
         << endl;

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
