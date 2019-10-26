/*
Copyright (C) 2019 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <chrono>

using namespace std;

struct VCoord           //for 3D points or normal vectors
{
    double x, y, z;
};
struct Faces            //for triangular mesh faces
{
    int a, b, c;
};
struct LMcoord          //for landmark coordinates
{
    int name;
    double x, y, z;
};
struct PCoord           //for 2D projected points
{
    double x, y;
};
struct Mesh             //for mesh object data
{
    double V1x, V1y, V1z;
    double V2x, V2y, V2z;
    double V3x, V3y, V3z;
};
struct Raster
{
    bool pixel;
    vector<int> pointer;
};
struct Pixel
{
    int row, col;
    bool valid;
};
struct HeightMap
{
    int height;
    bool pixel;
};    
struct Contour
{
    int row, col;
    int direction;
    vector<int> pointer;
};
struct FCCode
{
    int chain;
    int delta_x;
    int delta_y;
    double delta_t;
};
struct EFDcoef
{
    double a, b, c, d;
};

// define namespace for landmark presence
namespace landmark
{
    bool bregma = false;
    bool mastoidaleL = false;
    bool mastoidaleR = false;
    bool Frankfurt = false;
    bool nasion_optimal = false;
}
// define namespace for global variables
namespace global
{
    string obj_filename;
    string pp_filename;
    bool print_out = true;
    bool print_oct = false;
    bool print_exp = false;
    bool print_lmk = false;
    VCoord opisthocranion3D;
    PCoord opisthocranion2D;
    VCoord opisthion3D;
    PCoord opisthion2D;
    VCoord nasion3D;
    PCoord nasion2D;
    VCoord bregma3D;
    PCoord bregma2D;
    double resolution = 0.5;
    Pixel InitialPixel;
    double left_upper_bound_y;
    double right_upper_bound_y;
}
// define namespace for Freeman chain codes
namespace FreemanChainCode
{
    vector<FCCode> NasionBregma;
}
// define namespace for Elliptic Fourier Descriptors
namespace EFD
{
    vector<EFDcoef> NasionBregma;
}
// define namespace for Height Map Images
namespace HMI
{
    vector<vector<HeightMap>> opHeightMapImage;
    vector<vector<HeightMap>> srHeightMapImage;
    vector<vector<HeightMap>> lmLatHeightMapImage;
    vector<vector<HeightMap>> lmInfHeightMapImage;
    vector<vector<HeightMap>> rmLatHeightMapImage;
    vector<vector<HeightMap>> rmInfHeightMapImage;
}

// function for checking user input arguments and declaring specific global variables
void checkUserArguments(int argc, char *argv[])
{
    string parameter;
    if(argc == 1)
    {
        // store the first argument in a string and check for --citation and --help parameters
        parameter = argv[0];
        string license = "--license";
        string citation = "--citation";
        string help = "--help";
        if(!parameter.compare(citation))
        {
            cout << "skullanalyzer v1.0 - A cranial morphology analysis software" << endl;
            cout << "Copyright (C), Andreas Bertsatos <abertsatos@biol.uoa.gr>, All rights reserved." << endl;
            cout << "To use skullanalyzer in publications use:" << endl;
            cout << "Bertsatos A (2019) Skullanalyzer: a concrete way of extracting cranial geometric features." << endl;
            exit(EXIT_SUCCESS);
        }
        else if(!parameter.compare(help))
        {
            cout << "\nskullanalyzer procceses cranial 3D triangular mesh models saved in Wavefront Alias format and" << endl;
            cout << "automatically calculates the appropriate Frankfurt position according to specific landmarks" << endl;
            cout << "provided in a separate Meshlab PickedPoints file, while it also calculates the optimal position" << endl;
            cout << "of particular landmarks to minimize observer error. The landmarks provided in the sidecar .pp" << endl;
            cout << "file (text file in Meshlab PickedPoints format) should be named with the following convention:\n" << endl;
            cout << "   Nasion             name=1" << endl;
            cout << "   Bregma             name=2" << endl;
            cout << "   Opisthion          name=3" << endl;
            cout << "   left Mastoidale    name=4" << endl;
            cout << "   right Mastoidale   name=5" << endl;
            cout << "   left Orbitale      name=6" << endl;
            cout << "   left Porion        name=7" << endl;
            cout << "   right Porion       name=8" << endl;
            cout << "\nThe minimum required landmarks, which are used for defining the Frankfurt position of the analyzed" << endl;
            cout << "cranium are Nasion, Opisthion, left Orbitale, and left and right Porion. According to availability" << endl;
            cout << "of additional landmarks the 'skullanalyzer' will calculale a number of geometric features from the" << endl;
            cout << "provided cranium and extract each feature in respective Comma Separated Values files as well as an" << endl;
            cout << "appropriate .mat data text format, which can be loaded for post-processing in GNU Octave programming" << endl;
            cout << "environment. All geometric features extracted in csv files can be visualized with the accompanying" << endl;
            cout << "'plot_features.m' function within the GNU Octave environment.\n" << endl;
            cout << "The geometric features extracted by skullanalyzer are the Nasion-Bregma midsaggital segment, whose" << endl;
            cout << "2D polyline is saved in the respective csv file and the corresponding Eliptic Fourier Descriptors are," << endl;
            cout << "saved in the .mat container, and a number of Height Map Images, which represent depth projections of" << endl;
            cout << "specific areas onto specific planes and are saved both in respective csv files as well as in the .mat" << endl;
            cout << "data container. The HMIs concern the supraorbital ridge, the occipital protuberance, and the lateral" << endl;
            cout << "and inferior projections of the left and right mastoid processes. More detailed description about the" << endl;
            cout << "orientation and calculation of the HMIs and EFDs can be found in the supplementary report provided" << endl;
            cout << "with the skullanalyzer. The report also contains details of the optimization algorithms used for the" << endl;
            cout << "Nasion and Mastoidale landmarks to define their correct position on the cranial surface.\n" << endl;
            cout << "The program can take up to 4 discrete and independent parameters bundled in a single string argument" << endl;
            cout << "beginning with a dash (-) as shown below. The four letters, 'p', 'e', 'o', 's', can be set in any" << endl;
            cout << "random order after the initial (-) character.\n" << endl;
            cout << "   $ skullanalyzer -peos model.obj landmarks.pp" << endl;
            cout << "\nThe following table lists the functionality of each parameter as well as its default value. Note that" << endl;
            cout << "different characters are ignored by the skullanalyzer.\n" << endl;
            cout << "Character          Functionality when set in parameter string            Default operation when missing\n" << endl;
            cout << "    p       Save the updated list of landmarks including the               Does not export updated" << endl;
            cout << "            optimized points into a new separate .pp file.                 landmark coordinates.\n" << endl;
            cout << "    e       Export the geometric features into the relevant .csv           Does not export the .csv" << endl;
            cout << "            files. This includes the 2D polyline of the nasion-bregma      file of any geometric" << endl;
            cout << "            segment as well as all available HMIs.                         feature.\n" << endl;
            cout << "    o       Save all calculated features into GNU Octave text data         Does not export the .mat" << endl;
            cout << "            format. This includes both the EFDs of the nasion-bregma       file with all calculated" << endl;
            cout << "            segment as well as all available HMIs.                         results.\n" << endl;
            cout << "    s       Silence output to terminal. Suitable for batch processing      All requested operations" << endl;
            cout << "            several samples through a bash script or other programs.       are reported in terminal." << endl;
            cout << "\nA compmlete User Manual & Algorithm Description can be found at https://doi.org/10.5281/zenodo.3519248" << endl;
            cout << "To cite the skullanalyzer in publications or see its license enter:\n" << endl;
            cout << "   $ skullanalyzer --citation" << endl;
            cout << "   $ skullanalyzer --license\n" << endl;
            exit(EXIT_SUCCESS);
        }
        else if(!parameter.compare(license))
        {
            cout << "Copyright (C) 2019 Andreas Bertsatos <abertsatos@biol.uoa.gr>\n" << endl;
            cout << "This program is free software; you can redistribute it and/or modify it under" << endl;
            cout << "the terms of the GNU General Public License as published by the Free Software" << endl;
            cout << "Foundation; either version 3 of the License, or (at your option) any later" << endl;
            cout << "version.\n" << endl;
            cout << "This program is distributed in the hope that it will be useful, but WITHOUT" << endl;
            cout << "ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or" << endl;
            cout << "FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more" << endl;
            cout << "details.\n" << endl;
            cout << "You should have received a copy of the GNU General Public License along with" << endl;
            cout << "his program; if not, see <http://www.gnu.org/licenses/>." << endl;
            exit(EXIT_SUCCESS);
        }
        else
        {
            cout << "Type 'skullanalyzer -- help' to see usage." << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(argc < 2)
    {
        cout << "Type 'skullanalyzer -- help' to see usage." << endl;
        exit(EXIT_FAILURE);
    }
    if(argc == 2)
    {
        // check if first argument is an .obj file and second argument is a .pp file
        parameter = argv[0];
        string obj = ".obj";
        if(parameter.find(obj) == string::npos)
        {
            cout << "First filename must point to a triangular mesh .obj file." << endl;
            exit(EXIT_FAILURE);
        }
        parameter = argv[1];
        string pp = ".pp";
        if(parameter.find(pp) == string::npos)
        {
            cout << "Second filename must point to the corresponding Meshlab point .pp file." << endl;
            exit(EXIT_FAILURE);
        }
        // store filename string of source and target obj files
        global::obj_filename = argv[0];
        global::pp_filename = argv[1];
    }
    if(argc ==3)
    {
        // check if second argument is an .obj file and third argument is a .pp file
        parameter = argv[1];
        string obj = ".obj";
        if(parameter.find(obj) == string::npos)
        {
            cout << "Second argument must be a filename of a triangular mesh .obj file." << endl;
            exit(EXIT_FAILURE);
        }
        parameter = argv[2];
        string pp = ".pp";
        if(parameter.find(pp) == string::npos)
        {
            cout << "Third argument must be a filename of the corresponding Meshlab point .pp file." << endl;
            exit(EXIT_FAILURE);
        }
        // store the filenames of the triangular mesh OBJect and Meshlab PickedPoints 
        global::obj_filename = argv[1];
        global::pp_filename = argv[2];
        // store the first argument in a string and check for print_csv and print_out parameters
        parameter = argv[0];
        if(parameter.find('s') != string::npos)
        {
            global::print_out = false;              //silence standard output to terminal
        }
        if(parameter.find('o') != string::npos)
        {
            global::print_oct = true;               //save FCC, EFD, and HMI data into GNU octave text data format
        }
        if(parameter.find('e') != string::npos)
        {
            global::print_exp = true;               //export polyline segments and height map images into corresponding files
        }
        if(parameter.find('p') != string::npos)
        {
            global::print_lmk = true;               //save updated Meshlab points into a new .pp file
        }
    }
}
// function for reading Meshlab pp files
vector<LMcoord> readMeshlabPoints(string pp_filename)
{
    // inform user
    if(global::print_out){
        cout << "Reading landmarks from " << pp_filename.c_str() <<" file..." << flush;
    }
    // initialize vector and read file
    vector<LMcoord> LMpoints;
    vector<LMcoord> LMpoints_ordered;
    ifstream inputFile(pp_filename.c_str());
    if(inputFile)
    {
        string line;
        int LMpoint_counter = 0;
        while (getline(inputFile, line))
        {
            size_t valid = line.find("<point");
            if(valid != string::npos)                                //check for point coordinates
            {
                float tmpx, tmpy, tmpz; int point;
                size_t name_index = line.find("name=");
                point = atof(line.c_str() + name_index + 6);
                size_t x_index = line.find("x=");
                tmpx = atof(line.c_str() + x_index + 3);
                size_t y_index = line.find("y=");
                tmpy = atof(line.c_str() + y_index + 3);
                size_t z_index = line.find("z=");
                tmpz = atof(line.c_str() + z_index + 3);
                LMcoord tempLMpoint = {point, tmpx, tmpy, tmpz};
                LMpoints.push_back(tempLMpoint);
                LMpoint_counter++;
            }
        }
        // list the landmarks in appropriate order as follows:
        // [0] nasion => name=1, [1] bregma => name=2, [2] opisthion => name=3, [3] mastoidale left => name=4
        // [4] mastoidale right => name=5, [5] orbitale left => name=6, [6] porion left => name=7, [7] porion right => name=8
        bool ordered = false;
        int landmark = 1;       //start with nasion
        while(!ordered)
        {
            bool found = false;
            for(int i = 0; i < LMpoint_counter; ++i)
            {
                if(LMpoints[i].name == landmark)
                {
                    found = true;
                    LMpoints_ordered.push_back(LMpoints[i]);
                }
            }
            if(found)
            {
                landmark++;
            }
            else
            {
                LMpoints_ordered.push_back({0, 0, 0, 0});
                landmark++;
            }
            if(landmark > 8)
            {
                ordered = true;
            }
        }
        // check if at least nasion, bregma, opisthion or nasion, opisthion, left orbitale and left+right porion landmarks are present in the file
        if(!((LMpoints_ordered[0].name == 1 && LMpoints_ordered[1].name == 2 && LMpoints_ordered[2].name == 3) || (LMpoints_ordered[0].name == 1 && LMpoints_ordered[2].name == 3 && LMpoints_ordered[5].name == 6 && LMpoints_ordered[6].name == 7 && LMpoints_ordered[7].name == 8)))
        {
            cout << "Not enough landmarks present in the MeshlabPoints file." << endl;
            cout << "You need at least 3 landmarks to initialize automatic registration." << endl;
            terminate();
        }
        inputFile.close();
        if(global::print_out){
            cout << " OK" << endl;
        }
        return LMpoints_ordered;
    }
    else
    {
        cout << "Cannot read file or file not present." << endl;
        terminate();
    }
}
// function for exporting 3D points to Meshlab pp files.
void writeMeshlabPoints(string pp_filename, vector<LMcoord> LMpoints)
{
    ofstream outputFile(pp_filename.c_str());
    // check if writing to file is permitted
    if (!outputFile.is_open())
    {
      cout << "Error opening " << pp_filename.c_str() << "for write" << endl;
      terminate();
    }
    else
    {
        // get current local time and date
        time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
        // inform user about writing process
        if(global::print_out){
            cout << "Writing to " << pp_filename.c_str() << " file..." << flush;
        }
        // writing header to file
        outputFile << "<!DOCTYPE PickedPoints>\n<PickedPoints>\n <DocumentData>\n";
        outputFile << "  <DateTime time=\"" << put_time(localtime(&now), "%T") << "\" date=\"" << put_time(localtime(&now), "%F") << "\"/>\n";
        outputFile << "  <User name=\"skullanalyzer\"/>\n";
        outputFile << "  <DataFileName name=\"" << global::obj_filename.c_str() << "\"/>\n";
        outputFile << "  <templateName name=\"\"/>\n </DocumentData>\n";
        // write vertices to file
        for(vector<LMcoord>::iterator lm_it = LMpoints.begin(); lm_it != LMpoints.end(); ++lm_it)
        {
            if(lm_it->name != 0)
            {
                outputFile << " <point active=\"1\" name=\"" << lm_it->name <<"\" x=\"" << lm_it->x <<"\" y=\"" <<lm_it->y <<"\" z=\"" <<lm_it->z <<"\"/>\n";
            }
        }
        outputFile << "</PickedPoints>";
        outputFile.close();
        if(global::print_out){
            cout << " OK" << endl;
        }
    }
}
// function for loading triangular meshes from obj files
vector<Mesh> readObjMesh(string obj_filename)
{
    // inform user
    if(global::print_out){
        cout << "Loading Vertices and Faces from " << obj_filename.c_str() <<" file..." << flush;
    }
    // initialize vectors
    vector<VCoord> vertex;
    vector<Faces> face;
    vector<Mesh> MeshElements;
    // initialize counters for vertices and faces
    int vertex_counter = 0;
    int face_counter = 0;
    ifstream inputFile(obj_filename.c_str());
    if(inputFile)
    {
        // define string for reading obj file per line
        string line;
        while (getline(inputFile, line))
        {
            if(line[0] == 'v' && line[1] == ' ')                                //check for vertex coordinates
            {
                float tmpx, tmpy, tmpz;
                sscanf(line.c_str(), "v %f %f %f" ,&tmpx,&tmpy,&tmpz);
                VCoord temp3D = {tmpx, tmpy, tmpz};
                vertex.push_back(temp3D);
                vertex_counter++;
            }
            if(line[0] == 'f' && line[1] == ' ')                                //check for faces
            {
                // check for non triangular mesh by ensuring that only
                // triplets of vertices, normals and texture coordinates
                // are referenced in each face.
                int v1=0, v2=0, v3=0, v4=0, vt1=0, vt2=0, vt3=0, vt4=0, vn1=0, vn2=0, vn3=0, vn4=0;
                sscanf(line.c_str(), "f %d %d %d %d" ,&v1,&v2,&v3,&v4);
                if (v1 > 0 && v2 > 0 && v3 > 0 && v4 > 0)
                {
                    cout << "Mesh is not triangular." << endl;
                    terminate();
                }
                sscanf(line.c_str(), "f %d/%d %d/%d %d/%d %d/%d" ,&v1,&vt1,&v2,&vt2,&v3,&vt3,&v4,&vt4);
                if (v1 > 0 && v2 > 0 && v3 > 0 && v4 > 0 && vt1 > 0 && vt2 > 0 && vt3 > 0 && vt4 > 0)
                {
                    cout << "Mesh is not triangular." << endl;
                    terminate();
                }
                sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", &v1,&vt1,&vn1,&v2,&vt2,&vn2,&v3,&vt3,&vn3,&v4,&vt4,&vn4);
                if (v1 > 0 && v2 > 0 && v3 > 0 && v4 > 0 && vt1 > 0 && vt2 > 0
                    && vt3 > 0 && vt4 > 0 && vn1 > 0 && vn2 > 0 && vn3 > 0 && vn4 > 0)
                {
                    cout << "Mesh is not triangular." << endl;
                    terminate();
                }
                sscanf(line.c_str(), "f %d//%d %d//%d %d//%d %d//%d", &v1,&vn1,&v2,&vn2,&v3,&vn3,&v4,&vn4);
                if (v1 > 0 && v2 > 0 && v3 > 0 && v4 > 0 && vn1 > 0 && vn2 > 0 && vn3 > 0 && vn4 > 0)
                {
                    cout << "Mesh is not triangular." << endl;
                    terminate();
                }
                // increment face counter
                face_counter++;
                // Reset temp face variables
                v1=0;v2=0;v3=0,vt1=0,vt2=0,vt3=0,vn1=0,vn2=0,vn3=0;

                // scan for faces registering only vertices
                sscanf(line.c_str(), "f %d %d %d" ,&v1,&v2,&v3);
                if (v1 > 0 && v2 > 0 && v3 > 0)
                {
                    Faces temp_face_v = {v1, v2, v3};
                    face.push_back(temp_face_v);
                }
                else
                {
                    // scan for faces registering vertices and texture coordinates
                    sscanf(line.c_str(), "f %d/%d %d/%d %d/%d", &v1,&vt1,&v2,&vt2,&v3,&vt3);
                    if (v1 > 0 && v2 > 0 && v3 > 0 && vt1 > 0 && vt2 > 0 && vt3 > 0)
                    {
                        Faces temp_face_v = {v1, v2, v3};
                        face.push_back(temp_face_v);
                    }
                    else
                    {
                        // scan for faces registering vertices, texture coordinates and normals
                        sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d", &v1,&vt1,&vn1,&v2,&vt2,&vn2,&v3,&vt3,&vn3);
                        if (v1 > 0 && v2 > 0 && v3 > 0 && vt1 > 0 && vt2 > 0 && vt3 > 0 && vn1 > 0 && vn2 > 0 && vn3 > 0)
                        {
                            Faces temp_face_v = {v1, v2, v3};
                            face.push_back(temp_face_v);
                        }
                        else
                        {
                            // scan for faces registering vertices and normals
                            sscanf(line.c_str(), "f %d//%d %d//%d %d//%d", &v1,&vn1,&v2,&vn2,&v3,&vn3);
                            if (v1 > 0 && v2 > 0 && v3 > 0 && vn1 > 0 && vn2 > 0 && vn3 > 0)
                            {
                                Faces temp_face_v = {v1, v2, v3};
                                face.push_back(temp_face_v);
                            }
                            else
                            {
                                cout << "Mesh is not triangular." << endl;
                                terminate();
                            }
                        }
                    }
                }
            }
        }
        inputFile.close();
        if(global::print_out){
            cout << " OK" << endl;
            cout << "Mesh contains " << vertex_counter << " vertices and " << face_counter << " faces." << endl;
        }
        // populate return structure with face indices and vertex coordinates
        int f_V1, f_V2, f_V3; double V1x, V1y, V1z, V2x, V2y, V2z, V3x, V3y, V3z; Mesh tmpmesh;
        for(int i = 0; i < face_counter; ++i)
        {            
            f_V1 = face[i].a; f_V2 = face[i].b; f_V3 = face[i].c;
            V1x = vertex[f_V1 - 1].x; V1y = vertex[f_V1 - 1].y; V1z = vertex[f_V1 - 1].z;
            V2x = vertex[f_V2 - 1].x; V2y = vertex[f_V2 - 1].y; V2z = vertex[f_V2 - 1].z;
            V3x = vertex[f_V3 - 1].x; V3y = vertex[f_V3 - 1].y; V3z = vertex[f_V3 - 1].z;
            tmpmesh = {V1x, V1y, V1z, V2x, V2y, V2z, V3x, V3y, V3z};
            MeshElements.push_back(tmpmesh);
        }
        return MeshElements;
    }
    else
    {
        cout << "Cannot read file or file not present." << endl;
        terminate();
    }
}
// function for calculating DotProduct between two vectors
double dotProduct(VCoord A, VCoord B)
{
    double dotprod = A.x * B.x + A.y * B.y + A.z * B.z;
    return dotprod;
}
// function for calculating CrossProduct between two vector
VCoord crossProduct(VCoord A, VCoord B)
{
    double x, y, z;
    x = A.y * B.z - A.z * B.y;
    y = A.z * B.x - A.x * B.z;
    z = A.x * B.y - A.y * B.x;
    VCoord crossprod = {x, y, z};
    return crossprod;
}
// function for normalizing a vector
VCoord normalizeVector(VCoord A)
{
    double v_len;
    v_len = sqrt(A.x * A.x + A.y * A.y + A.z * A.z);
    VCoord normal = {A.x / v_len, A.y / v_len, A.z / v_len};
    return normal;
}
// function for calculating distance between two 3D points
double distancePoints3D(VCoord A, VCoord B)
{
    double distance;
    double x, y, z;
    x = (A.x - B.x) * (A.x - B.x);
    y = (A.y - B.y) * (A.y - B.y);
    z = (A.z - B.z) * (A.z - B.z);
    distance = sqrt(x + y + z);
    return distance;
}
// function for calculating the midpoint between two 3D points
VCoord midPoint3D(VCoord A, VCoord B)
{
    VCoord midPoint;
    double x, y, z;
    x = (A.x + B.x) / 2;
    y = (A.y + B.y) / 2;
    z = (A.z + B.z) / 2;
    midPoint = {x, y, z};
    return midPoint;
}
// function for calculating distance between two 2D points
double distancePoints2D(PCoord A, PCoord B)
{
    double distance;
    double x, y;
    x = (A.x - B.x) * (A.x - B.x);
    y = (A.y - B.y) * (A.y - B.y);
    distance = sqrt(x + y);
    return distance;
}
// function for defining anatomical planes according to present landmarks. Nominally the
// transverse place normal faces upwards (superiorly) and the coronal plane normal faces
// to the front (anteriorly) and the sagittal plane normal towards the left hand side.
vector<VCoord> anatomicalPlanes(vector<LMcoord> LMpoints)
{
    // initialize vectors and variables
    vector<VCoord> anatomicalElements;
    VCoord nasion;
    VCoord bregma;
    VCoord opisthion;
    VCoord mastoidaleL;
    VCoord mastoidaleR;
    VCoord orbitale;
    VCoord porionL;
    VCoord porionR;
    int LM_Fra = 0;
    for(int i = 0; i < LMpoints.size(); ++i)
    {
        if(LMpoints[i].name == 1)
        {
            nasion = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            global::nasion3D = nasion;
            LM_Fra++;
        }
        if(LMpoints[i].name == 2)
        {
            bregma = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            global::bregma3D = bregma;
            landmark::bregma = true;
        }
        if(LMpoints[i].name == 3)
        {
            opisthion = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            global::opisthion3D = opisthion;
            LM_Fra++;
        }
        if(LMpoints[i].name == 4)
        {
            mastoidaleL = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            landmark::mastoidaleL = true;
        }
        if(LMpoints[i].name == 5)
        {
            mastoidaleR = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            landmark::mastoidaleR = true;
        }
        if(LMpoints[i].name == 6)
        {
            orbitale = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            LM_Fra++;
        }
        if(LMpoints[i].name == 7)
        {
            porionL = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            LM_Fra++;
        }
        if(LMpoints[i].name == 8)
        {
            porionR = {LMpoints[i].x, LMpoints[i].y, LMpoints[i].z};
            LM_Fra++;
        }
    }
    // check that minimum requirements are met
    if(LM_Fra < 5)
    {
        cout << "Not enough landmarks to calculale cranial orientation" << endl;
        terminate();
    }
    if(LM_Fra == 5)                         //use nasion, opisthion, orbitale, left and right porion
    {
        landmark::Frankfurt = true;
        if(global::print_out){
            cout << "Defining anatomical position by Frankfurt plane." << endl;
        }
        // get orbitale - left porion and orbitale - right porion vectors
        VCoord orb_pL = {porionL.x - orbitale.x, porionL.y - orbitale.y, porionL.z - orbitale.z};
        VCoord orb_pR = {porionR.x - orbitale.x, porionR.y - orbitale.y, porionR.z - orbitale.z};
        // calculale their crossproduct and normalize it
        VCoord crossprod = crossProduct(orb_pL, orb_pR);
        VCoord normal_T = normalizeVector(crossprod);             //transverse plane normal
        // get nasion-opisthion vector
        VCoord n_o = {opisthion.x - nasion.x, opisthion.y - nasion.y, opisthion.z - nasion.z};
        // calculale crossproduct with transverse normal
        crossprod = crossProduct(n_o, normal_T);
        VCoord normal_S = normalizeVector(crossprod);             //sagittal plane normal
        // calculale crossproduct with sagittal normal
        crossprod = crossProduct(normal_S, normal_T);
        VCoord normal_C = normalizeVector(crossprod);             //coronal plane normal
        //
        // apend anatomical plane normals to anatomicalElements vector
        anatomicalElements.push_back(normal_S);         //sagittal plane normal     [0]
        anatomicalElements.push_back(normal_T);         //transverse plane normal   [1]
        anatomicalElements.push_back(normal_C);         //coronal plane normal      [2]
        //
        // apend available landmark coordinates to anatomicalElements vector
        anatomicalElements.push_back(nasion);           // [3]
        if(landmark::bregma){
            anatomicalElements.push_back(bregma);       // [4]
        }
        else{
            anatomicalElements.push_back({0, 0, 0});    //set to zero to align index [4]
        }
        anatomicalElements.push_back(opisthion);        // [5]
        if(landmark::mastoidaleL){
            anatomicalElements.push_back(mastoidaleL);  // [6]
        }
        else{
            anatomicalElements.push_back({0, 0, 0});    //set to zero to align index [6]
        }
        if(landmark::mastoidaleR){
            anatomicalElements.push_back(mastoidaleR);  // [7]
        }
        else{
            anatomicalElements.push_back({0, 0, 0});    //set to zero to align index [7]
        }
        anatomicalElements.push_back(orbitale);         // [8]
        anatomicalElements.push_back(porionL);          // [9]
        anatomicalElements.push_back(porionR);          // [10]
    }
    return anatomicalElements;
}
// function for specifying landmark along with midsagittal plane normal: reminder that if
// coronal plane is chosen, then it is inverted so that the subsequent 2D projections are
// positioned appropriately
vector<VCoord> chooseSlicingPlane(vector<VCoord> anatomicalElements, string landmark, string plane)
{
    // initialize vectors and variables
    vector<VCoord> intersectionPlane;
    int LM;
    VCoord normal;
    if(plane == "sagittal")
    {
        if(global::print_out){
            cout << "Slicing along sagittal plane " << flush;
        }
        normal = {anatomicalElements[0].x, anatomicalElements[0].y, anatomicalElements[0].z};
    }
    if(plane == "transverse")
    {
        if(global::print_out){
            cout << "Slicing along transverse plane " << flush;
        }
        normal = {anatomicalElements[1].x, anatomicalElements[1].y, anatomicalElements[1].z};
    }
    if(plane == "coronal")
    {
        if(global::print_out){
            cout << "Slicing along coronal plane " << flush;
        }
        normal = {-anatomicalElements[2].x, -anatomicalElements[2].y, -anatomicalElements[2].z};
    }
    if(landmark == "nasion")
    {
        LM = 3;
        if(global::print_out){
            cout << "at nasion:" << endl;
        }
    }
    if(landmark == "bregma")
    {
        LM = 4;
        if(global::print_out){
            cout << "at bregma:" << endl;
        }
    }
    if(landmark == "mastoidale left")
    {
        LM = 6;
        if(global::print_out){
            cout << "at left mastoidale:" << endl;
        }
    }
    if(landmark == "mastoidale right")
    {
        LM = 7;
        if(global::print_out){
            cout << "at right mastoidale:" << endl;
        }
    }
    VCoord point = {anatomicalElements[LM].x, anatomicalElements[LM].y, anatomicalElements[LM].z};
    intersectionPlane.push_back(normal);
    intersectionPlane.push_back(point);
    return(intersectionPlane);
}
// function for extracting intersection points ginen a plane specified by its normal and a point
// three equally spaced intersection points per face are calculated to minimize the effect of subsequent rasterization
vector<VCoord> sliceMesh(vector<Mesh> MeshElements, vector<VCoord> intersectionPlane)
{
    // inform user
    if(global::print_out){
        cout << "Calculating intersection points..." << flush;
    }
    // initialize vectors and variables
    int iFaces = 0;
    vector<VCoord> iPoints;
    VCoord V1, V2, V3;
    double dotV1, dotV2, dotV3;
    VCoord normal;
    VCoord point;
    normal = {intersectionPlane[0].x, intersectionPlane[0].y, intersectionPlane[0].z};
    point = {intersectionPlane[1].x, intersectionPlane[1].y, intersectionPlane[1].z};
    // for each triplet of vertices find which face is intersected by specified plane
    for(int i = 0; i < MeshElements.size(); ++i)
    {
        V1 = {MeshElements[i].V1x - point.x, MeshElements[i].V1y - point.y, MeshElements[i].V1z - point.z};
        V2 = {MeshElements[i].V2x - point.x, MeshElements[i].V2y - point.y, MeshElements[i].V2z - point.z};
        V3 = {MeshElements[i].V3x - point.x, MeshElements[i].V3y - point.y, MeshElements[i].V3z - point.z};
        dotV1 = dotProduct(V1, normal);
        dotV2 = dotProduct(V2, normal);
        dotV3 = dotProduct(V3, normal);
        if(!((dotV1 < 0 && dotV2 < 0 && dotV3 < 0) || (dotV1 > 0 && dotV2 > 0 && dotV3 > 0)))   // true if intersecting
        {
            VCoord CSv1; int v1 = 0;
            VCoord CSv2; int v2 = 0;
            VCoord CSv3; int v3 = 0;
            if(dotV1 * dotV2 < 0)       // vertices 1 and 2 of ith face lie on opposite sides
            {
                VCoord Pv1 = {MeshElements[i].V1x - dotV1 * normal.x, MeshElements[i].V1y - dotV1 * normal.y, MeshElements[i].V1z - dotV1 * normal.z};
                VCoord Pv2 = {MeshElements[i].V2x - dotV2 * normal.x, MeshElements[i].V2y - dotV2 * normal.y, MeshElements[i].V2z - dotV2 * normal.z};
                double Dv1 = abs(dotV1);
                double Dv2 = abs(dotV2);
                double sum = Dv1 + Dv2;
                CSv1 = {(Dv1 * Pv1.x + Dv2 * Pv2.x)/sum, (Dv1 * Pv1.y + Dv2 * Pv2.y)/sum, (Dv1 * Pv1.z + Dv2 * Pv2.z)/sum};
                v1++;
            }
            if(dotV1 * dotV2 == 0)      // vertices 1 and 2 of ith face lie on the cross-section plane
            {
                CSv1 = {MeshElements[i].V1x, MeshElements[i].V1y, MeshElements[i].V1z};
                CSv2 = {MeshElements[i].V2x, MeshElements[i].V2y, MeshElements[i].V2z};
                v1++;
                v2++;
            }
            if(dotV1 * dotV3 < 0)       // vertices 1 and 3 of ith face lie on opposite sides
            {
                VCoord Pv1 = {MeshElements[i].V1x - dotV1 * normal.x, MeshElements[i].V1y - dotV1 * normal.y, MeshElements[i].V1z - dotV1 * normal.z};
                VCoord Pv3 = {MeshElements[i].V3x - dotV3 * normal.x, MeshElements[i].V3y - dotV3 * normal.y, MeshElements[i].V3z - dotV3 * normal.z};
                double Dv1 = abs(dotV1);
                double Dv3 = abs(dotV3);
                double sum = Dv1 + Dv3;
                CSv2 = {(Dv1 * Pv1.x + Dv3 * Pv3.x)/sum, (Dv1 * Pv1.y + Dv3 * Pv3.y)/sum, (Dv1 * Pv1.z + Dv3 * Pv3.z)/sum};
                v2++;
            }
            if(dotV1 * dotV3 == 0)      // vertices 1 and 3 of ith face lie on the cross-section plane
            {
                CSv1 = {MeshElements[i].V1x, MeshElements[i].V1y, MeshElements[i].V1z};
                CSv3 = {MeshElements[i].V3x, MeshElements[i].V3y, MeshElements[i].V3z};
                v1++;
                v3++;
            }
            if(dotV2 * dotV3 < 0)       // vertices 2 and 3 of ith face lie on opposite sides
            {
                VCoord Pv2 = {MeshElements[i].V2x - dotV2 * normal.x, MeshElements[i].V2y - dotV2 * normal.y, MeshElements[i].V2z - dotV2 * normal.z};
                VCoord Pv3 = {MeshElements[i].V3x - dotV3 * normal.x, MeshElements[i].V3y - dotV3 * normal.y, MeshElements[i].V3z - dotV3 * normal.z};
                double Dv2 = abs(dotV2);
                double Dv3 = abs(dotV3);
                double sum = Dv2 + Dv3;
                CSv3 = {(Dv2 * Pv2.x + Dv3 * Pv3.x)/sum, (Dv2 * Pv2.y + Dv3 * Pv3.y)/sum, (Dv2 * Pv2.z + Dv3 * Pv3.z)/sum};
                v3++;
            }
            if(dotV2 * dotV3 == 0)      // vertices 1 and 3 of ith face lie on the cross-section plane
            {
                CSv2 = {MeshElements[i].V2x, MeshElements[i].V2y, MeshElements[i].V2z};
                CSv3 = {MeshElements[i].V3x, MeshElements[i].V3y, MeshElements[i].V3z};
                v2++;
                v3++;
            }
            if(v1 == 1 && v2 == 1)
            {
                VCoord tmpuPoint;
                tmpuPoint = {(9 * CSv1.x + CSv2.x)/10, (9 * CSv1.y + CSv2.y)/10, (9 * CSv1.z + CSv2.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(7 * CSv1.x + 3 * CSv2.x)/10, (7 * CSv1.y + 3 * CSv2.y)/10, (7 * CSv1.z + 3 * CSv2.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(CSv1.x + CSv2.x)/2, (CSv1.y + CSv2.y)/2, (CSv1.z + CSv2.z)/2};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(3 * CSv1.x + 7 * CSv2.x)/10, (3 * CSv1.y + 7 * CSv2.y)/10, (3 * CSv1.z + 7 * CSv2.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(CSv1.x + 9 * CSv2.x)/10, (CSv1.y + 9 * CSv2.y)/10, (CSv1.z + 9 * CSv2.z)/10};
                iPoints.push_back(tmpuPoint);
            }
            else if(v1 == 1 && v3 == 1)
            {
                VCoord tmpuPoint;
                tmpuPoint = {(9 * CSv1.x + CSv3.x)/10, (9 * CSv1.y + CSv3.y)/10, (9 * CSv1.z + CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(7 * CSv1.x + 3 * CSv3.x)/10, (7 * CSv1.y + 3 * CSv3.y)/10, (7 * CSv1.z + 3 * CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(CSv1.x + CSv3.x)/2, (CSv1.y + CSv3.y)/2, (CSv1.z + CSv3.z)/2};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(3 * CSv1.x + 7 * CSv3.x)/10, (3 * CSv1.y + 7 * CSv3.y)/10, (3 * CSv1.z + 7 * CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(CSv1.x + 9 * CSv3.x)/10, (CSv1.y + 9 * CSv3.y)/10, (CSv1.z + 9 * CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
            }
            else
            {
                VCoord tmpuPoint;
                tmpuPoint = {(9 * CSv2.x + CSv3.x)/10, (9 * CSv2.y + CSv3.y)/10, (9 * CSv2.z + CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(7 * CSv2.x + 3 * CSv3.x)/10, (7 * CSv2.y + 3 * CSv3.y)/10, (7 * CSv2.z + 3 * CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(CSv2.x + CSv3.x)/2, (CSv2.y + CSv3.y)/2, (CSv2.z + CSv3.z)/2};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(3 * CSv2.x + 7 * CSv3.x)/10, (3 * CSv2.y + 7 * CSv3.y)/10, (3 * CSv2.z + 7 * CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
                tmpuPoint = {(CSv2.x + 9 * CSv3.x)/10, (CSv2.y + 9 * CSv3.y)/10, (CSv2.z + 9 * CSv3.z)/10};
                iPoints.push_back(tmpuPoint);
            }
            iFaces++;
        }
    }
    if(global::print_out){
        cout << " OK" << endl;
        cout << "   " << iFaces << " faces were intersected." << endl;
    }
    return iPoints;
}
// function for finding nasion-opisthocranion max distance and opisthocranion 3D coordinates
void horMaxDistanceNasionOpisthocranion(vector<VCoord> MidSagPoints, vector<VCoord> anatomicalElements)
{
    // define variable
    double NasionOpCranionLength;
    VCoord normal_T = {anatomicalElements[1].x, anatomicalElements[1].y, anatomicalElements[1].z};
    VCoord nasion = {anatomicalElements[3].x, anatomicalElements[3].y, anatomicalElements[3].z};
    VCoord opisthion = {anatomicalElements[5].x, anatomicalElements[5].y, anatomicalElements[5].z};
    vector<VCoord> upper_points;
    vector<VCoord> lower_points;
    // calculale two threshold points +-1mm appart from nasion landmark pointing upwards and downwards with respect to transverse plane
    VCoord Upper_Threshold = {nasion.x + normal_T.x, nasion.y + normal_T.y, nasion.z + normal_T.z};
    VCoord Lower_Threshold = {nasion.x - normal_T.x, nasion.y - normal_T.y, nasion.z - normal_T.z};
    // counters for points lying between the threshold distance above and below the transverse plane passing through the nasion landmark
    int u_p = 0;
    int l_p = 0;
    // calculale minimum length from nasion landmark as the distance between nasion - opisthion
    double nasion_opisthion_dist = distancePoints3D(nasion, opisthion);
    // iterate once throught the entire midsaggital curve
    for(int i = 0; i < MidSagPoints.size(); ++i)
    {
        // calculale dot product of each point to the upper and lower threshold points
        VCoord VecUT = {MidSagPoints[i].x - Upper_Threshold.x, MidSagPoints[i].y - Upper_Threshold.y, MidSagPoints[i].z - Upper_Threshold.z};
        double dotVecUT = dotProduct(VecUT, normal_T);
        VCoord VecLT = {MidSagPoints[i].x - Lower_Threshold.x, MidSagPoints[i].y - Lower_Threshold.y, MidSagPoints[i].z - Lower_Threshold.z};
        double dotVecLT = dotProduct(VecLT, normal_T);
        // calculale distance of point to nasion
        double nasion_point_dist = distancePoints3D(nasion, MidSagPoints[i]);
        //keeping what lies between the upper and lower threshold planes and is more distant than nasion-opisthion
        if(dotVecUT < 0 && dotVecLT > 0 && nasion_point_dist > nasion_opisthion_dist)
        {
            // check if point lies above or below nasion transverse plane
            VCoord Vec = {MidSagPoints[i].x - nasion.x, MidSagPoints[i].y - nasion.y, MidSagPoints[i].z - nasion.z};
            double dotVec = dotProduct(Vec, normal_T);
            if(dotVec > 0)
            {
                u_p++;
                //upper_dist = nasion_point_dist; cout << upper_dist << endl;
                upper_points.push_back(MidSagPoints[i]);
            }
            if(dotVec < 0)
            {
                l_p++;
                //lower_dist = nasion_point_dist; cout << lower_dist << endl;
                lower_points.push_back(MidSagPoints[i]);
            }
        }
    }
    // calculate the intersection point between upper and lower points closer to the transverse plane passing through nasion
    // at the back fo the skull. Check the number of points found above and below the transverse plane. If more than one point
    // is found then keep the closest to the tranverse plane passing through nasion, if only one is present use that one, if
    // no points were found keep the point closest to transverse plane from the opposite side
    if(u_p > 0 && l_p > 0)              //both upper and lower points are present
    {
        double min_upper_dist = 1;
        VCoord u_point;
        for(int i = 0; i < u_p; ++i)    //for upper points
        {
            VCoord v = {upper_points[i].x - nasion.x, upper_points[i].y - nasion.y, upper_points[i].z - nasion.z};
            double dist = abs(dotProduct(v, normal_T));
            if(dist <= min_upper_dist)  //only if closer
            {
                min_upper_dist = dist;
                u_point = upper_points[i];
            }
        }
        double min_lower_dist = 1;
        VCoord l_point;
        for(int i = 0; i < l_p; ++i)    //for lower points
        {
            VCoord v = {lower_points[i].x - nasion.x, lower_points[i].y - nasion.y, lower_points[i].z - nasion.z};
            double dist = abs(dotProduct(v, normal_T));
            if(dist <= min_lower_dist)   //only if closer
            {
                min_lower_dist = dist;
                l_point = lower_points[i];
            }
        }
        // find the orthogonal projections of upper and lower points onto the transverse plane passing through nasion
        // and calculate the point of intersection based on the ratio of their disctances from the transverse plane
        VCoord u_v = {u_point.x - nasion.x, u_point.y - nasion.y, u_point.z - nasion.z};
        double dot_uvT = dotProduct(u_v, normal_T);
        double uvT = abs(dot_uvT);
        VCoord u_proj = {u_point.x - dot_uvT * normal_T.x, u_point.y - dot_uvT * normal_T.y, u_point.z - dot_uvT * normal_T.z};
        
        VCoord l_v = {l_point.x - nasion.x, l_point.y - nasion.y, l_point.z - nasion.z};
        double dot_lvT = dotProduct(l_v, normal_T);
        double lvT = abs(dot_lvT);
        VCoord l_proj = {l_point.x - dot_lvT * normal_T.x, l_point.y - dot_lvT * normal_T.y, l_point.z - dot_lvT * normal_T.z};
                
        double sum = uvT + lvT;
        global::opisthocranion3D = {(uvT * u_proj.x + lvT * l_proj.x)/sum, (uvT * u_proj.y + lvT * l_proj.y)/sum, (uvT * u_proj.z + lvT * l_proj.z)/sum};
        // calculate the distance between nasion and opisthocranion
        NasionOpCranionLength = distancePoints3D(nasion, global::opisthocranion3D);
    }
    if(u_p > 0 && l_p == 0)              //only upper points are present
    {
        double min_upper_dist = 1;
        VCoord u_point;
        for(int i = 0; i < u_p; ++i)    //for upper points
        {
            VCoord v = {upper_points[i].x - nasion.x, upper_points[i].y - nasion.y, upper_points[i].z - nasion.z};
            double dist = abs(dotProduct(v, normal_T));
            if(dist <= min_upper_dist)   //only if closer
            {
                min_upper_dist = dist;
                u_point = upper_points[i];
            }
        }
        global::opisthocranion3D = u_point;
        NasionOpCranionLength = distancePoints3D(nasion, u_point);
    }
    if(u_p == 0 && l_p > 0)              //only lower points are present
    {
        double min_lower_dist = 1;
        VCoord l_point;
        for(int i = 0; i < l_p; ++i)    //for lower points
        {
            VCoord v = {lower_points[i].x - nasion.x, lower_points[i].y - nasion.y, lower_points[i].z - nasion.z};
            double dist = abs(dotProduct(v, normal_T));
            if(dist <= min_lower_dist)   //only if closer
            {
                min_lower_dist = dist;
                l_point = lower_points[i];
            }
        }
        global::opisthocranion3D = l_point;
        NasionOpCranionLength = distancePoints3D(nasion, l_point);
    }

    if(global::print_out){
        cout << "Nasion - pseudoOpisthocranion horizontal length is " << NasionOpCranionLength << endl;
    }
}
// function for optimizing the location of nasion landmark
void optimizeNasion(vector<Contour> boundaryPixels, vector<PCoord> MidSagSection, vector<VCoord> intersectionPlane, VCoord X2D_vector)
{
    // calculate local x, y axis and origin to transform the optimized point from the 2D midsagittal section back to original 3D space
    // define local y axis
    VCoord SP_normal = {intersectionPlane[0].x, intersectionPlane[0].y, intersectionPlane[0].z};    //normal
    VCoord XY_origin = {intersectionPlane[1].x, intersectionPlane[1].y, intersectionPlane[1].z};    //point
    VCoord Y2D_vector = crossProduct(X2D_vector, SP_normal);
    // normalize local axes
    X2D_vector = normalizeVector(X2D_vector);
    Y2D_vector = normalizeVector(Y2D_vector);
    if(global::print_out)
    {
        cout << "Optimizing landmark coordinates for nasion" << endl;
        cout << "User coordinates are: x=" << XY_origin.x << " y=" << XY_origin.y << " z=" << XY_origin.z << endl;
    }
    // search the corresponding points of the boundary pixels near the original nasion location (beginning of vector) and pick the point
    // with the most negative value on the local 2D x-axis. The target points are expected to be on the same column or to the right (<col)
    // of the initial pixel. It is assumed that the user defined nasion landmark if located within 5mm of its optimal position
    //
    // set optimal point to current nasion location
    PCoord optimal_point = global::nasion2D;
    // start searching towards the lower side of the nasion
    bool lower_found = false;
    int it = boundaryPixels.size() - 1;     //start from penultimate point in vector to compare with ultimate which coinsidess with starting point
    while(!lower_found)
    {
        Contour pixel = boundaryPixels[it];
        if(pixel.col <= boundaryPixels[it + 1].col && it > boundaryPixels.size() - 11)
        {
            for(int p_it = 0; p_it < pixel.pointer.size() - 1; ++p_it)
            {
                if(MidSagSection[pixel.pointer[p_it]].x < optimal_point.x)
                {
                    optimal_point = MidSagSection[pixel.pointer[p_it]];
                }
            }
        }
        else
        {
            lower_found = true;
        }
        --it;
    }
    // if lower optimal was NOT found start searching towards the upper side of the nasion
    // make sure that the most posterior point found is a local minimum, which means that subsequent points move anteriorly
    // in opposite case (all pixels move posteriorly, i.e. check_index increases to 10) use user defined nasion as optimal
    int check_index = 0;
    if(optimal_point.x == global::nasion2D.x && optimal_point.y == global::nasion2D.y)
    {
        bool upper_found = false;
        it = 1;                         //start from second point to compare with previous point
        while(!upper_found)
        {
            Contour pixel = boundaryPixels[it];
            if(pixel.col <= boundaryPixels[it - 1].col && it < 11)      //until ten pixels are considered
            {
                for(int p_it = 0; p_it < pixel.pointer.size() - 1; ++p_it)
                {
                    if(MidSagSection[pixel.pointer[p_it]].x < optimal_point.x)
                    {
                        optimal_point = MidSagSection[pixel.pointer[p_it]];
                        // if lower point found 
                    }
                }
                // update the check index for every pixel lying at the same vertical axis or posteriorly of the previous pixel
                check_index++;
            }
            else
            {
                upper_found = true;
            }
            ++it;
        }
    }
    // otherwise set landmark::nasion_optimal to true
    else
    {
        landmark::nasion_optimal = true;
    }
    // if all superior pixels were posteriorly positioned then use the initial nasion2D user coordinates as optimal location
    if(check_index == 10)
    {
        optimal_point = global::nasion2D;
        landmark::nasion_optimal = false;
    }
    // otherwise set landmark::nasion_optimal to true
    else
    {
        landmark::nasion_optimal = true;
    }
    // convert the optimal 2D coordinates into original 3D coordinates for the optimized nasion landmark
    double opt_x = XY_origin.x + optimal_point.x * X2D_vector.x + optimal_point.y * Y2D_vector.x;
    double opt_y = XY_origin.y + optimal_point.x * X2D_vector.y + optimal_point.y * Y2D_vector.y;
    double opt_z = XY_origin.z + optimal_point.x * X2D_vector.z + optimal_point.y * Y2D_vector.z;
    if(global::print_out)
    {
        cout << "Optimal coordinates are: x=" <<opt_x << " y=" << opt_y << " z=" << opt_z << endl;
    }
    global::nasion3D = {opt_x, opt_y, opt_z};
}
// function for projecting 3D points onto slicing plane 2D coordinate system
vector<PCoord> project3Dto2Dplane(vector<VCoord> SectionPoints, vector<VCoord> intersectionPlane, VCoord X2D_vector)
{
    vector<PCoord> crossSection;
    // define local y axis
    VCoord SP_normal = {intersectionPlane[0].x, intersectionPlane[0].y, intersectionPlane[0].z};    //normal
    VCoord XY_origin = {intersectionPlane[1].x, intersectionPlane[1].y, intersectionPlane[1].z};    //point
    VCoord Y2D_vector = crossProduct(X2D_vector, SP_normal);
    // normalize local axes
    X2D_vector = normalizeVector(X2D_vector);
    Y2D_vector = normalizeVector(Y2D_vector);
    // transform cross-sectional 3D points into planar local 2D coordinates
    for(int i = 0; i < SectionPoints.size(); ++i)
    {
        VCoord XY = {SectionPoints[i].x - XY_origin.x, SectionPoints[i].y - XY_origin.y, SectionPoints[i].z - XY_origin.z};
        double x = XY.x * X2D_vector.x + XY.y * X2D_vector.y + XY.z * X2D_vector.z;
        double y = XY.x * Y2D_vector.x + XY.y * Y2D_vector.y + XY.z * Y2D_vector.z;
        PCoord point2D = {x, y};
        crossSection.push_back(point2D);
    }
    // transform nasion, opisthion and opisthocranion landmarks to specified 2D local coordinate system
    VCoord XY_nasion = {global::nasion3D.x - XY_origin.x, global::nasion3D.y - XY_origin.y, global::nasion3D.z - XY_origin.z};
    double x_nas = XY_nasion.x * X2D_vector.x + XY_nasion.y * X2D_vector.y + XY_nasion.z * X2D_vector.z;
    double y_nas = XY_nasion.x * Y2D_vector.x + XY_nasion.y * Y2D_vector.y + XY_nasion.z * Y2D_vector.z;
    global::nasion2D = {x_nas, y_nas};
    VCoord XY_opisthion = {global::opisthion3D.x - XY_origin.x, global::opisthion3D.y - XY_origin.y, global::opisthion3D.z - XY_origin.z};
    double x_op = XY_opisthion.x * X2D_vector.x + XY_opisthion.y * X2D_vector.y + XY_opisthion.z * X2D_vector.z;
    double y_op = XY_opisthion.x * Y2D_vector.x + XY_opisthion.y * Y2D_vector.y + XY_opisthion.z * Y2D_vector.z;
    global::opisthion2D = {x_op, y_op};
    VCoord XY_opcranion = {global::opisthocranion3D.x - XY_origin.x, global::opisthocranion3D.y - XY_origin.y, global::opisthocranion3D.z - XY_origin.z};
    double x_opcranion = XY_opcranion.x * X2D_vector.x + XY_opcranion.y * X2D_vector.y + XY_opcranion.z * X2D_vector.z;
    double y_opcranion = XY_opcranion.x * Y2D_vector.x + XY_opcranion.y * Y2D_vector.y + XY_opcranion.z * Y2D_vector.z;
    global::opisthocranion2D = {x_opcranion, y_opcranion};
    // transform bregma (if present) to specified 2D local coordinate system
    if(landmark::bregma)
    {
        VCoord XY_bregma = {global::bregma3D.x - XY_origin.x, global::bregma3D.y - XY_origin.y, global::bregma3D.z - XY_origin.z};
        double x_bregma = XY_bregma.x * X2D_vector.x + XY_bregma.y * X2D_vector.y + XY_bregma.z * X2D_vector.z;
        double y_bregma = XY_bregma.x * Y2D_vector.x + XY_bregma.y * Y2D_vector.y + XY_bregma.z * Y2D_vector.z;
        global::bregma2D = {x_bregma, y_bregma};
    }
    // return 2D cross-sectional vector
    return crossSection;
}
// function for storing cross-sectional 2D vectors in cvs files
void csvWrite(vector<PCoord> crossSection, string csv_filename)
{
    ofstream outputFile(csv_filename.c_str());
    // check if writing to file is permitted
    if (!outputFile.is_open())
    {
      cout << "Error opening " << csv_filename.c_str() << "for write" << endl;
      terminate();
    }
    else
    {
        for(vector<PCoord>::iterator it = crossSection.begin(); it != crossSection.end(); ++it)
        {
            outputFile << it->x << "," << it->y << "\n";
        }
        outputFile.close();
    }
}
// function of optimizing the location of left and right mastoidale landmarks
vector<VCoord> optimizeMastoidale(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements, string landmark)
{
    vector<VCoord> reducedPoints;
    VCoord normal_T = anatomicalElements[1];
    VCoord mastoidale_new;
    VCoord mastoidale_old;
    if(landmark == "mastoidale left")
    {
        mastoidale_old = anatomicalElements[6];
        if(global::print_out){
            cout << "Optimizing landmark coordinates for " << landmark.c_str() << endl;
            cout << "User coordinates are: x=" << mastoidale_old.x << " y=" << mastoidale_old.y << " z=" << mastoidale_old.z << endl;
        }
    }
    if(landmark == "mastoidale right")
    {
        mastoidale_old = anatomicalElements[7];
        if(global::print_out){
            cout << "Optimizing landmark coordinates for " << landmark.c_str() << endl;
            cout << "User coordinates are: x=" << mastoidale_old.x << " y=" << mastoidale_old.y << " z=" << mastoidale_old.z << endl;
        }
    }
    // keep vertices that are within 8mm from user defined landmark
    for(vector<Mesh>::iterator m_it = MeshElements.begin(); m_it != MeshElements.end(); ++m_it)
    {
        VCoord V1 = {m_it->V1x, m_it->V1y, m_it->V1z};
        VCoord V2 = {m_it->V2x, m_it->V2y, m_it->V2z};
        VCoord V3 = {m_it->V3x, m_it->V3y, m_it->V3z};
        double V1dist = distancePoints3D(mastoidale_old, V1);
        double V2dist = distancePoints3D(mastoidale_old, V2);
        double V3dist = distancePoints3D(mastoidale_old, V3);
        if(V1dist < 8)
        {
            reducedPoints.push_back(V1);
        }
        if(V2dist < 8)
        {
            reducedPoints.push_back(V2);
        }
        if(V3dist < 8)
        {
            reducedPoints.push_back(V3);
        }
    }
    // detect any remaining points lying below original mastoidale landmark and pick the lowest
    // normal_C = {anatomicalElements[2].x, anatomicalElements[2].y, anatomicalElements[2].z};
    double maxdist = 0;
    for(vector<VCoord>::iterator r_it = reducedPoints.begin(); r_it != reducedPoints.end(); ++r_it)
    {
        VCoord point = {r_it->x, r_it->y, r_it->z};
        VCoord pvector = {point.x - mastoidale_old.x, point.y - mastoidale_old.y, point.z - mastoidale_old.z};
        double Pdist = dotProduct(pvector, normal_T);
        if(Pdist <= maxdist)
        {
            maxdist = Pdist;
            mastoidale_new = point;
        }
    }
    // report new landmark coordinates
    if(landmark == "mastoidale left")
    {
        anatomicalElements[6] = mastoidale_new;
        if(global::print_out){
            cout << "Optimal coordinates are: x=" << mastoidale_new.x << " y=" << mastoidale_new.y << " z=" << mastoidale_new.z << endl;
        }
    }
    if(landmark == "mastoidale right")
    {
        anatomicalElements[7] = mastoidale_new;
        if(global::print_out){
            cout << "Optimal coordinates are: x=" << mastoidale_new.x << " y=" << mastoidale_new.y << " z=" << mastoidale_new.z << endl;
        }
    }
    return anatomicalElements;
}
// function for saving cross-sectional segments in csv files
void exportPolylineSegment(vector<PCoord> Section, string csv_filename)
{
    if(global::print_exp)
    {
        string csv = global::obj_filename;
        csv.erase(csv.end() - 4, csv.end());
        csv += csv_filename;
        if(global::print_out){
            cout << "   Writing to " << csv.c_str() << " file..." << flush;
        }
        csvWrite(Section, csv);
        if(global::print_out){
            cout << " OK" << endl;
        }
    }
    
}
// function for rasterizing the 2D points of a planar cross-section into a b&w bitmap. The initial pixel is defined as a
// global variable according to the pixel indices of the point closest to the specified origin of the planarCrossSection
vector<vector<Raster>> rasterizeCrossSection(vector<PCoord> planarCrossSection)
{
    // declare variables
    PCoord origin = {0, 0};
    // find max edge limits in x and y axis
    double left_side = 0;
    double right_side = 0;
    double upper_side = 0;
    double lower_side = 0;
    for(vector<PCoord>::iterator it = planarCrossSection.begin(); it != planarCrossSection.end(); ++it)
    {
        if(it->x < left_side)
        {
            left_side = it->x;
        }
        if(it->x > right_side)
        {
            right_side = it->x;
        }
        if(it->y < lower_side)
        {
            lower_side = it->y;
        }
        if(it->y > upper_side)
        {
            upper_side = it->y;
        }
    }
    // find maximum horizontal and vertical dimensions
    double hor_max_dist = right_side - left_side;
    double ver_max_dist = upper_side - lower_side;
    // get number of necessary pixels for x and y axis by dividing with the given resolution
    int hor_pixels = int(hor_max_dist / global::resolution) + 8;
    int ver_pixels = int(ver_max_dist / global::resolution) + 8;
    if(global::print_out){
        cout << "Raster grid size is " << ver_pixels << " rows by " << hor_pixels << " columns at a resolution of " << global::resolution << "mm/pixel." << endl;
    }
    double left_boundary = left_side - (4 * global::resolution);
    double top_boundary = upper_side + (4 * global::resolution);
    // create a matrix container for pixels and associated pointers to the planar cross section vector
    // the matrix follows rows x columns indexing => rasterMatrix[row_pixel][col_pixel].
    // and each matrix element corresponds to a structure containing a bool pixel and a vector<int> pointer
    vector< vector<Raster> > rasterMatrix;
    // resize matrix according to required pixels, which are set to false by default
    rasterMatrix.resize(ver_pixels, vector<Raster>(hor_pixels, {false, {-1}}));
    // iterate though the planar cross section and allocate each 2D point to rasterMatrix, if a pixel is already set to true
    // then append the 2D point's pointer in the planarCrossSection vector
    for(int it = 0; it < planarCrossSection.size(); ++it)
    {
        // measure distance from left and top boundaries and calculate the respective pixel coordinates
        // starting from the top left corner of the raster matrix
        double dist_2_left = abs(left_boundary) + planarCrossSection[it].x;
        double dist_2_top = top_boundary - planarCrossSection[it].y;
        int row_index = ceil(dist_2_top / global::resolution);
        int col_index = ceil(dist_2_left / global::resolution);
        // check if pixel is already set: if it is (true) then just push_back the pointer of the 2D point being checked
        // if not (false), then set it to true and write the pointer in place of the default 0 value
        if(!rasterMatrix[row_index][col_index].pixel)
        {
            rasterMatrix[row_index][col_index].pixel = true;
            rasterMatrix[row_index][col_index].pointer = {it};
        }
        else
        {
            rasterMatrix[row_index][col_index].pointer.push_back(it);
        }
        // check if point closest to origin and set the global Pixel variable struct
        double distance = distancePoints2D(origin, {planarCrossSection[it].x, planarCrossSection[it].y});
        if(distance < hor_max_dist && distance >= 0)
        {
            hor_max_dist = distance;
            global::InitialPixel = {row_index, col_index, true};
        }
    }
    return rasterMatrix;
}
// function for finding pixel coordinate of specified landmark in the raster image of a planar cross-section
Pixel findLandmarkRasterLocation(vector<PCoord> planarCrossSection, PCoord landmark_coord)
{
    // declare variables
    Pixel landmark_pixel;
    // find max edge limits in x and y axis
    double left_side = 0;
    double right_side = 0;
    double upper_side = 0;
    double lower_side = 0;
    for(vector<PCoord>::iterator it = planarCrossSection.begin(); it != planarCrossSection.end(); ++it)
    {
        if(it->x < left_side)
        {
            left_side = it->x;
        }
        if(it->x > right_side)
        {
            right_side = it->x;
        }
        if(it->y < lower_side)
        {
            lower_side = it->y;
        }
        if(it->y > upper_side)
        {
            upper_side = it->y;
        }
    }
    // calculale horizontal and vertical boundary limits
    double left_boundary = left_side - (4 * global::resolution);
    double top_boundary = upper_side + (4 * global::resolution);
    // find nearest point in the planarCrossSection vector and use it (just in case!)
    PCoord closestPoint;
    double max_distance = abs(right_side) + abs(left_side);
    for(vector<PCoord>::iterator it = planarCrossSection.begin(); it != planarCrossSection.end(); ++it)
    {
        PCoord point = {it->x, it->y};
        double distance = distancePoints2D(landmark_coord, point);
        if(distance < max_distance && distance >= 0)
        {
            max_distance = distance;
            closestPoint = point;
        }
    }
    landmark_coord = closestPoint;
    // measure distance from left and top boundaries and calculate the respective pixel dimensions
    // starting from the top left corner of the raster matrix
    double dist_2_left = abs(left_boundary) + landmark_coord.x;
    double dist_2_top = top_boundary - landmark_coord.y;
    int row_index = ceil(dist_2_top / global::resolution);
    int col_index = ceil(dist_2_left / global::resolution);
    // return landmark 
    landmark_pixel = {row_index, col_index, true};
    return landmark_pixel;
}
// function for calculating the Moore Neighbor algorithm on a b&w raster image given an initial pixel for the contour of interest
// additionally it calculales the respective normal for each pixel in the contour so the direction can be easily traced forward.
// the result is returned in a struct containing pixel coordinates <int>[row_index] <int>[col_index], the associated <int>direction
// and the vector<int> pointer(s) of the original 2D points present in the contour pixel. The Initial pixel location is passed as a
// global variable.
// The implementation follows a counter clockwise rotation so that the outcome vector of the boundary pixels will follow a counter
// clockwise direction around the contour. The initial direction towards the initial pixed is given as a string to avoid tracing
// an inside contour!
vector<Contour> calculateMooreNeighbor(vector<vector<Raster>> rasterMatrix, string init_direction)
{
    // define required variables
    bool contour_open = true;
    vector<Contour> boundaryPixels;
    // populate Neighbor vectors with row,column directions
    // directions in rasterMatrix       [row], [column]
    //              // top -              -1,       0
    //              // top -left          -1,      -1
    //              //     -left           0,      -1
    //              /bottom-left           1,      -1
    //              /bottom-               1,       0
    //              /bottom-right          1,       1
    //              //     -right          0,       1
    //              // top -right         -1,       1
    struct Neighbor
    {
        vector<int> row = {-1, -1, 0, 1, 1, 1, 0, -1};
        vector<int> col = {0, -1, -1, -1, 0, 1, 1, 1};
    };
    Neighbor Mp;
    int direction;
    Pixel p = global::InitialPixel;         //current boundary pixel
    Pixel b;                                //backtrack of current boundary pixel
    // find backtrack pixel of initial entry direction as defined in <string>direction
    // note that initial directions are defined opposite to the Neighbor searching vectors
    // in case the initial pixel is not at the outside border but is masked by another valid
    // pixel preceding from the direction of entry, then the initial pixel is shifted to that
    
    if(init_direction == "upwards")                 //from bottom to top
    {
        int rows = rasterMatrix.size();       //get number of rows
        int row_dist = rows - p.row;          //get number of pixels until raster bottom boundary including initial pixel
        for(int i = 0; i < row_dist; ++i)           //from initial pixel until bottom boundary
        {
            if(rasterMatrix[global::InitialPixel.row + i][p.col].pixel)         //check for active pixels
            {
                p.row = global::InitialPixel.row + i;                           //update initial point's vertical location
                b.col = p.col;
                b.row = p.row + 1;
                direction = 0;
            }
        }
    }
    if(init_direction == "downwards")               //from top to bottom
    {
        for(int i = p.row; i > 0; --i)              //from initial pixel until top boundary
        {
            if(rasterMatrix[i][p.col].pixel)        //check for active pixels
            {
                p.row = i;                          //update initial point's vertical location
                b.col = p.col;
                b.row = p.row - 1;
                direction = 4;
            }
        }
    }
    if(init_direction == "leftwards")               //from right to left
    {
        int cols = rasterMatrix[0].size();    //get number of columns
        int col_dist = cols - p.col;          //get number of pixels until raster right boundary including initial pixel
        for(int i = 0; i < col_dist; ++i)        //from initial pixel until right boundary
        {
            //cout << it << "   ";
            if(rasterMatrix[p.row][global::InitialPixel.col + i].pixel)        //check for active pixels
            {
                
                p.col = global::InitialPixel.col + i;                          //update initial point's horizontal location
                b.row = p.row;
                b.col = p.col + 1;
                direction = 2;
            }
        }
    }
    if(init_direction == "rightwards")              //from left to right
    {
        for(int i = p.col; i > 0; --i)          //from left boundary until initial pixel included
        {
            if(rasterMatrix[p.row][i].pixel)        //check for active pixels
            {
                p.col = i;                          //update initial point's horizontal location
                b.row = p.row;
                b.col = p.col - 1;
                direction = 6;
            }
        }
    }
    // save initial pixel to global variable in case it is altered
    global::InitialPixel = p;
    // insert current boundary pixel to boundary vector
    vector<int> pointers = rasterMatrix[p.row][p.col].pointer;
    Contour valid_pixel = {p.row, p.col, direction, pointers};
    boundaryPixels.push_back(valid_pixel);
    int it = 0;
    // start loop until the contour is closed or maximum number of pixels in rasterMatrix is reached
    int max_it = rasterMatrix.size() * rasterMatrix[0].size();
    while(contour_open && it < max_it)
    {
        Pixel c;                                        //current pixel under consideration
        bool dir_found = false;
        int dir = boundaryPixels[it].direction + 4;     //inverse direction for the next search
        while(!dir_found)
        {
            // start from initial searching direction based on the previous valid pixel
            c.valid = rasterMatrix[p.row + Mp.row[(dir+1)%8]][p.col + Mp.col[(dir+1)%8]].pixel;
            if(c.valid)     //next pixel found!
            {
                c.row = p.row + Mp.row[(dir+1)%8];      //put indices on the validated point under consideration
                c.col = p.col + Mp.col[(dir+1)%8];
                dir_found = true;                       //check found flag to stop rotating aroung the Neighborhood
                b = p;
                p = c;
                // load initial pixels indices with associated pointers into returning vector struct
                int direction = (dir+1)%8;
                vector<int> pointers = rasterMatrix[c.row][c.col].pointer;
                Contour valid_pixel = {c.row, c.col, direction, pointers};
                boundaryPixels.push_back(valid_pixel);
            }
            dir++;
            // check if the contour has been closed. However, this test is initialized after the second loop
            // in order to compare that the closing pixel with the second found along with the direction from
            // initial pixel towards the second one in the contour.
            if(boundaryPixels.size() > 2)
            {
                int row_0 = boundaryPixels[0].row;
                int col_0 = boundaryPixels[0].col;
                
                int row_1 = boundaryPixels[1].row;
                int col_1 = boundaryPixels[1].col;
                int dir_1 = boundaryPixels[1].direction;
                
                int size = boundaryPixels.size() - 1;
                int row_l = boundaryPixels[size].row;
                int col_l = boundaryPixels[size].col;
                int dir_l = boundaryPixels[size].direction;
                
                int row_bl = boundaryPixels[size - 1].row;
                int col_bl = boundaryPixels[size - 1].col;
                int dir_bl = boundaryPixels[size - 1].direction;
                if(row_0 == row_bl && col_0 == col_bl && row_1 == row_l && col_1 == col_l && dir_1 == dir_l)
                {
                    // if contour closed then copy the direction of the initial pixel to its instance at the
                    // beginning of the vector and erase the last element (also second pixel) so that the contour
                    // is complete in the sense of a closed polygon vector.
                    boundaryPixels[0].direction = dir_bl;
                    boundaryPixels.erase(boundaryPixels.end() - 1);
                    // assert the closed contour flag to break the loop
                    contour_open = false;
                }
            }
        }
        it++;
    }
    return boundaryPixels;
}
// function for extracting 2D polyline from the contour vector after boundary tracing the rasterized cross-section
// for each pixel the centroid of the present points is calculaled in order to minimize quantization error due to
// pixelization of the original coordinates
vector<PCoord> extract2Dpolyline(vector<Contour> boundaryPixels, vector<PCoord> planarCrossSection)
{
    // initialize variables
    vector<PCoord> Polyline;
    int poly_length = boundaryPixels.size();
    // start iterating through the boundaryPixels vector
    for(int pixel_it = 0; pixel_it < poly_length; ++pixel_it)
    {
        // check that boundary pixel pointer contains a valid pointer to 2D coordinates (it is set to -1 by default)
        if(boundaryPixels[pixel_it].pointer[0] >= 0)
        {
            // find the corresponding 2D points from the planar cross-section
            int num_pointers = boundaryPixels[pixel_it].pointer.size();
            // initialize 2D centroid for each pixel
            PCoord pixel_centroid = {0, 0};
            for(int pointer_it = 0; pointer_it < num_pointers; ++pointer_it)
            {
                pixel_centroid.x += planarCrossSection[boundaryPixels[pixel_it].pointer[pointer_it]].x;
                pixel_centroid.y += planarCrossSection[boundaryPixels[pixel_it].pointer[pointer_it]].y;
            }
            pixel_centroid.x = pixel_centroid.x / num_pointers;
            pixel_centroid.y = pixel_centroid.y / num_pointers;
            Polyline.push_back(pixel_centroid);
        }
    }
    if(global::print_out){
        cout << "   " << Polyline.size() << " polyline 2D coordinates calculated." << endl;
    }
    return Polyline;
}
// function for extracting specific segments from the contour boundary pixel vector after boundary tracing of the midsagittal curve
vector<Contour> extractMidlineBoundarySegment(vector<Contour> boundaryPixels, Pixel start, Pixel stop)
{
    // initialize variables
    vector<Contour> Segment;
    // make sure that the start and stop pixels are present in the boundaryPixels vector, otherwise find their
    // closest correspondence.
    int start_index;
    int stop_index;
    double distance_i = 5;
    double distance_t = 5;
    for(int it = 0; it < boundaryPixels.size(); ++it)
    {
        double row_i = boundaryPixels[it].col - start.col;
        double col_i = boundaryPixels[it].row - start.row;
        double dist_i = sqrt(row_i * row_i + col_i * col_i);
        if(dist_i < distance_i && dist_i >= 0)
        {
            distance_i = dist_i;
            start_index = it;           //vector index for initial point
        }
        double row_t = boundaryPixels[it].col - stop.col;
        double col_t = boundaryPixels[it].row - stop.row;
        double dist_t = sqrt(row_t * row_t + col_t * col_t);
        if(dist_t < distance_t && dist_t >= 0)
        {
            distance_t = dist_t;
            stop_index = it;            //vector index for terminal point
        }
    }
    // write the closest corresponding pixels back to pixel variables.
    start = {boundaryPixels[start_index].row, boundaryPixels[start_index].col, true};
    stop = {boundaryPixels[stop_index].row, boundaryPixels[stop_index].col, true};
    // start searching for different instances of the pixels
    boundaryPixels.erase(boundaryPixels.end() - 1);     //do not reuse the first point which is also at the end
    int boundary_length = boundaryPixels.size();
    // start iterating through the boundaryPixels vector until the first occurences of the starting and ending pixel coordinates are found
    bool beg_found = false;
    bool end_found = false;
    int pixel_it = 0;
    while((!beg_found || !end_found) && pixel_it < boundary_length - 2)
    {
        // check for starting boundary pixel pointer
        // if the distance to current pixel is less that the distance to the next pixel then select the current as a starting point
        double row_i = boundaryPixels[pixel_it].col - start.col;
        double col_i = boundaryPixels[pixel_it].row - start.row;
        double dist_current = sqrt(row_i * row_i + col_i * col_i);
        row_i = boundaryPixels[pixel_it + 1].col - start.col;
        col_i = boundaryPixels[pixel_it + 1].row - start.row;
        double dist_next = sqrt(row_i * row_i + col_i * col_i);
        if(!beg_found && ((boundaryPixels[pixel_it].row == start.row && boundaryPixels[pixel_it].col == start.col) || (dist_current < 5 && dist_current < dist_next)))
        {
            start_index = pixel_it;
            beg_found = true;
        }
        // check for ending boundary pixel pointer
        // if the distance to current pixel is less that the distance to the next pixel then select the current as a starting point
        row_i = boundaryPixels[pixel_it].col - stop.col;
        col_i = boundaryPixels[pixel_it].row - stop.row;
        dist_current = sqrt(row_i * row_i + col_i * col_i);
        row_i = boundaryPixels[pixel_it + 1].col - stop.col;
        col_i = boundaryPixels[pixel_it + 1].row - stop.row;
        dist_next = sqrt(row_i * row_i + col_i * col_i);
        if(!end_found && ((boundaryPixels[pixel_it].row == stop.row && boundaryPixels[pixel_it].col == stop.col) || (dist_current < 5 && dist_current < dist_next)))
        {
            stop_index = pixel_it;
            end_found = true;
        }
        pixel_it++;
    }
    // copy the appropriate segment into a new Contour vector
    for(int copy_it = start_index; copy_it < stop_index + 1; ++copy_it)
    {
        Segment.push_back(boundaryPixels[copy_it]);
    }
    return Segment;
}
// function for extracting Freeman chain code from contour vectors. Since the chain coding in Moore's Neighbor algorithm
// does not use the standard direction pattern (0 points upwards instead towards the left, the present function makes the
// necessary shift so that the final chain code follows the standard implementation of Freeman chain coding (0 points to
// the positive side of x axis and coding follows a counter clockwise direction.
// Furthermore, the respective displacement along the x and y axes are calculated along with the length of each step.
// Note that boundary pixels' directions regard the direction from previous to current point, whereas Freeman chain code
// regard the direction to the following point. As a result, the Freeman chain code will have N-1 elements (N = number of
// boundary pixels in the segment) corresponding to the k links or steps along the line segment.
vector<FCCode> extractFreemanChainCode(vector<Contour> boundaryPixels)
{
    vector<FCCode> chaincode;
    // start with the 2nd element of the segment
    double sqroot2 = sqrt(2);
    for(int it = 1; it < boundaryPixels.size(); ++it)
    {
        int FCCode = (boundaryPixels[it].direction + 2) % 8;
        double displacement;
        int x, y;
        if(FCCode == 0)
        {
            x = 1; y = 0; displacement = 1;
        }
        else if(FCCode == 1)
        {
            x = 1; y = 1; displacement = sqroot2;
        }
        else if(FCCode == 2)
        {
            x = 0; y = 1; displacement = 1;
        }
        else if(FCCode == 3)
        {
           x = -1; y = 1; displacement = sqroot2;
        }
        else if(FCCode == 4)
        {
            x = -1; y = 0; displacement = 1;
        }
        else if(FCCode == 5)
        {
            x = -1; y = -1; displacement = sqroot2;
        }
        else if(FCCode == 6)
        {
            x = 0; y = -1; displacement = 1;
        }
        else
        {
            x = 1; y = -1; displacement = sqroot2;
        }
        chaincode.push_back({FCCode, x, y, displacement});
    }    
    return chaincode;
}
// function for closing the Freeman Chain Codes of open lines to closed contours by a appending at the end
// an inverted chain of the original chain code
vector<FCCode> closeContourFreemanChainCode(vector<FCCode> FreemanChainCode)
{
    for(int i = FreemanChainCode.size() - 1; i > -1; --i)
    {
        int chain = (FreemanChainCode[i].chain + 4) % 8;
        int delta_x = - FreemanChainCode[i].delta_x;
        int delta_y = - FreemanChainCode[i].delta_y;
        double delta_t = FreemanChainCode[i].delta_t;
        FreemanChainCode.push_back({chain, delta_x, delta_y, delta_t});
    }
    return FreemanChainCode;
}
// function for calculating Elliptic Fourier Descriptors for a closed contour or an open line
// full range Fourier expansion is applied on closed contours and half range expansion to open lines
vector<EFDcoef> calculaleEFD(vector<FCCode> FreemanChainCode)
{
    vector<EFDcoef> EFDs;
    // in case of open line, the closed contour is defined by transversing the outline from final link back to the first
    // by inversing the delta_x and delta_x components.
    double PI = 3.14159265358979323846;
    // calculale the DC components as well as total displacement T along the contour
    double A_0 = 0;
    double C_0 = 0;
    double T = 0;
    double zeta_A = 0;
    double zeta_B = 0;
    double deltaA = 0;
    double deltaB = 0;
    double t_p = 0;
    double t_prev = 0;
    for(int p = 0; p < FreemanChainCode.size(); ++p)
    {
        int delta_x = FreemanChainCode[p].delta_x;
        int delta_y = FreemanChainCode[p].delta_y;
        double delta_t = FreemanChainCode[p].delta_t;
        t_p += delta_t; 
        A_0 += (delta_x / (2 * delta_t)) * (t_p * t_p - t_prev * t_prev) + (zeta_A - (delta_x / delta_t) * zeta_B) * (t_p - t_prev);
        C_0 += (delta_y / (2 * delta_t)) * (t_p * t_p - t_prev * t_prev) + (deltaA - (delta_y / delta_t) * deltaB) * (t_p - t_prev);
        t_prev = t_p;
        if(p > 0)
        {
            zeta_A += FreemanChainCode[p - 1].delta_x;
            zeta_B += FreemanChainCode[p - 1].delta_t;
            deltaA += FreemanChainCode[p - 1].delta_y;
            deltaB += FreemanChainCode[p - 1].delta_t;
        }
        T += FreemanChainCode[p].delta_t;
    }
    A_0 = A_0 / T;
    C_0 = C_0 / T;
    EFDs.push_back({A_0, 0, C_0, 0});
    // calculate the EFD coefficients for a number of harmonics up to half the nummber of links present in the contour
    // for each harmonic calulcate also the respective Fourier power
    vector<EFDcoef> harmonics;
    vector<double> H_power;
    double power = 0;
    for(int n = 1; n < 1 + (FreemanChainCode.size() / 2); ++n)
    {
        // calculate constant terms for each harmonic
        double T_2npi = T / (2 * n * n * PI * PI);
        double pi2n_T = (2 * PI * n) / T;
        double A_n = 0;
        double B_n = 0;
        double C_n = 0;
        double D_n = 0;
        double t_p = 0;         //accumulated length of step segments until current point
        double t_prev = 0;      //accumulated length of step segments until previous point
        for(int p = 0; p < FreemanChainCode.size(); ++p)
        {
            int delta_x = FreemanChainCode[p].delta_x;
            int delta_y = FreemanChainCode[p].delta_y;
            double delta_t = FreemanChainCode[p].delta_t;
            t_p += delta_t;                 //add pth links length to accumulated length
            A_n += (delta_x / delta_t) * (cos(pi2n_T * t_p) - cos(pi2n_T * t_prev));
            B_n += (delta_x / delta_t) * (sin(pi2n_T * t_p) - sin(pi2n_T * t_prev));
            C_n += (delta_y / delta_t) * (cos(pi2n_T * t_p) - cos(pi2n_T * t_prev));
            D_n += (delta_y / delta_t) * (sin(pi2n_T * t_p) - sin(pi2n_T * t_prev));
            t_prev = t_p;                   //save current accumulated length to previous for next iteration
        }
        A_n = A_n * T_2npi;
        B_n = B_n * T_2npi;
        C_n = C_n * T_2npi;
        D_n = D_n * T_2npi;
        harmonics.push_back({A_n, B_n, C_n, D_n});
        power += (A_n*A_n + B_n*B_n + C_n*C_n + D_n*D_n) / 2;
        H_power.push_back(power);
    }
    // for every harmonic check if the average cumulative Fourier power is less than 99.99% and append it in the EFDs vector
    for(int i = 0; i < harmonics.size(); ++i)
    {
        if(H_power[i] < 0.9999 * (H_power[H_power.size() - 1]))
        {
            EFDs.push_back(harmonics[i]);
        }
    }
    if(global::print_out)
    {
        cout << "Number of harmonics up to 99.99% cumulative power spectral density: " << EFDs.size() - 1 << endl;
    }
    return EFDs;
}
// function for reducing a mesh into a point cloud based on the barycenters of each triangular face
vector<VCoord> mesh2PointCloud(vector<Mesh> MeshElements)
{
    // declare variables
    vector<VCoord> pointCloud;
    // for each face in the 3D model, calculate the barycenter and another three points at 2/3 distance from barycenter to each vertex
    // for each face, 4 3D points are returned in the point cloud
    for(vector<Mesh>::iterator it = MeshElements.begin(); it != MeshElements.end(); ++it)
    {
        // calculate barycenter
        double x = (it->V1x + it->V2x + it->V3x) / 3;
        double y = (it->V1y + it->V2y + it->V3y) / 3;
        double z = (it->V1z + it->V2z + it->V3z) / 3;
        VCoord barycenter = {x, y, z};
        // calculate peri-barycentric points
        VCoord peri_BC_V1 = {(barycenter.x + 2 * it->V1x) / 3, (barycenter.y + 2 * it->V1y) / 3,(barycenter.z + 2 * it->V1z) / 3};
        VCoord peri_BC_V2 = {(barycenter.x + 2 * it->V2x) / 3, (barycenter.y + 2 * it->V2y) / 3,(barycenter.z + 2 * it->V2z) / 3};
        VCoord peri_BC_V3 = {(barycenter.x + 2 * it->V3x) / 3, (barycenter.y + 2 * it->V3y) / 3,(barycenter.z + 2 * it->V3z) / 3};
        pointCloud.push_back(barycenter);
        pointCloud.push_back(peri_BC_V1);
        pointCloud.push_back(peri_BC_V2);
        pointCloud.push_back(peri_BC_V3);
    }
    return pointCloud;
}
// function for exporting height map image into csv format
void exportHeightMapImage(vector<vector<HeightMap>> heightMap, string csv_filename)
{
    if(global::print_exp)
    {
        string csv = global::obj_filename;
        csv.erase(csv.end() - 4, csv.end());
        csv += csv_filename;
        if(global::print_out){
            cout << "   Writing to " << csv.c_str() << " file..." << flush;
        }
        ofstream outputFile(csv.c_str());
        // check if writing to file is permitted
        if (!outputFile.is_open())
        {
            cout << "Error opening " << csv.c_str() << "for write" << endl;
            terminate();
        }
        else
        {
            for(int r = 0; r < heightMap.size(); ++r)
            {
                for(int c = 0; c < heightMap[0].size(); ++c)
                {
                    outputFile << heightMap[r][c].height;
                    if(c < heightMap[0].size() - 1)
                    {
                        outputFile << ", ";
                    }
                }
                outputFile << "\n";
            }
            outputFile.close();
        }
        if(global::print_out){
            cout << " OK" << endl;
        }
    }
}
// function for extracting the occipital protuberance area from the cranial mesh based on locations of opisthion and opisthocranion landmarks
vector<Mesh> extractOccipitalMesh(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements)
{
    // declare variables
    vector<Mesh> occipitalMesh;
    VCoord opisthion = global::opisthion3D;
    VCoord opCranion = global::opisthocranion3D;
    // find center of cylinder base as the midpoint of the two landmarks and its normal facings backwards parallel to midsaggital plane
    VCoord Cyl_midPoint = midPoint3D(opisthion, opCranion);
    double Cyl_radius = distancePoints3D(opisthion, opCranion) / 2;
    // cylinder diameter vector facing downwards
    VCoord Cyl_diam_vector = {opisthion.x - opCranion.x, opisthion.y - opCranion.y, opisthion.z - opCranion.z};
    VCoord sagittal_normal = anatomicalElements[0];
    VCoord Cyl_base_normal = crossProduct(sagittal_normal, Cyl_diam_vector);
    Cyl_base_normal = normalizeVector(Cyl_base_normal);
    // for each face in the 3D model, keep only those circumscribed by cylinder facing backwards
    for(vector<Mesh>::iterator it = MeshElements.begin(); it != MeshElements.end(); ++it)
    {
        double dist1 = distancePoints3D(Cyl_midPoint, {it->V1x, it->V1y, it->V1z});
        double dist2 = distancePoints3D(Cyl_midPoint, {it->V2x, it->V2y, it->V2z});
        double dist3 = distancePoints3D(Cyl_midPoint, {it->V3x, it->V3y, it->V3z});
        double dot1 = dotProduct({it->V1x - Cyl_midPoint.x, it->V1y - Cyl_midPoint.y, it->V1z - Cyl_midPoint.z}, Cyl_base_normal);
        double dot2 = dotProduct({it->V2x - Cyl_midPoint.x, it->V2y - Cyl_midPoint.y, it->V2z - Cyl_midPoint.z}, Cyl_base_normal);
        double dot3 = dotProduct({it->V3x - Cyl_midPoint.x, it->V3y - Cyl_midPoint.y, it->V3z - Cyl_midPoint.z}, Cyl_base_normal);
        if(dist1 < Cyl_radius && dist2 < Cyl_radius && dist3 < Cyl_radius && dot1 > 0 && dot2 > 0 && dot3 > 0)
        {
            Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
            occipitalMesh.push_back(tmpmesh);
        }
    }
    if(global::print_out){
        cout << "   " << occipitalMesh.size() << " faces were retained ";
    }
    return occipitalMesh;
}
// function for calculating height map of the occipital protuberance area with 1mm resolution containing the maximum projected
// height of the corresponding vertices with a 0.1mm resolution coded in 8-bit grayscale (0-255).
vector<vector<HeightMap>> calculateOccipitalHeightMap(vector<VCoord> occipitalPointCloud, vector<VCoord> anatomicalElements)
{
    // declare variables
    vector<vector<HeightMap>> occipitalHeightMap;
    VCoord opisthion = global::opisthion3D;
    VCoord opCranion = global::opisthocranion3D;
    // find center of cylinder base as the midpoint of the two landmarks and its normal facings backwards parallel to midsaggital plane
    VCoord Cyl_midPoint = midPoint3D(opisthion, opCranion);
    double Cyl_radius = distancePoints3D(opisthion, opCranion) / 2;
    // cylinder diameter vector facing downwards
    VCoord Cyl_diam_vector = {opisthion.x - opCranion.x, opisthion.y - opCranion.y, opisthion.z - opCranion.z};
    VCoord sagittal_normal = anatomicalElements[0];
    VCoord Cyl_base_normal = crossProduct(sagittal_normal, Cyl_diam_vector);
    Cyl_base_normal = normalizeVector(Cyl_base_normal);
    // re-orient local 2D axes and define 2D origin
    VCoord X2D_vector = {- sagittal_normal.x, - sagittal_normal.y, - sagittal_normal.z};
    VCoord Y2D_vector = {- Cyl_diam_vector.x, - Cyl_diam_vector.y, - Cyl_diam_vector.z};
    Y2D_vector = normalizeVector(Y2D_vector);
    VCoord XY_origin = Cyl_midPoint;
    // calculale top and left boundaries for height map image
    double boundary = Cyl_radius;
    //int hmi_size = ceil(boundary * 2);
    // define fixed raster size for image container
    int img_size = 90;  //////////////////////////ceil(boundary * 2);
    // report circular diameter of projected area
    if(global::print_out){
        cout << "from a " << boundary << "mm diameter circular projection area." << endl;
    }
    // resize matrix according to required pixels, which are set to false by default
    occipitalHeightMap.resize(img_size, vector<HeightMap>(img_size, {255, false}));
    //
    // for every 3D point in the generated point cloud
    for(int i = 0; i < occipitalPointCloud.size(); ++i)
    {
        // transform the 3D points of the point cloud into planar local 2D coordinates and calculale their distance (height)
        // from the projection plane as defined by the base of the cylinder
        VCoord XY = {occipitalPointCloud[i].x - XY_origin.x, occipitalPointCloud[i].y - XY_origin.y, occipitalPointCloud[i].z - XY_origin.z};
        double x = XY.x * X2D_vector.x + XY.y * X2D_vector.y + XY.z * X2D_vector.z;
        double y = XY.x * Y2D_vector.x + XY.y * Y2D_vector.y + XY.z * Y2D_vector.z;
        double height = abs(dotProduct(XY, Cyl_base_normal));
        // convert height to 8 bit grayscale by a resolution of 0.2mm. if height exceeds 255 * 0.1 = 25.5mm, then limit grayscale to 255
        int grayscale = ceil(height / 0.1);
        if(grayscale > 255)
        {
            grayscale = 255;
        }
        // invert grayscale values so that height map gets darker as height increases
        grayscale = 255 - grayscale;
        // calculate raster image location of each projected 3D point with 1mm^2 per pixel area
        double dist_2_left = (img_size / 2) + x;
        double dist_2_top = (img_size / 2) - y;
        int row_index = ceil(dist_2_top);
        int col_index = ceil(dist_2_left);
        // check if row or column indices exceed boundaries of image size and limit them to min-max values
        if(row_index > img_size)
        {
            row_index = img_size;
        }
        if(col_index > img_size)
        {
            col_index = img_size;
        }
        if(row_index < 1)
        {
            row_index = 1;
        }
        if(col_index < 1)
        {
            col_index = 1;
        }
        // check if pixel is already set: if it is (true) then just compare the already set height with the newly calculated one and if
        // the new one is higher, then update the height's value with the larger one. If pixel is not set (false), then set it to true
        // and write the height value in place of the default 0 value.
        if(!occipitalHeightMap[row_index-1][col_index-1].pixel)
        {
            occipitalHeightMap[row_index-1][col_index-1].pixel = true;
            occipitalHeightMap[row_index-1][col_index-1].height = {grayscale};
        }
        else
        {
            if(occipitalHeightMap[row_index-1][col_index-1].height > grayscale)     // to account for the inversion
            {
                occipitalHeightMap[row_index-1][col_index-1].height = {grayscale};
            }
        }
    }
    // scan through the heigh map image for missing pixels due to discretization of 3D points and append an average
    // grayscale value based on the sourrounding active pixels. Missing pixels are considered only those that are still inactive
    // and are surrounded by at least 6 active pixels with non-zero height (i.e. grayscale < 255)
    for(int r = 1; r < img_size - 2; ++r)
    {
        for(int c = 1; c < img_size -2; ++c)
        {
            bool n0 = occipitalHeightMap[r][c].pixel;
            int active = 0;
            int average = 0;
            bool n1 = occipitalHeightMap[r-1][c-1].pixel;
            if(n1){
                active++;
                average += occipitalHeightMap[r-1][c-1].height;
            }
            bool n2 = occipitalHeightMap[r-1][c].pixel;
            if(n2){
                active++;
                average += occipitalHeightMap[r-1][c].height;
            }
            bool n3 = occipitalHeightMap[r-1][c+1].pixel;
            if(n3){
                active++;
                average += occipitalHeightMap[r-1][c+1].height;
            }
            bool n4 = occipitalHeightMap[r][c+1].pixel;
            if(n4){
                active++;
                average += occipitalHeightMap[r][c+1].height;
            }
            bool n5 = occipitalHeightMap[r+1][c+1].pixel;
            if(n5){
                active++;
                average += occipitalHeightMap[r+1][c+1].height;
            }
            bool n6 = occipitalHeightMap[r+1][c].pixel;
            if(n6){
                active++;
                average += occipitalHeightMap[r+1][c].height;
            }
            bool n7 = occipitalHeightMap[r+1][c-1].pixel;
            if(n7){
                active++;
                average += occipitalHeightMap[r+1][c-1].height;
            }
            bool n8 = occipitalHeightMap[r][c-1].pixel;
            if(n8){
                active++;
                average += occipitalHeightMap[r][c-1].height;
            }
            if(!n0 && active >=6)
            {
                int average_grayscale = average / active;
                // append average grayscale value
                occipitalHeightMap[r][c].height = average_grayscale;
                occipitalHeightMap[r][c].pixel = true;
            }
        }
    }
    if(global::print_out){
        cout << "   Height Map Image centered on " << img_size << " by " << img_size << " pixels of 1mm^2 square area each." << endl;
    }
    return occipitalHeightMap;
}
// function for extracting the supraorbital area from the cranial mesh that lies anteriorly of nasion landmark
vector<Mesh> extractSupraOrbitalMesh(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements)
{
    // declare variables
    vector<Mesh> supraorbitalMesh;
    VCoord nasion = global::nasion3D;
    VCoord transverse_normal = anatomicalElements[1];
    VCoord coronal_normal = anatomicalElements[2];
    
    // scan through the entire mesh and find all faces superiorly and anteriorly of the nasion landmark
    for(vector<Mesh>::iterator it = MeshElements.begin(); it != MeshElements.end(); ++it)
    {
        double sup1 = dotProduct({it->V1x - nasion.x, it->V1y - nasion.y, it->V1z - nasion.z}, transverse_normal);
        double sup2 = dotProduct({it->V2x - nasion.x, it->V2y - nasion.y, it->V2z - nasion.z}, transverse_normal);
        double sup3 = dotProduct({it->V3x - nasion.x, it->V3y - nasion.y, it->V3z - nasion.z}, transverse_normal);
        double ant1 = dotProduct({it->V1x - nasion.x, it->V1y - nasion.y, it->V1z - nasion.z}, coronal_normal);
        double ant2 = dotProduct({it->V2x - nasion.x, it->V2y - nasion.y, it->V2z - nasion.z}, coronal_normal);
        double ant3 = dotProduct({it->V3x - nasion.x, it->V3y - nasion.y, it->V3z - nasion.z}, coronal_normal);
        if(sup1 > 0 && sup2 > 0 && sup3 > 0 && ant1 > 0 && ant2 > 0 && ant3 > 0)
        {
            Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
            supraorbitalMesh.push_back(tmpmesh);
        }
    }
    if(global::print_out){
        cout << "   " << supraorbitalMesh.size() << " faces were retained ";
    }
    return supraorbitalMesh;
}
// function for calculating height map of the supraorbital ridge with 1mm resolution containing the maximum projected
// height of the corresponding vertices with a 0.05mm resolution coded in 8-bit grayscale (0-255).
vector<vector<HeightMap>> calculateSupraOrbitalHeightMap(vector<VCoord> supraorbitalPointCloud, vector<VCoord> anatomicalElements)
{
    // declare variables
    vector<vector<HeightMap>> supraorbitalHeightMap;
    VCoord nasion = global::nasion3D;
    VCoord sagittal_normal = anatomicalElements[0];
    VCoord transverse_normal = anatomicalElements[1];
    VCoord coronal_normal = anatomicalElements[2];
    // re-orient local 2D axes and define 2D origin
    VCoord X2D_vector = sagittal_normal;
    VCoord Y2D_vector = transverse_normal;
    VCoord XY_origin = nasion;
    // project 3D points onto nasion coronal plane to find the most superior point to calculale the top boundary for height map image
    // and the most lateral points (left being positive and right being negative on 2D x axis) to find left boundary and width (number
    // of columns) for the height map image
    double top_boundary = 0;
    double left_boundary = 0;
    double right_boundary = 0;
    for(int i = 0; i < supraorbitalPointCloud.size(); ++i)
    {
        // calculale image boundaries
        VCoord XY = {supraorbitalPointCloud[i].x - XY_origin.x, supraorbitalPointCloud[i].y - XY_origin.y, supraorbitalPointCloud[i].z - XY_origin.z};
        double x = XY.x * X2D_vector.x + XY.y * X2D_vector.y + XY.z * X2D_vector.z;
        double y = XY.x * Y2D_vector.x + XY.y * Y2D_vector.y + XY.z * Y2D_vector.z;
        if(y > top_boundary)
        {
            top_boundary = y;
        }
        if(x < left_boundary)
        {
            left_boundary = x;
        }
        if(x > right_boundary)
        {
            right_boundary = x;
        }
    }
    double image_width = right_boundary - left_boundary;
    // calculale height map image dimensions in pixels rows x columns
    int img_width = 60;     ///////ceil(image_width);
    int img_height = 50;    ///////ceil(top_boundary);
    // report rectangular dimensions of projected area
    if(global::print_out){
        cout << "from a " << image_width << "mm x " << top_boundary << "mm rectangular projection area." << endl;
    }
    // resize matrix according to required pixels, which are set to false by default
    supraorbitalHeightMap.resize(img_height, vector<HeightMap>(img_width, {255, false}));
    //
    // for every 3D point in the generated point cloud
    for(int i = 0; i < supraorbitalPointCloud.size(); ++i)
    {
        // transform the 3D points of the point cloud into planar local 2D coordinates and calculale their distance (height)
        // from the projection plane as defined by the nasion coronal plane
        VCoord XY = {supraorbitalPointCloud[i].x - XY_origin.x, supraorbitalPointCloud[i].y - XY_origin.y, supraorbitalPointCloud[i].z - XY_origin.z};
        double x = XY.x * X2D_vector.x + XY.y * X2D_vector.y + XY.z * X2D_vector.z;
        double y = XY.x * Y2D_vector.x + XY.y * Y2D_vector.y + XY.z * Y2D_vector.z;
        double height = abs(dotProduct(XY, coronal_normal));
        // convert height to 8 bit grayscale by a resolution of 0.05mm. if height exceeds 255 * 0.05 = 12.75mm, then limit grayscale to 255
        int grayscale = ceil(height / 0.05);
        if(grayscale > 255)
        {
            grayscale = 255;
        }
        // invert grayscale values so that height map gets darker as height increases
        grayscale = 255 - grayscale;
        // calculate raster image location of each projected 3D point with 1mm^2 per pixel area
        double dist_2_left = (img_width / 2) + x;
        double dist_2_top = img_height - y;
        int row_index = ceil(dist_2_top);
        int col_index = ceil(dist_2_left);
        // check if row or column indices exceed boundaries of image size and limit them to min-max values
        if(row_index > img_height)
        {
            row_index = img_height;
        }
        if(col_index > img_width)
        {
            col_index = img_width;
        }
        if(row_index < 1)
        {
            row_index = 1;
        }
        if(col_index < 1)
        {
            col_index = 1;
        }
        // check if pixel is already set: if it is (true) then just compare the already set height with the newly calculated one and if
        // the new one is higher, then update the height's value with the larger one. If pixel is not set (false), then set it to true
        // and write the height value in place of the default 0 value.
        if(!supraorbitalHeightMap[row_index-1][col_index-1].pixel)
        {
            supraorbitalHeightMap[row_index-1][col_index-1].pixel = true;
            supraorbitalHeightMap[row_index-1][col_index-1].height = {grayscale};
        }
        else
        {
            if(supraorbitalHeightMap[row_index-1][col_index-1].height > grayscale)     // to account for the inversion
            {
                supraorbitalHeightMap[row_index-1][col_index-1].height = {grayscale};
            }
        }
    }
    // scan through the heigh map image for missing pixels due to discretization of 3D points and append an average
    // grayscale value based on the sourrounding active pixels. Missing pixels are considered only those that are still inactive
    // and are surrounded by at least 6 active pixels with non-zero height (i.e. grayscale < 255)
    for(int r = 1; r < img_height - 2; ++r)
    {
        for(int c = 1; c < img_width -2; ++c)
        {
            bool n0 = supraorbitalHeightMap[r][c].pixel;
            int active = 0;
            int average = 0;
            bool n1 = supraorbitalHeightMap[r-1][c-1].pixel;
            if(n1){
                active++;
                average += supraorbitalHeightMap[r-1][c-1].height;
            }
            bool n2 = supraorbitalHeightMap[r-1][c].pixel;
            if(n2){
                active++;
                average += supraorbitalHeightMap[r-1][c].height;
            }
            bool n3 = supraorbitalHeightMap[r-1][c+1].pixel;
            if(n3){
                active++;
                average += supraorbitalHeightMap[r-1][c+1].height;
            }
            bool n4 = supraorbitalHeightMap[r][c+1].pixel;
            if(n4){
                active++;
                average += supraorbitalHeightMap[r][c+1].height;
            }
            bool n5 = supraorbitalHeightMap[r+1][c+1].pixel;
            if(n5){
                active++;
                average += supraorbitalHeightMap[r+1][c+1].height;
            }
            bool n6 = supraorbitalHeightMap[r+1][c].pixel;
            if(n6){
                active++;
                average += supraorbitalHeightMap[r+1][c].height;
            }
            bool n7 = supraorbitalHeightMap[r+1][c-1].pixel;
            if(n7){
                active++;
                average += supraorbitalHeightMap[r+1][c-1].height;
            }
            bool n8 = supraorbitalHeightMap[r][c-1].pixel;
            if(n8){
                active++;
                average += supraorbitalHeightMap[r][c-1].height;
            }
            if(!n0 && active >=6)
            {
                int average_grayscale = average / active;
                // append average grayscale value
                supraorbitalHeightMap[r][c].height = average_grayscale;
                supraorbitalHeightMap[r][c].pixel = true;
            }
        }
    }
    if(global::print_out){
        cout << "   Height Map Image horizontally aligned on " << img_width << " by " << img_height << " pixels of 1mm^2 square area each." << endl;
    }
    return supraorbitalHeightMap;
}
// function for extracting the mastoid process area from the cranial mesh that lies inferioposteriorly of porion landmark
vector<Mesh> extractLateralMastoidProcessMesh(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements, string side)
{
    // declare variables
    vector<Mesh> remainingMesh;
    vector<Mesh> mastoidProcessMesh;
    VCoord mastoidale;
    VCoord porion;
    VCoord medial_margin;
    VCoord sagittal_normal = anatomicalElements[0];
    VCoord transverse_normal = anatomicalElements[1];
    VCoord coronal_normal = anatomicalElements[2];
    // predefine mastoidale displacement at 5mm and porion displacement at 3mm
    double md = 5;
    double pd = 3;
    // get landmark coordinates based on side
    if(side == "left")
    {
        mastoidale = anatomicalElements[6];
        porion = anatomicalElements[9];
        // displace left mastoidale coordinates towards the midsaggital plane (sagittal normal faces to the left),
        // so apply negative displacemt by 5mm
        medial_margin = {mastoidale.x - md * sagittal_normal.x, mastoidale.y - md * sagittal_normal.y, mastoidale.z - md * sagittal_normal.z};
    }
    if(side == "right")
    {
        mastoidale = anatomicalElements[7];
        porion = anatomicalElements[10];
        // displace right mastoidale coordinates towards the midsaggital plane (sagittal normal faces to the left),
        // so apply positive displacemt by 5mm
        medial_margin = {mastoidale.x + md * sagittal_normal.x, mastoidale.y + md * sagittal_normal.y, mastoidale.z + md * sagittal_normal.z};
    }
    // displace porion coordinates anteriorly by 5mm
    porion = {porion.x + pd * coronal_normal.x, porion.y + pd * coronal_normal.y, porion.z + pd * coronal_normal.z};
    // scan through the entire mesh and find all faces inferioposteriorly of the porion landmark
    for(vector<Mesh>::iterator it = MeshElements.begin(); it != MeshElements.end(); ++it)
    {
        double inf1 = dotProduct({it->V1x - porion.x, it->V1y - porion.y, it->V1z - porion.z}, transverse_normal);
        double inf2 = dotProduct({it->V2x - porion.x, it->V2y - porion.y, it->V2z - porion.z}, transverse_normal);
        double inf3 = dotProduct({it->V3x - porion.x, it->V3y - porion.y, it->V3z - porion.z}, transverse_normal);
        double pos1 = dotProduct({it->V1x - porion.x, it->V1y - porion.y, it->V1z - porion.z}, coronal_normal);
        double pos2 = dotProduct({it->V2x - porion.x, it->V2y - porion.y, it->V2z - porion.z}, coronal_normal);
        double pos3 = dotProduct({it->V3x - porion.x, it->V3y - porion.y, it->V3z - porion.z}, coronal_normal);
        if(inf1 < 0 && inf2 < 0 && inf3 < 0 && pos1 < 0 && pos2 < 0 && pos3 < 0)
        {
            Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
            remainingMesh.push_back(tmpmesh);
        }
    }
    // scan through the remaining mesh and and find all faces laterally of the mastoid landmark offset by 5mm towards
    // the midsaggital plane of the cranium
    if(side == "left")
    {
        for(vector<Mesh>::iterator it = remainingMesh.begin(); it != remainingMesh.end(); ++it)
        {
            double lat1 = dotProduct({it->V1x - medial_margin.x, it->V1y - medial_margin.y, it->V1z - medial_margin.z}, sagittal_normal);
            double lat2 = dotProduct({it->V2x - medial_margin.x, it->V2y - medial_margin.y, it->V2z - medial_margin.z}, sagittal_normal);
            double lat3 = dotProduct({it->V3x - medial_margin.x, it->V3y - medial_margin.y, it->V3z - medial_margin.z}, sagittal_normal);
            if(lat1 > 0 && lat2 > 0 && lat3 > 0)
            {
                Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
                mastoidProcessMesh.push_back(tmpmesh);
            }
        }
    }
    if(side == "right")
    {
        for(vector<Mesh>::iterator it = remainingMesh.begin(); it != remainingMesh.end(); ++it)
        {
            double lat1 = dotProduct({it->V1x - medial_margin.x, it->V1y - medial_margin.y, it->V1z - medial_margin.z}, sagittal_normal);
            double lat2 = dotProduct({it->V2x - medial_margin.x, it->V2y - medial_margin.y, it->V2z - medial_margin.z}, sagittal_normal);
            double lat3 = dotProduct({it->V3x - medial_margin.x, it->V3y - medial_margin.y, it->V3z - medial_margin.z}, sagittal_normal);
            if(lat1 < 0 && lat2 < 0 && lat3 < 0)
            {
                Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
                mastoidProcessMesh.push_back(tmpmesh);
            }
        }
    }
    if(global::print_out){
        cout << "   " << mastoidProcessMesh.size() << " faces were retained ";
    }
    return mastoidProcessMesh;
}
// function for calculating lateral height map of the mastoid process with 0.5mm resolution containing the minimum projected
// height of the corresponding vertices with a 0.1mm resolution coded in 8-bit grayscale (0-255). The projection plane
// is set laterally of the mastoidale landmark by 20mm of the mastoid process being analyzed
vector<vector<HeightMap>> calculateLateralMastoidProcessHeightMap(vector<VCoord> mastoidProcessPointCloud, vector<VCoord> anatomicalElements, string side)
{
    // declare variables
    vector<vector<HeightMap>> mastoidProcessHeightMap;
    VCoord mastoidale;
    VCoord lateral_margin;
    VCoord X2D_vector;
    VCoord Y2D_vector;
    VCoord XY_origin;
    VCoord sagittal_normal = anatomicalElements[0];
    VCoord transverse_normal = anatomicalElements[1];
    VCoord coronal_normal = anatomicalElements[2];
    // predefine mastoidale displacement at 20mm
    double md = 20;
    // re-orient local 2D axes and define 2D origin acoording to which mastoid process (left or right) is analyzed
    if(side == "left")
    {
        mastoidale = anatomicalElements[6];
        // positive x axis faces posteriorly
        X2D_vector = {-coronal_normal.x, -coronal_normal.y, -coronal_normal.z};
        // positive y axis faces superiorly
        Y2D_vector = transverse_normal;
        // define XY origin by displacing mastoidale along the sagittal normal by 20mm laterally and 1mm inferiorly along the transverse normal
        lateral_margin = {mastoidale.x + md * sagittal_normal.x, mastoidale.y + md * sagittal_normal.y, mastoidale.z + md * sagittal_normal.z};
        XY_origin = {lateral_margin.x - transverse_normal.x, lateral_margin.y - transverse_normal.y, lateral_margin.z - transverse_normal.z};
    }
    if(side == "right")
    {
        mastoidale = anatomicalElements[7];
        // positive x axis faces posteriorly
        X2D_vector = coronal_normal;
        // positive y axis faces superiorly
        Y2D_vector = transverse_normal;
        // define XY origin by displacing mastoidale along the sagittal normal by 20mm laterally and 1mm inferiorly along the transverse normal
        lateral_margin = {mastoidale.x - md * sagittal_normal.x, mastoidale.y - md * sagittal_normal.y, mastoidale.z - md * sagittal_normal.z};
        XY_origin = {lateral_margin.x - transverse_normal.x, lateral_margin.y - transverse_normal.y, lateral_margin.z - transverse_normal.z};
    }
    //
    // project 3D points onto lateral sagittal plane to find the most superior point to calculale the top boundary for height map image
    // and the most lateral points (anterior being positive and posterior being negative on 2D x axis) to find left boundary and width (number
    // of columns) for the height map image (X2D_vector has already been assigned according to side)
    double top_boundary = 0;
    double left_boundary = 0;
    double right_boundary = 0;
    for(vector<VCoord>::iterator it = mastoidProcessPointCloud.begin(); it != mastoidProcessPointCloud.end(); ++it)
    {
        // calculale image boundaries
        VCoord XY = {it->x - XY_origin.x, it->y - XY_origin.y, it->z - XY_origin.z};
        double x = XY.x * X2D_vector.x + XY.y * X2D_vector.y + XY.z * X2D_vector.z;
        double y = XY.x * Y2D_vector.x + XY.y * Y2D_vector.y + XY.z * Y2D_vector.z;
        if(y > top_boundary)
        {
            top_boundary = y;
        }
        if(x < left_boundary)
        {
            left_boundary = x;
        }
        if(x > right_boundary)
        {
            right_boundary = x;
        }
    }
    double image_width = right_boundary - left_boundary;
    // calculale height map image dimensions in pixels rows x columns
    int img_width = 100;     ///////ceil(image_width);
    int img_height = 80;    ///////ceil(top_boundary);
    // report rectangular dimensions of projected area
    if(global::print_out){
        cout << "from a " << image_width << "mm x " << top_boundary << "mm rectangular projection area." << endl;
    }
    // resize matrix according to required pixels, which are set to false by default
    mastoidProcessHeightMap.resize(img_height, vector<HeightMap>(img_width, {0, false}));
    //
    // for every 3D point in the generated point cloud
    for(vector<VCoord>::iterator it = mastoidProcessPointCloud.begin(); it != mastoidProcessPointCloud.end(); ++it)
    {
        // transform the 3D points of the point cloud into planar local 2D coordinates and calculale their distance (height)
        // from the projection plane as defined by the lateral sagittal plane
        VCoord XY = {it->x - XY_origin.x, it->y - XY_origin.y, it->z - XY_origin.z};
        double x = XY.x * X2D_vector.x + XY.y * X2D_vector.y + XY.z * X2D_vector.z;
        double y = XY.x * Y2D_vector.x + XY.y * Y2D_vector.y + XY.z * Y2D_vector.z;
        double height = abs(dotProduct(XY, sagittal_normal));
        // convert height to 8 bit grayscale by a resolution of 0.1mm. if height exceeds 255 * 0.1 = 25.5mm, then limit grayscale to 255
        int grayscale = ceil(height / 0.1);
        if(grayscale > 255)
        {
            grayscale = 255;
        }
        // invert grayscale values so that height map gets darker as height increases
        grayscale = 255 - grayscale;
        // declare row and column indices
        int row_index;
        int col_index;
        // calculate raster image location of each projected 3D point with 0.5mm^2 per pixel area
        // use left boundary for left mastoid process and right boundary for right mastoid process
        if(side == "left")
        {
            double dist_2_left = abs(left_boundary) + x;
            double dist_2_top = top_boundary - y;
            row_index = ceil(dist_2_top * 2);
            col_index = ceil(dist_2_left * 2);
            // check if row or column indices exceed boundaries of image size and limit them to min-max values
            if(row_index > img_height)
            {
                row_index = img_height;
            }
            if(col_index > img_width)
            {
                col_index = img_width;
            }
            if(row_index < 1)
            {
                row_index = 1;
            }
            if(col_index < 1)
            {
                col_index = 1;
            }
        }
        if(side == "right")
        {
            double dist_2_right = right_boundary - x;
            double dist_2_top = top_boundary - y;
            row_index = ceil(dist_2_top * 2);
            col_index = img_width - floor(dist_2_right * 2);
            // check if row or column indices exceed boundaries of image size and limit them to min-max values
            if(row_index > img_height)
            {
                row_index = img_height;
            }
            if(col_index > img_width)
            {
                col_index = img_width;
            }
            if(row_index < 1)
            {
                row_index = 1;
            }
            if(col_index < 1)
            {
                col_index = 1;
            }
        }
        // check if pixel is already set: if it is (true) then just compare the already set height with the newly calculated one and if
        // the new one is higher, then update the height's value with the larger one. If pixel is not set (false), then set it to true
        // and write the height value in place of the default 0 value.
        if(!mastoidProcessHeightMap[row_index-1][col_index-1].pixel)
        {
            mastoidProcessHeightMap[row_index-1][col_index-1].pixel = true;
            mastoidProcessHeightMap[row_index-1][col_index-1].height = {grayscale};
        }
        else
        {
            if(mastoidProcessHeightMap[row_index-1][col_index-1].height < grayscale)     // to account for the inversion
            {
                mastoidProcessHeightMap[row_index-1][col_index-1].height = {grayscale};
            }
        }
    }
    // scan through the heigh map image for missing pixels due to discretization of 3D points and append an average
    // grayscale value based on the sourrounding active pixels. Missing pixels are considered only those that are still inactive
    // and are surrounded by at least 6 active pixels with non-zero height (i.e. grayscale < 255)
    for(int r = 1; r < img_height - 1; ++r)
    {
        for(int c = 1; c < img_width -1; ++c)
        {
            bool n0 = mastoidProcessHeightMap[r][c].pixel;
            int active = 0;
            int average = 0;
            bool n1 = mastoidProcessHeightMap[r-1][c-1].pixel;
            if(n1){
                active++;
                average += mastoidProcessHeightMap[r-1][c-1].height;
            }
            bool n2 = mastoidProcessHeightMap[r-1][c].pixel;
            if(n2){
                active++;
                average += mastoidProcessHeightMap[r-1][c].height;
            }
            bool n3 = mastoidProcessHeightMap[r-1][c+1].pixel;
            if(n3){
                active++;
                average += mastoidProcessHeightMap[r-1][c+1].height;
            }
            bool n4 = mastoidProcessHeightMap[r][c+1].pixel;
            if(n4){
                active++;
                average += mastoidProcessHeightMap[r][c+1].height;
            }
            bool n5 = mastoidProcessHeightMap[r+1][c+1].pixel;
            if(n5){
                active++;
                average += mastoidProcessHeightMap[r+1][c+1].height;
            }
            bool n6 = mastoidProcessHeightMap[r+1][c].pixel;
            if(n6){
                active++;
                average += mastoidProcessHeightMap[r+1][c].height;
            }
            bool n7 = mastoidProcessHeightMap[r+1][c-1].pixel;
            if(n7){
                active++;
                average += mastoidProcessHeightMap[r+1][c-1].height;
            }
            bool n8 = mastoidProcessHeightMap[r][c-1].pixel;
            if(n8){
                active++;
                average += mastoidProcessHeightMap[r][c-1].height;
            }
            if(!n0 && active >=6)
            {
                int average_grayscale = average / active;
                // append average grayscale value
                mastoidProcessHeightMap[r][c].height = average_grayscale;
                mastoidProcessHeightMap[r][c].pixel = true;
            }
        }
    }
    // scan through the first row of the heigh map image for missing pixels due to discretization of 3D points and append an average
    // grayscale value based on the sourrounding active pixels.
    for(int c = 1; c < img_width -1; ++c)
    {
        bool n0 = mastoidProcessHeightMap[0][c].pixel;
        int active = 0;
        int average = 0;
        bool n1 = mastoidProcessHeightMap[0][c-1].pixel;
        if(n1){
            active++;
            average += mastoidProcessHeightMap[0][c-1].height;
        }
        bool n2 = mastoidProcessHeightMap[0][c+1].pixel;
        if(n2){
            active++;
            average += mastoidProcessHeightMap[0][c+1].height;
        }
        bool n3 = mastoidProcessHeightMap[1][c+1].pixel;
        if(n3){
            active++;
            average += mastoidProcessHeightMap[1][c+1].height;
        }
        bool n4 = mastoidProcessHeightMap[1][c].pixel;
        if(n4){
            active++;
            average += mastoidProcessHeightMap[1][c].height;
        }
        bool n5 = mastoidProcessHeightMap[1][c-1].pixel;
        if(n5){
            active++;
            average += mastoidProcessHeightMap[1][c-1].height;
        }
        if(!n0 && active >3)
        {
            int average_grayscale = average / active;
            // append average grayscale value
            mastoidProcessHeightMap[0][c].height = average_grayscale;
            mastoidProcessHeightMap[0][c].pixel = true;
        }
    }
    // scan the first (left side) or last (right side) column for inactive pixels and register the grayscale value of the laterally
    // adjacent pixel. Just for illustrative purposes of the height map image. Not much of actual usefulness.
    if(side == "left")
    {
        for(int r = 0; r < img_height; ++r)
        {
            bool n0 = mastoidProcessHeightMap[r][0].pixel;
            if(!n0)
            {
                // append grayscale value from adjacent pixel
                mastoidProcessHeightMap[r][0].height = mastoidProcessHeightMap[r][1].height;
                mastoidProcessHeightMap[r][0].pixel = true;
            }
        }
    }
    if(side == "right")
    {
        for(int r = 0; r < img_height; ++r)
        {
            bool n0 = mastoidProcessHeightMap[r][img_width - 1].pixel;
            if(!n0)
            {
                // append grayscale value from adjacent pixel
                mastoidProcessHeightMap[r][img_width - 1].height = mastoidProcessHeightMap[r][img_width - 2].height;
                mastoidProcessHeightMap[r][img_width - 1].pixel = true;
            }
        }
    }
    if(global::print_out){
        cout << "   Height Map Image extracted on " << img_width << " by " << img_height << " pixels of 0.5mm^2 square area each." << endl;
    }
    return mastoidProcessHeightMap;
}
// function for extracting the mastoid process area from the cranial mesh that lies inferiorly of poriom landmark and it is bounded by a rectangular
// box with dimensions 30mm (medial lateral) x 30mm (anterior posterior). The box superior face is coplanar with the transverse plane passing at porion
// landmark and its center is vertically located above mastoidale landmark
vector<Mesh> extractInferiorMastoidProcessMesh(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements, string side)
{
    // declare variables
    vector<Mesh> remainingMesh;
    vector<Mesh> mastoidProcessMesh;
    VCoord mastoidale;
    VCoord porion;
    VCoord medial_margin;
    VCoord lateral_margin;
    VCoord anterior_margin;
    VCoord posterior_margin;
    VCoord sagittal_normal = anatomicalElements[0];
    VCoord transverse_normal = anatomicalElements[1];
    VCoord coronal_normal = anatomicalElements[2];
    // predefine rectangular base dimensions at 30mm x 30mm and vertical distance from porion at 40mm
    double mel = 15; //from center to sides
    // get landmark coordinates based on side
    if(side == "left")
    {
        mastoidale = anatomicalElements[6];
        porion = anatomicalElements[9];
        // displace left mastoidale coordinates towards the midsaggital plane (sagittal normal faces to the left)
        medial_margin = {mastoidale.x - mel * sagittal_normal.x, mastoidale.y - mel * sagittal_normal.y, mastoidale.z - mel * sagittal_normal.z};
        // displace left mastoidale coordinates mediolaterally (sagittal normal faces to the left)
        lateral_margin = {mastoidale.x + mel * sagittal_normal.x, mastoidale.y + mel * sagittal_normal.y, mastoidale.z + mel * sagittal_normal.z};        
    }
    if(side == "right")
    {
        mastoidale = anatomicalElements[7];
        porion = anatomicalElements[10];
        // displace right mastoidale coordinates towards the midsaggital plane (sagittal normal faces to the left)
        medial_margin = {mastoidale.x + mel * sagittal_normal.x, mastoidale.y + mel * sagittal_normal.y, mastoidale.z + mel * sagittal_normal.z};
        // displace right mastoidale coordinates mediolaterally (sagittal normal faces to the left)
        lateral_margin = {mastoidale.x - mel * sagittal_normal.x, mastoidale.y - mel * sagittal_normal.y, mastoidale.z - mel * sagittal_normal.z};        
    }
    // displace mastoidale coordinates anteriorly (coronal normal faces forward)
    anterior_margin = {mastoidale.x + mel * coronal_normal.x, mastoidale.y + mel * coronal_normal.y, mastoidale.z + mel * coronal_normal.z};
    // displace mastoidale coordinates posteriorly (coronal normal faces forward)
    posterior_margin = {mastoidale.x - mel * coronal_normal.x, mastoidale.y - mel * coronal_normal.y, mastoidale.z - mel * coronal_normal.z};
    //
    // scan through the entire mesh and find all faces inferiorly of the porion landmark and within the box's anteroposterior edges
    for(vector<Mesh>::iterator it = MeshElements.begin(); it != MeshElements.end(); ++it)
    {
        double inf1 = dotProduct({it->V1x - porion.x, it->V1y - porion.y, it->V1z - porion.z}, transverse_normal);
        double inf2 = dotProduct({it->V2x - porion.x, it->V2y - porion.y, it->V2z - porion.z}, transverse_normal);
        double inf3 = dotProduct({it->V3x - porion.x, it->V3y - porion.y, it->V3z - porion.z}, transverse_normal);
        
        double ant1 = dotProduct({it->V1x - anterior_margin.x, it->V1y - anterior_margin.y, it->V1z - anterior_margin.z}, coronal_normal);
        double ant2 = dotProduct({it->V2x - anterior_margin.x, it->V2y - anterior_margin.y, it->V2z - anterior_margin.z}, coronal_normal);
        double ant3 = dotProduct({it->V3x - anterior_margin.x, it->V3y - anterior_margin.y, it->V3z - anterior_margin.z}, coronal_normal);
        
        double pos1 = dotProduct({it->V1x - posterior_margin.x, it->V1y - posterior_margin.y, it->V1z - posterior_margin.z}, coronal_normal);
        double pos2 = dotProduct({it->V2x - posterior_margin.x, it->V2y - posterior_margin.y, it->V2z - posterior_margin.z}, coronal_normal);
        double pos3 = dotProduct({it->V3x - posterior_margin.x, it->V3y - posterior_margin.y, it->V3z - posterior_margin.z}, coronal_normal);
        if(inf1 < 0 && inf2 < 0 && inf3 < 0 && ant1 < 0 && ant2 < 0 && ant3 < 0 && pos1 > 0 && pos2 > 0 && pos3 > 0)
        {
            Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
            remainingMesh.push_back(tmpmesh);
        }
    }
    // scan through the remaining mesh and and find all faces within the box's mediolateral edges
    if(side == "left")
    {
        for(vector<Mesh>::iterator it = remainingMesh.begin(); it != remainingMesh.end(); ++it)
        {
            double lat1 = dotProduct({it->V1x - lateral_margin.x, it->V1y - lateral_margin.y, it->V1z - lateral_margin.z}, sagittal_normal);
            double lat2 = dotProduct({it->V2x - lateral_margin.x, it->V2y - lateral_margin.y, it->V2z - lateral_margin.z}, sagittal_normal);
            double lat3 = dotProduct({it->V3x - lateral_margin.x, it->V3y - lateral_margin.y, it->V3z - lateral_margin.z}, sagittal_normal);
            
            double med1 = dotProduct({it->V1x - medial_margin.x, it->V1y - medial_margin.y, it->V1z - medial_margin.z}, sagittal_normal);
            double med2 = dotProduct({it->V2x - medial_margin.x, it->V2y - medial_margin.y, it->V2z - medial_margin.z}, sagittal_normal);
            double med3 = dotProduct({it->V3x - medial_margin.x, it->V3y - medial_margin.y, it->V3z - medial_margin.z}, sagittal_normal);
            if(lat1 < 0 && lat2 < 0 && lat3 < 0 && med1 > 0 && med2 > 0 && med3 > 0)
            {
                Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
                mastoidProcessMesh.push_back(tmpmesh);
            }
        }
    }
    if(side == "right")
    {
        for(vector<Mesh>::iterator it = remainingMesh.begin(); it != remainingMesh.end(); ++it)
        {
            double lat1 = dotProduct({it->V1x - lateral_margin.x, it->V1y - lateral_margin.y, it->V1z - lateral_margin.z}, sagittal_normal);
            double lat2 = dotProduct({it->V2x - lateral_margin.x, it->V2y - lateral_margin.y, it->V2z - lateral_margin.z}, sagittal_normal);
            double lat3 = dotProduct({it->V3x - lateral_margin.x, it->V3y - lateral_margin.y, it->V3z - lateral_margin.z}, sagittal_normal);
            
            double med1 = dotProduct({it->V1x - medial_margin.x, it->V1y - medial_margin.y, it->V1z - medial_margin.z}, sagittal_normal);
            double med2 = dotProduct({it->V2x - medial_margin.x, it->V2y - medial_margin.y, it->V2z - medial_margin.z}, sagittal_normal);
            double med3 = dotProduct({it->V3x - medial_margin.x, it->V3y - medial_margin.y, it->V3z - medial_margin.z}, sagittal_normal);
            if(lat1 > 0 && lat2 > 0 && lat3 > 0 && med1 < 0 && med2 < 0 && med3 < 0)
            {
                Mesh tmpmesh = {it->V1x, it->V1y, it->V1z, it->V2x, it->V2y, it->V2z, it->V3x, it->V3y, it->V3z};
                mastoidProcessMesh.push_back(tmpmesh);
            }
        }
    }
    if(global::print_out){
        cout << "   " << mastoidProcessMesh.size() << " faces were retained ";
    }
    return mastoidProcessMesh;
}
// function for calculating inferior height map of the mastoid process with 0.5mm resolution containing the minimum projected
// height of the corresponding vertices with a 0.15mm resolution coded in 8-bit grayscale (0-255). The projection plane
// is set inferiorly of the respective porion landmark by 40mm of the mastoid process being analyzed
vector<vector<HeightMap>> calculateInferiorMastoidProcessHeightMap(vector<VCoord> mastoidProcessPointCloud, vector<VCoord> anatomicalElements, string side)
{
    // declare variables
    vector<vector<HeightMap>> mastoidProcessHeightMap;
    VCoord mastoidale;
    VCoord porion;
    VCoord superior_margin;
    VCoord X2D_vector;
    VCoord Y2D_vector;
    VCoord XY_origin;
    VCoord sagittal_normal = anatomicalElements[0];
    VCoord transverse_normal = anatomicalElements[1];
    VCoord coronal_normal = anatomicalElements[2];
    // predefine porion displacement at 40mm
    double pd = 40;
    // re-orient local 2D axes and define 2D origin acoording to which mastoid process (left or right) is analyzed
    if(side == "left")
    {
        mastoidale = anatomicalElements[6];
        porion = anatomicalElements[9];
        // positive x axis faces to the left hand side of the cranium
        X2D_vector = sagittal_normal;
        // positive y axis faces anteriorly
        Y2D_vector = coronal_normal;
        // define XY origin by mastoidale along the transverse normal by such amount so that its distance from the transverse plane passing at
        // porion is 40mm. First find the mastoidale projection to the transverse plane defined by porion
        VCoord v = {mastoidale.x - porion.x, mastoidale.y - porion.y, mastoidale.z - porion.z};
        double dp = dotProduct(v, transverse_normal);
        superior_margin = {mastoidale.x - dp * transverse_normal.x, mastoidale.y - dp * transverse_normal.y, mastoidale.z - dp * transverse_normal.z};
        XY_origin = {superior_margin.x - pd * transverse_normal.x, superior_margin.y - pd * transverse_normal.y, superior_margin.z - pd * transverse_normal.z};
    }
    if(side == "right")
    {
        mastoidale = anatomicalElements[7];
        porion = anatomicalElements[10];
        // positive x axis faces posteriorly
        X2D_vector = sagittal_normal;
        // positive y axis faces anteriorly
        Y2D_vector = coronal_normal;
        // define XY origin by mastoidale along the transverse normal by such amount so that its distance from the transverse plane passing at
        // porion is 40mm.
        VCoord v = {mastoidale.x - porion.x, mastoidale.y - porion.y, mastoidale.z - porion.z};
        double dp = dotProduct(v, transverse_normal);
        superior_margin = {mastoidale.x - dp * transverse_normal.x, mastoidale.y - dp * transverse_normal.y, mastoidale.z - dp * transverse_normal.z};
        XY_origin = {superior_margin.x - pd * transverse_normal.x, superior_margin.y - pd * transverse_normal.y, superior_margin.z - pd * transverse_normal.z};
    }
    //
    // calculale height map image dimensions in pixels rows x columns
    int img_size = 60;      //30mm @ 0.5mm resolution = 60 pixels
    // report rectangular dimensions of projected area
    if(global::print_out){
        cout << "from a 30mm x 30mm rectangular projection area." << endl;
    }
    // resize matrix according to required pixels, which are set to false by default
    mastoidProcessHeightMap.resize(img_size, vector<HeightMap>(img_size, {0, false}));
    //
    // for every 3D point in the generated point cloud
    for(vector<VCoord>::iterator it = mastoidProcessPointCloud.begin(); it != mastoidProcessPointCloud.end(); ++it)
    {
        // transform the 3D points of the point cloud into planar local 2D coordinates and calculale their distance (height)
        // from the projection plane as defined by the lateral sagittal plane
        VCoord XY = {it->x - XY_origin.x, it->y - XY_origin.y, it->z - XY_origin.z};
        double x = XY.x * X2D_vector.x + XY.y * X2D_vector.y + XY.z * X2D_vector.z;
        double y = XY.x * Y2D_vector.x + XY.y * Y2D_vector.y + XY.z * Y2D_vector.z;
        double height = abs(dotProduct(XY, transverse_normal));
        // convert height to 8 bit grayscale by a resolution of 0.15mm. if height exceeds 255 * 0.15 = 38.25mm, then limit grayscale to 255
        int grayscale = ceil(height / 0.15);
        if(grayscale > 255)
        {
            grayscale = 255;
        }
        // invert grayscale values so that height map gets darker as height increases
        grayscale = 255 - grayscale;
        //
        // calculate raster image location of each projected 3D point with 0.5mm^2 per pixel area
        double dist_2_left = (img_size / 4) + x;
        double dist_2_top = (img_size / 4) - y;
        int row_index = ceil(dist_2_top * 2);
        int col_index = ceil(dist_2_left * 2);
        // check if row or column indices exceed boundaries of image size and limit them to min-max values
        if(row_index > img_size)
        {
            row_index = img_size;
        }
        if(col_index > img_size)
        {
            col_index = img_size;
        }
        if(row_index < 1)
        {
            row_index = 1;
        }
        if(col_index < 1)
        {
            col_index = 1;
        }
        // check if pixel is already set: if it is (true) then just compare the already set height with the newly calculated one and if
        // the new one is higher, then update the height's value with the larger one. If pixel is not set (false), then set it to true
        // and write the height value in place of the default 0 value.
        if(!mastoidProcessHeightMap[row_index-1][col_index-1].pixel)
        {
            mastoidProcessHeightMap[row_index-1][col_index-1].pixel = true;
            mastoidProcessHeightMap[row_index-1][col_index-1].height = {grayscale};
        }
        else
        {
            if(mastoidProcessHeightMap[row_index-1][col_index-1].height < grayscale)     // to account for the inversion
            {
                mastoidProcessHeightMap[row_index-1][col_index-1].height = {grayscale};
            }
        }
    }
    // scan through the heigh map image for missing pixels due to discretization of 3D points and append an average
    // grayscale value based on the sourrounding active pixels. Missing pixels are considered only those that are still inactive
    // and are surrounded by at least 6 active pixels with non-zero height (i.e. grayscale < 255)
    for(int r = 1; r < img_size - 1; ++r)
    {
        for(int c = 1; c < img_size -1; ++c)
        {
            bool n0 = mastoidProcessHeightMap[r][c].pixel;
            int active = 0;
            int average = 0;
            bool n1 = mastoidProcessHeightMap[r-1][c-1].pixel;
            if(n1){
                active++;
                average += mastoidProcessHeightMap[r-1][c-1].height;
            }
            bool n2 = mastoidProcessHeightMap[r-1][c].pixel;
            if(n2){
                active++;
                average += mastoidProcessHeightMap[r-1][c].height;
            }
            bool n3 = mastoidProcessHeightMap[r-1][c+1].pixel;
            if(n3){
                active++;
                average += mastoidProcessHeightMap[r-1][c+1].height;
            }
            bool n4 = mastoidProcessHeightMap[r][c+1].pixel;
            if(n4){
                active++;
                average += mastoidProcessHeightMap[r][c+1].height;
            }
            bool n5 = mastoidProcessHeightMap[r+1][c+1].pixel;
            if(n5){
                active++;
                average += mastoidProcessHeightMap[r+1][c+1].height;
            }
            bool n6 = mastoidProcessHeightMap[r+1][c].pixel;
            if(n6){
                active++;
                average += mastoidProcessHeightMap[r+1][c].height;
            }
            bool n7 = mastoidProcessHeightMap[r+1][c-1].pixel;
            if(n7){
                active++;
                average += mastoidProcessHeightMap[r+1][c-1].height;
            }
            bool n8 = mastoidProcessHeightMap[r][c-1].pixel;
            if(n8){
                active++;
                average += mastoidProcessHeightMap[r][c-1].height;
            }
            if(!n0 && active >=6)
            {
                int average_grayscale = average / active;
                // append average grayscale value
                mastoidProcessHeightMap[r][c].height = average_grayscale;
                mastoidProcessHeightMap[r][c].pixel = true;
            }
        }
    }
    // scan through the first and last rows of the heigh map image for missing pixels due to discretization of 3D points
    // and append an average grayscale value based on the sourrounding active pixels.
    for(int c = 1; c < img_size -1; ++c)
    {
        bool n0 = mastoidProcessHeightMap[0][c].pixel;
        int active = 0;
        int average = 0;
        bool n1 = mastoidProcessHeightMap[0][c-1].pixel;
        if(n1){
            active++;
            average += mastoidProcessHeightMap[0][c-1].height;
        }
        bool n2 = mastoidProcessHeightMap[0][c+1].pixel;
        if(n2){
            active++;
            average += mastoidProcessHeightMap[0][c+1].height;
        }
        bool n3 = mastoidProcessHeightMap[1][c+1].pixel;
        if(n3){
            active++;
            average += mastoidProcessHeightMap[1][c+1].height;
        }
        bool n4 = mastoidProcessHeightMap[1][c].pixel;
        if(n4){
            active++;
            average += mastoidProcessHeightMap[1][c].height;
        }
        bool n5 = mastoidProcessHeightMap[1][c-1].pixel;
        if(n5){
            active++;
            average += mastoidProcessHeightMap[1][c-1].height;
        }
        if(!n0 && active >=3)
        {
            int average_grayscale = average / active;
            // append average grayscale value
            mastoidProcessHeightMap[0][c].height = average_grayscale;
            mastoidProcessHeightMap[0][c].pixel = true;
        }
    }
    for(int c = 1; c < img_size -1; ++c)
    {
        bool n0 = mastoidProcessHeightMap[img_size - 1][c].pixel;
        int active = 0;
        int average = 0;
        bool n1 = mastoidProcessHeightMap[img_size - 1][c-1].pixel;
        if(n1){
            active++;
            average += mastoidProcessHeightMap[img_size - 1][c-1].height;
        }
        bool n2 = mastoidProcessHeightMap[img_size - 1][c+1].pixel;
        if(n2){
            active++;
            average += mastoidProcessHeightMap[img_size - 1][c+1].height;
        }
        bool n3 = mastoidProcessHeightMap[img_size - 2][c+1].pixel;
        if(n3){
            active++;
            average += mastoidProcessHeightMap[img_size - 2][c+1].height;
        }
        bool n4 = mastoidProcessHeightMap[img_size - 2][c].pixel;
        if(n4){
            active++;
            average += mastoidProcessHeightMap[img_size - 2][c].height;
        }
        bool n5 = mastoidProcessHeightMap[img_size - 2][c-1].pixel;
        if(n5){
            active++;
            average += mastoidProcessHeightMap[img_size - 2][c-1].height;
        }
        if(!n0 && active >=3)
        {
            int average_grayscale = average / active;
            // append average grayscale value
            mastoidProcessHeightMap[img_size - 1][c].height = average_grayscale;
            mastoidProcessHeightMap[img_size - 1][c].pixel = true;
        }
    }
    // scan the first (left side) or last (right side) column for inactive pixels and register the grayscale value of the laterally
    // adjacent pixel. Just for illustrative purposes of the height map image. Not much of actual usefulness.
    if(side == "left")
    {
        for(int r = 0; r < img_size; ++r)
        {
            bool n0 = mastoidProcessHeightMap[r][0].pixel;
            if(!n0)
            {
                // append grayscale value from adjacent pixel
                mastoidProcessHeightMap[r][0].height = mastoidProcessHeightMap[r][1].height;
                mastoidProcessHeightMap[r][0].pixel = true;
            }
        }
    }
    if(side == "right")
    {
        for(int r = 0; r < img_size; ++r)
        {
            bool n0 = mastoidProcessHeightMap[r][img_size - 1].pixel;
            if(!n0)
            {
                // append grayscale value from adjacent pixel
                mastoidProcessHeightMap[r][img_size - 1].height = mastoidProcessHeightMap[r][img_size - 2].height;
                mastoidProcessHeightMap[r][img_size - 1].pixel = true;
            }
        }
    }
    if(global::print_out){
        cout << "   Height Map Image extracted on " << img_size << " by " << img_size << " pixels of 0.5mm^2 square area each." << endl;
    }
    return mastoidProcessHeightMap;
}
// function for saving results data into an appropriate GNU Octave data text format.
// The data includes the Elliptic Fourier Descriptors along with their respective Freeman chaincodes
// and all the measurements stored in the measurements::variables namespace
void saveResultsData2Octave()
{
    if(global::print_oct)
    {
        string extension = ".mat";
        string mat = global::obj_filename;
        mat.erase(mat.end() - 4, mat.end());
        mat += extension;
        ofstream outputFile(mat.c_str());
        // check if writing to file is permitted
        if (!outputFile.is_open())
        {
            cout << "Error opening " << mat.c_str() << "for write" << endl;
            terminate();
        }
        else
        {
            // get current local time and date
            time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
            // writing header to file
            outputFile << "# Created by skullanalyzer 1.0.0, ";
            outputFile << "at " << put_time(localtime(&now), "%T") << " on " << put_time(localtime(&now), "%F") << "\n";
            //////
            // FreemanChainCode structure
            outputFile << "# name: FCC\n# type: scalar struct\n# ndims: 2\n 1 1\n";
            outputFile << "# length: 1\n";
            outputFile << "# name: NasionBregma\n# type: matrix\n# rows: 1\n# columns: " << FreemanChainCode::NasionBregma.size() << "\n";
            for(vector<FCCode>::iterator it = FreemanChainCode::NasionBregma.begin(); it != FreemanChainCode::NasionBregma.end(); ++it)
            {
                outputFile << " " << it->chain;
            }
            outputFile << "\n\n\n";
            outputFile << "\n\n";
            //////
            // Eliptic Fourier Descriptors structure
            outputFile << "# name: EFD\n# type: scalar struct\n# ndims: 2\n 1 1\n";
            outputFile << "# length: 1\n";
            outputFile << "# name: NasionBregma\n# type: matrix\n# rows: " << EFD::NasionBregma.size() << "\n# columns: 4\n";
            for(vector<EFDcoef>::iterator it = EFD::NasionBregma.begin(); it != EFD::NasionBregma.end(); ++it)
            {
                outputFile << " " << it->a << " " << it->b << " " << it->c << " " << it->d << "\n";
            }
            outputFile << "\n\n";
            outputFile << "\n\n";
            //////
            // Height Map Image structure
            outputFile << "# name: HMI\n# type: scalar struct\n# ndims: 2\n 1 1\n";
            if(landmark::mastoidaleL && landmark::mastoidaleR && landmark::nasion_optimal)
            {
                outputFile << "# length: 6\n";
            }
            else if(landmark::mastoidaleL && landmark::mastoidaleR && !landmark::nasion_optimal)
            {
                outputFile << "# length: 5\n";
            }
            else if((landmark::mastoidaleL || landmark::mastoidaleR) && landmark::nasion_optimal)
            {
                outputFile << "# length: 4\n";
            }
            else if((landmark::mastoidaleL || landmark::mastoidaleR) && !landmark::nasion_optimal)
            {
                outputFile << "# length: 3\n";
            }
            else if(landmark::nasion_optimal)
            {
                outputFile << "# length: 2\n";
            }
            else
            {
                outputFile << "# length: 1\n";
            }
            outputFile << "# name: OccipitalProtuberance\n# type: matrix\n# rows: " << HMI::opHeightMapImage.size() << "\n";
            outputFile << "# columns: " << HMI::opHeightMapImage[0].size() << "\n";
            for(int r = 0; r < HMI::opHeightMapImage.size(); ++r)
            {
                for(int c = 0; c < HMI::opHeightMapImage[0].size(); ++c)
                {
                    outputFile << HMI::opHeightMapImage[r][c].height;
                    if(c < HMI::opHeightMapImage[0].size() - 1)
                    {
                        outputFile << " ";
                    }
                }
                outputFile << "\n";
            }
            outputFile << "\n\n";
            if(landmark::nasion_optimal)
            {
                outputFile << "# name: SupraorbitalRidge\n# type: matrix\n# rows: " << HMI::srHeightMapImage.size() << "\n";
                outputFile << "# columns: " << HMI::srHeightMapImage[0].size() << "\n";
                for(int r = 0; r < HMI::srHeightMapImage.size(); ++r)
                {
                    for(int c = 0; c < HMI::srHeightMapImage[0].size(); ++c)
                    {
                        outputFile << HMI::srHeightMapImage[r][c].height;
                        if(c < HMI::srHeightMapImage[0].size() - 1)
                        {
                            outputFile << " ";
                        }
                    }
                    outputFile << "\n";
                }
            }
            outputFile << "\n\n";
            if(landmark::mastoidaleL)
            {
                outputFile << "# name: LeftMastoidLateral\n# type: matrix\n# rows: " << HMI::lmLatHeightMapImage.size() << "\n";
                outputFile << "# columns: " << HMI::lmLatHeightMapImage[0].size() << "\n";
                for(int r = 0; r < HMI::lmLatHeightMapImage.size(); ++r)
                {
                    for(int c = 0; c < HMI::lmLatHeightMapImage[0].size(); ++c)
                    {
                        outputFile << HMI::lmLatHeightMapImage[r][c].height;
                        if(c < HMI::lmLatHeightMapImage[0].size() - 1)
                        {
                            outputFile << " ";
                        }
                    }
                    outputFile << "\n";
                }
                outputFile << "\n\n";
                outputFile << "# name: LeftMastoidInferior\n# type: matrix\n# rows: " << HMI::lmInfHeightMapImage.size() << "\n";
                outputFile << "# columns: " << HMI::lmInfHeightMapImage[0].size() << "\n";
                for(int r = 0; r < HMI::lmInfHeightMapImage.size(); ++r)
                {
                    for(int c = 0; c < HMI::lmInfHeightMapImage[0].size(); ++c)
                    {
                        outputFile << HMI::lmInfHeightMapImage[r][c].height;
                        if(c < HMI::lmInfHeightMapImage[0].size() - 1)
                        {
                            outputFile << " ";
                        }
                    }
                    outputFile << "\n";
                }
                outputFile << "\n\n";
            }
            if(landmark::mastoidaleR)
            {
                outputFile << "# name: RightMastoidLateral\n# type: matrix\n# rows: " << HMI::rmLatHeightMapImage.size() << "\n";
                outputFile << "# columns: " << HMI::rmLatHeightMapImage[0].size() << "\n";
                for(int r = 0; r < HMI::rmLatHeightMapImage.size(); ++r)
                {
                    for(int c = 0; c < HMI::rmLatHeightMapImage[0].size(); ++c)
                    {
                        outputFile << HMI::rmLatHeightMapImage[r][c].height;
                        if(c < HMI::rmLatHeightMapImage[0].size() - 1)
                        {
                            outputFile << " ";
                        }
                    }
                    outputFile << "\n";
                }
                outputFile << "\n\n";
                outputFile << "# name: RightMastoidInferior\n# type: matrix\n# rows: " << HMI::rmInfHeightMapImage.size() << "\n";
                outputFile << "# columns: " << HMI::rmInfHeightMapImage[0].size() << "\n";
                for(int r = 0; r < HMI::rmInfHeightMapImage.size(); ++r)
                {
                    for(int c = 0; c < HMI::rmInfHeightMapImage[0].size(); ++c)
                    {
                        outputFile << HMI::rmInfHeightMapImage[r][c].height;
                        if(c < HMI::rmInfHeightMapImage[0].size() - 1)
                        {
                            outputFile << " ";
                        }
                    }
                    outputFile << "\n";
                }
                outputFile << "\n\n";
            }
            outputFile << "\n\n\n";
            // inform user
            outputFile.close();
            if(global::print_out)
            {
                cout << "Data saved to " << mat.c_str() << endl;
            }
        }
    }
}
//
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
//
// group of functions for specific geometric features
//
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
//
// function for calculating parameters of the nasion midsagittal section
void nasionBregmaSegment(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements, vector<LMcoord> LMpoints)
{
    // initialize 3D point vectors
    vector<VCoord> intersectionPlane;
    vector<VCoord> MidSagPoints;
    // initialize 2D cross-section vectors
    vector<PCoord> MidSagSection;
    vector<PCoord> MidSagPolyline;
    
    // mute console output during this part of the function (if applicable)
    bool toggle = false;
    if(global::print_out)
    {
        global::print_out = false;
        toggle = true;
    }
    // take cross-section at nasion along the sagittal plane
    intersectionPlane = chooseSlicingPlane(anatomicalElements, "nasion", "sagittal");
    MidSagPoints = sliceMesh(MeshElements, intersectionPlane);
    // define 2D local x axis colinear with the coronal normal
    VCoord X2D_vector = {anatomicalElements[2].x, anatomicalElements[2].y, anatomicalElements[2].z};
    // project to 2D coordinates
    MidSagSection = project3Dto2Dplane(MidSagPoints, intersectionPlane, X2D_vector);
    // rasterize mid sagittal planar section
    vector<vector<Raster>> rasterMatrix = rasterizeCrossSection(MidSagSection);
    // apply Moore Neighbor tracing algorithm
    vector<Contour> boundaryPixels = calculateMooreNeighbor(rasterMatrix, "leftwards");
    // unmute console output for the rest of the function (if applicable)
    if(toggle)
    {
        global::print_out = true;
    }
    //
    // optimize nasion landmark location
    optimizeNasion(boundaryPixels, MidSagSection, intersectionPlane, X2D_vector);
    // update anatomicalElements vector with optimized nasion 3D coordinates
    anatomicalElements[3] = global::nasion3D;
    //
    // go through the whole process again to get the proper cross-sectional coordinates
    //
    // take cross-section at nasion along the sagittal plane
    intersectionPlane = chooseSlicingPlane(anatomicalElements, "nasion", "sagittal");
    MidSagPoints = sliceMesh(MeshElements, intersectionPlane);
    // find optimized nasion-opisthocranion max distance and opisthocranion 3D coordinates
    horMaxDistanceNasionOpisthocranion(MidSagPoints, anatomicalElements);
    //
    // write the optimized nasion and opisthocranion landmark coordinates to a new Meshlab .pp file
    if(global::print_lmk)
    {
        string pp = global::pp_filename;
        pp.erase(pp.end() - 3, pp.end());
        pp += "_a.pp";
        double LMx = global::nasion3D.x;
        double LMy = global::nasion3D.y;
        double LMz = global::nasion3D.z;
        LMpoints[0] = {1, LMx, LMy, LMz};
        LMx = global::opisthocranion3D.x;
        LMy = global::opisthocranion3D.y;
        LMz = global::opisthocranion3D.z;
        LMpoints.push_back({9, LMx, LMy, LMz});
        writeMeshlabPoints(pp , LMpoints);
    }
    // extract the Nasion - Bregma segment (if applicable) of the boundary pixels as separate 2D polyline
    if(landmark::bregma)
    {
        if(global::print_out){
            cout << "Extracting Nasion - Bregma midsagittal curve:" << endl;
        }
        // take new projection to 2D plane
        MidSagSection = project3Dto2Dplane(MidSagPoints, intersectionPlane, X2D_vector);
        // rasterize mid sagittal planar section and save it to csv
        rasterMatrix = rasterizeCrossSection(MidSagSection);
        // apply Moore Neighbor tracing algorithm
        boundaryPixels = calculateMooreNeighbor(rasterMatrix, "leftwards");
        // find bregma location on raster image
        Pixel bregma = findLandmarkRasterLocation(MidSagSection, global::bregma2D);
        // extract nasion-bregma segment and convert it to polyline
        vector<Contour> NasionBregma = extractMidlineBoundarySegment(boundaryPixels, global::InitialPixel, bregma);
        vector<PCoord> NasionBregmaPoly = extract2Dpolyline(NasionBregma, MidSagSection);
        exportPolylineSegment(NasionBregmaPoly, "_NasionBregmaSegment.csv");
        // calculale Freeman Chain Code and Elliptic Fourier Descriptors
        FreemanChainCode::NasionBregma = extractFreemanChainCode(NasionBregma);
        FreemanChainCode::NasionBregma = closeContourFreemanChainCode(FreemanChainCode::NasionBregma);
        EFD::NasionBregma = calculaleEFD(FreemanChainCode::NasionBregma);
    }
    else
    {
        if(global::print_out){
            cout << "Nasion - Bregma segment omitted due landmark absense." << endl;
        }
        FreemanChainCode::NasionBregma.push_back({0, 0, 0, 0});
        EFD::NasionBregma.push_back({0, 0, 0, 0});
    }
}
// function for calculating height map of the occipital protuberance area
void occipitalHeightMap(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements)
{
    if(global::print_out){
        cout << "Extracting occipital bone height map image:" << endl;
    }
    vector<Mesh> occipitalMesh = extractOccipitalMesh(MeshElements, anatomicalElements);
    vector<VCoord> occipitalPointCloud = mesh2PointCloud(occipitalMesh);
    HMI::opHeightMapImage = calculateOccipitalHeightMap(occipitalPointCloud, anatomicalElements);
    exportHeightMapImage(HMI::opHeightMapImage, "_occipitalHeightMapImage.csv");
}
// function for calculating height map of the supraorbital ridge
void supraorbitalHeightMap(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements)
{
    if(global::print_out){
        cout << "Extracting supraorbital ridge height map image:" << endl;
    }
    vector<Mesh> supraorbitalMesh = extractSupraOrbitalMesh(MeshElements, anatomicalElements);
    vector<VCoord> supraorbitalPointCloud = mesh2PointCloud(supraorbitalMesh);
    HMI::srHeightMapImage = calculateSupraOrbitalHeightMap(supraorbitalPointCloud, anatomicalElements);
    exportHeightMapImage(HMI::srHeightMapImage, "_supraorbitalHeightMapImage.csv");
}
// function for calculating parameters of the left mastoidale landmark
void leftMastoidProcessHeightMaps(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements, vector<LMcoord> LMpoints)
{    
    // write the optimized landmark coordinates to a new Meshlab .pp file
    if(global::print_lmk)
    {
        string pp = global::pp_filename;
        pp.erase(pp.end() - 3, pp.end());
        pp += "_a.pp";
        double LMx = global::nasion3D.x;
        double LMy = global::nasion3D.y;
        double LMz = global::nasion3D.z;
        LMpoints[0] = {1, LMx, LMy, LMz};
        LMx = global::opisthocranion3D.x;
        LMy = global::opisthocranion3D.y;
        LMz = global::opisthocranion3D.z;
        LMpoints.push_back({9, LMx, LMy, LMz});
        LMx = anatomicalElements[6].x;
        LMy = anatomicalElements[6].y;
        LMz = anatomicalElements[6].z;
        LMpoints[3] = {4, LMx, LMy, LMz};
        writeMeshlabPoints(pp , LMpoints);
    }
    
    // lateral height map image generation
    if(global::print_out){
        cout << "Extracting left mastoid process lateral height map image:" << endl;
    }
    vector<Mesh> mastoidProcessMesh = extractLateralMastoidProcessMesh(MeshElements, anatomicalElements, "left");
    vector<VCoord> mastoidProcessPointCloud = mesh2PointCloud(mastoidProcessMesh);
    HMI::lmLatHeightMapImage = calculateLateralMastoidProcessHeightMap(mastoidProcessPointCloud, anatomicalElements, "left");
    exportHeightMapImage(HMI::lmLatHeightMapImage, "_mastoidLeftLateralHeightMapImage.csv");
    // inferior height map image generation
    if(global::print_out){
        cout << "Extracting left mastoid process inferior height map image:" << endl;
    }
    mastoidProcessMesh = extractInferiorMastoidProcessMesh(MeshElements, anatomicalElements, "left");
    mastoidProcessPointCloud = mesh2PointCloud(mastoidProcessMesh);
    HMI::lmInfHeightMapImage = calculateInferiorMastoidProcessHeightMap(mastoidProcessPointCloud, anatomicalElements, "left");
    exportHeightMapImage(HMI::lmInfHeightMapImage, "_mastoidLeftInferiorHeightMapImage.csv");
}
// function for calculating parameters of the right mastoidale landmark
void rightMastoidProcessHeightMaps(vector<Mesh> MeshElements, vector<VCoord> anatomicalElements, vector<LMcoord> LMpoints)
{    
    // write the optimized landmark coordinates to a new Meshlab .pp file
    if(global::print_lmk)
    {
        string pp = global::pp_filename;
        pp.erase(pp.end() - 3, pp.end());
        pp += "_a.pp";
        double LMx = global::nasion3D.x;
        double LMy = global::nasion3D.y;
        double LMz = global::nasion3D.z;
        LMpoints[0] = {1, LMx, LMy, LMz};
        LMx = global::opisthocranion3D.x;
        LMy = global::opisthocranion3D.y;
        LMz = global::opisthocranion3D.z;
        LMpoints.push_back({9, LMx, LMy, LMz});
        if(landmark::mastoidaleL)
        {
            LMx = anatomicalElements[6].x;
            LMy = anatomicalElements[6].y;
            LMz = anatomicalElements[6].z;
            LMpoints[3] = {4, LMx, LMy, LMz};
        }
        LMx = anatomicalElements[7].x;
        LMy = anatomicalElements[7].y;
        LMz = anatomicalElements[7].z;
        LMpoints[4] = {5, LMx, LMy, LMz};
        writeMeshlabPoints(pp , LMpoints);
    }
    
    // lateral height map image generation
    if(global::print_out){
        cout << "Extracting right mastoid process lateral height map image:" << endl;
    }
    vector<Mesh> mastoidProcessMesh = extractLateralMastoidProcessMesh(MeshElements, anatomicalElements, "right");
    vector<VCoord> mastoidProcessPointCloud = mesh2PointCloud(mastoidProcessMesh);
    HMI::rmLatHeightMapImage = calculateLateralMastoidProcessHeightMap(mastoidProcessPointCloud, anatomicalElements, "right");
    exportHeightMapImage(HMI::rmLatHeightMapImage, "_mastoidRightLateralHeightMapImage.csv");
    // inferior height map image generation
    if(global::print_out){
        cout << "Extracting right mastoid process inferior height map image:" << endl;
    }
    mastoidProcessMesh = extractInferiorMastoidProcessMesh(MeshElements, anatomicalElements, "right");
    mastoidProcessPointCloud = mesh2PointCloud(mastoidProcessMesh);
    HMI::rmInfHeightMapImage = calculateInferiorMastoidProcessHeightMap(mastoidProcessPointCloud, anatomicalElements, "right");
    exportHeightMapImage(HMI::rmInfHeightMapImage, "_mastoidRightInferiorHeightMapImage.csv");
}
//
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
//
// main function
int main(int argc, char **argv)
{
    // check user input arguments: exclude the program's name
    checkUserArguments(argc - 1, argv + 1);
    //
    // read vertices and faces from the mesh file
    vector<Mesh> MeshElements = readObjMesh(global::obj_filename);
    // read landmark coordinates from pp file
    vector<LMcoord> LMpoints = readMeshlabPoints(global::pp_filename);
    // calculate normals of anatomical planes and order available cranial landmarks
    vector<VCoord> anatomicalElements = anatomicalPlanes(LMpoints);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////
    //// calculate all cranial geometric features according to which landmarks are present
    ////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //
    //////
    nasionBregmaSegment(MeshElements, anatomicalElements, LMpoints);
    //
    // update anatomicalElements vector with optimized nasion 3D coordinates
    anatomicalElements[3] = global::nasion3D;
    
    //////
    occipitalHeightMap(MeshElements, anatomicalElements);
    
    // check if nasion_optimal landmark is true before calculating the supraorbital Height Map Image
    if(landmark::nasion_optimal)
    {
        //////
        supraorbitalHeightMap(MeshElements, anatomicalElements);
    }
    else
    {
        if(global::print_out){
            cout << "No protruding supraorbital ridge to extract." << endl;
        }
    }
    
    // check if mastoidale landmarks are present and calculate the respective cross-sections
    // based on their optimized location at the very lowest end of the mastoid process
    if(landmark::mastoidaleL)
    {
        // find optimal position of the mastoidale landmark positioned at the very bottom of the mastoid process
        anatomicalElements = optimizeMastoidale(MeshElements, anatomicalElements, "mastoidale left");
        leftMastoidProcessHeightMaps(MeshElements, anatomicalElements, LMpoints);
    }
    if(landmark::mastoidaleR)
    {
        // find optimal position of the mastoidale landmark positioned at the very bottom of the mastoid process
        anatomicalElements = optimizeMastoidale(MeshElements, anatomicalElements, "mastoidale right");
        rightMastoidProcessHeightMaps(MeshElements, anatomicalElements, LMpoints);
    }
    saveResultsData2Octave();
    return 0;
}
