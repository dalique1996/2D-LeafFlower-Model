//
// This file is part of MorphoDynamX - http://www.MorphoDynamX.org
// Copyright (C) 2012-2017 Richard S. Smith and collaborators.
//
// If you use MorphoDynamX in your work, please cite:
//   http://dx.doi.org/10.7554/eLife.05864
//
// MorphoDynamX is free software, and is licensed under under the terms of the
// GNU General (GPL) Public License version 2.0, http://www.gnu.org/licenses.
//
#ifndef FLOWER_HPP
#define FLOWER_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <Attributes.hpp>
#include <Function.hpp>
#include <MDXProcessTissue.hpp>
#include <MeshProcessSystem.hpp>
#include <Contour.hpp>
#include <CCTopology.hpp>
#include <Mesh.hpp>
#include <MeshBuilder.hpp>
#include <Process.hpp>
#include <MDXProcessCellDivide.hpp>
#include <MeshProcessStructure.hpp>
#include <ToolsProcessSystem.hpp>
#include <CCVerify.hpp>
#include <CCIndex.hpp>
#include <CellMakerMesh2D.hpp>
#include <GeomMathUtils.hpp>
#include <Geometry.hpp>
#include <Matrix.hpp>
#include <Triangulate.hpp>
#include <limits>
#include <math.h>
#include <CCDivideCell.hpp>
#include <cstdlib>
#include "tbb/concurrent_unordered_map.h"
#include <fstream>
#include <chrono>
#include <complex>

//GENERAL PARAMETERS
using namespace mdx;
double INF = std::numeric_limits<double>::max();
const double EPS = 1e-6; //const constant value, cannot be modified;
const double BIG_VAL = 100000;
int LABEL=1;
Point3d refPosition;
Point3d ref_pos;
CCIndex refVertex;
double FinalStep=72000;

//GENERAL STRUCTURES
struct edgeStruct
      {
    double intercellularAuxin;
    double posPin1;
    double negPin1;
    double posPrevExport;
    double negPrevExport;

    Point3d fBound; //position
    Point3d sBound;
       };


struct edgeParData
      {
      int label;
      double pin;
      double edge_export_prev;
       };


//MIX FUNCTIONS
//for polarity calculation
void MBR(std::vector<Point3d> points, std::vector<Point3d>& rect) {
    std::vector<int> permutation;
    sortPolygonPoints(points, Point3d(0, 0, -1), permutation);

    std::vector<Point3d> p;
    for(auto i : permutation)
        p.push_back(points[i]);
    p.push_back(p[0]);

    std::vector<Point3d> e;
    for(uint i = 1; i < p.size(); i++)
        e.push_back(p[i] - p[i - 1]);

    std::vector<double> norms;
    for(auto i : e)
        norms.push_back(norm(i));

    std::vector<Point3d> v;
    for(uint i = 0; i < e.size(); i++)
        v.push_back(e[i] / norms[i]);

    std::vector<Point3d> w;
    for(uint i = 0; i < v.size(); i++)
        w.push_back(Point3d(-v[i][1], v[i][0], 0));

    std::vector<double> a, b, areas;
    std::vector<std::vector<double>> x, y;
    for(uint i = 0; i < v.size(); i++) {
        a.clear();
        b.clear();
        for(uint j = 0; j < p.size(); j++) {
            a.push_back(p[j] * v[i]);
            b.push_back(p[j] * w[i]);
        }

        double max_a = *max_element(a.begin(), a.end());
        double min_a = *min_element(a.begin(), a.end());
        double max_b = *max_element(b.begin(), b.end());
        double min_b = *min_element(b.begin(), b.end());
        areas.push_back((min_b - max_b) * (min_a - max_a));

        x.push_back({min_a, max_a});
        y.push_back({min_b, max_b});
    }

    int k = std::min_element(areas.begin(), areas.end()) - areas.begin();
    int qx[4] = {0, 1, 1, 0};
    int qy[4] = {0, 0, 1, 1};
    Matrix<4, 2, double> M1;
    for(int i = 0; i < 4; i++)
        M1[i] = Point2d(x[k][qx[i]], y[k][qy[i]]);
    Matrix<2, 3, double> M2;
    M2[0] = v[k];
    M2[1] = w[k];
    Matrix<4, 3, double> M = M1 * M2;
    rect.clear();
    for(int i = 0; i < 4; i++)
        rect.push_back(M[i]);
}

Point3d findClosestLineToLine(Point3d targetLine,
                           Point3d line1, Point3d line2) {
    Point3d finalLine;
    double minAxis1;
    if(mdx::angle(targetLine, line1) < mdx::angle(targetLine, -line1))
         minAxis1 = mdx::angle(targetLine, line1);
    else
         minAxis1 = mdx::angle(targetLine, -line1);
    double minAxis2;
    if(mdx::angle(targetLine, line2) < mdx::angle(targetLine, -line2))
         minAxis2 = mdx::angle(targetLine, line2);
    else
         minAxis2 = mdx::angle(targetLine, -line2);
    if(minAxis1 < minAxis2)
        finalLine = line1;
    else
        finalLine = line2;
    return finalLine;
}

Point3d shortenLength(Point3d A, float reductionLength) {
    Point3d B = A;
    B *= (1 - reductionLength / A.norm());
    return B;
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%T", &tstruct);

    return buf;
}

Point3d Centroid(std::vector<Point3d> vertices, int size)
{
    Point3d centroid {0.0, 0.0, 0.0};
    double signedArea = 0.0;
    double x0 = 0.0; // Current vertex X
    double y0 = 0.0; // Current vertex Y
    double x1 = 0.0; // Next vertex X
    double y1 = 0.0; // Next vertex Y
    double a = 0.0;  // Partial signed area

    // For all vertices except last
    int i=0;
    for (i=0; i<size-1; ++i)
    {
        x0 = vertices[i][0];
        y0 = vertices[i][1];
        x1 = vertices[i+1][0];
        y1 = vertices[i+1][1];
        a = x0*y1 - x1*y0;
        signedArea += a;
        centroid[0] += (x0 + x1)*a;
        centroid[1] += (y0 + y1)*a;
    }

    // Do last vertex separately to avoid performing an expensive
    // modulus operation in each iteration.
    x0 = vertices[i][0];
    y0 = vertices[i][1];
    x1 = vertices[0][0];
    y1 = vertices[0][1];
    a = x0*y1 - x1*y0;
    signedArea += a;
    centroid[0] += (x0 + x1)*a;
    centroid[1] += (y0 + y1)*a;

    signedArea *= 0.5;
    centroid[0] /= (6.0*signedArea);
    centroid[1] /= (6.0*signedArea);

    return centroid;
}

Point3d RotateClock(Point3d v, double angle)
{
    Point3d p = v;
    double Tx=p[0],Ty=p[1];

    double X=(Tx*cos(angle))+(Ty*sin(angle));
    double Y=(Ty*cos(angle))-(Tx*sin(angle));

    p[0] = X;
    p[1] = Y;

    return p;
}

Point3d RotateCounterClock(Point3d v, double angle)
{
    Point3d p = v;
    double Tx=p[0],Ty=p[1];

    double X=(Tx*cos(angle))-(Ty*sin(angle));
    double Y=(Ty*cos(angle))+(Tx*sin(angle));

    p[0] = X;
    p[1] = Y;

    return p;
}

Point3d Rotate90(Point3d v) //y>0?¿
{
    Point3d p = v;
    if (v[0]>=0){
        p[0]=-v[1];
        p[1]=v[0];
    }
    else if (v[0]<0){
        p[0]=v[1];
        p[1]=-v[0];
    }

    return p;
}

//check if the intersect Point is in the segment defined by p2q2
bool onSegment(Point3d p2, Point3d q2, Point3d intersect)
{
  double maxX, minX, maxY, minY;
  if(p2[0]>=q2[0]){
      maxX=p2[0];
      minX=q2[0];
  } else{
      maxX=q2[0];
      minX=p2[0];
   }

  if(p2[1]>=q2[1]){
      maxY=p2[1];
      minY=q2[1];
  } else{
      maxY=q2[1];
      minY=p2[1];
   }

  if(intersect[0]<=maxX && intersect[0]>=minX && intersect[1]<=maxY && intersect[1]>=minY)
      return true;
  else
      return false;
}

//intersect line p1q1, segement p2q2
bool lineSegmentIntersection(Point3d p1, Point3d q1, Point3d p2, Point3d q2, Point3d& intersect) //FIRST lineLineIntersection
{
    // Line p1q1 represented as a1x + b1y = c1
    double a1 = q1[1] - p1[1];
    double b1 = p1[0] - q1[0];
    double c1 = a1*(p1[0]) + b1*(p1[1]);
    // Line p2q2 represented as a2x + b2y = c2
    double a2 = q2[1] - p2[1];
    double b2 = p2[0] - q2[0];
    double c2 = a2*(p2[0])+ b2*(p2[1]);

    double determinant = a1*b2 - a2*b1;

    if (determinant == 0) //lines parallel
        return false;

    else //get intersection point
    {
        intersect[0] = (b2*c1 - b1*c2)/determinant;
        intersect[1] = (a1*c2 - a2*c1)/determinant;
        if(onSegment(p2, q2, intersect))//SECOND check if the intersection point is in the edge (it's a segment)
            return true;
        else
            return false;
    }
}

//sort vector and keep the track of indices
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}


template < typename T>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, const T  & element)
{
        std::pair<bool, int > result;

        // Find given element in vector
        auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

        if (it != vecOfElements.end())
        {
                result.second = distance(vecOfElements.begin(), it);
                result.first = true;
        }
        else
        {
                result.first = false;
                result.second = -1;
        }

        return result;
}


namespace Flower
{
struct VerticeData;
struct EdgeData;
struct FaceData;
typedef AttrMap<CCIndex, VerticeData> VerticeDataAttr;
typedef AttrMap<CCIndex, EdgeData> EdgeDataAttr;
typedef AttrMap<CCIndex, FaceData> FaceDataAttr;

struct Data {
    int dim = -1; // initialize dimension
};


////----------------FACE DATA----------------------//
enum CellSide {
    UndefSide,
    left,
    right
};
static CellSide stringToCellSide(const QString& str) {
    if(str == "UndefSide")
        return (UndefSide);
    else if(str == "Left")
        return (left);
    else if(str == "Right")
        return (right);
    else
        throw(QString("Bad cell side %1").arg(str));
}

enum CellLayer {
    UndefLayer,
    L1,
    L2,
    L3
};
static CellLayer stringToCellLayer(const QString& str) {
    if(str == "UndefLayer")
        return (UndefLayer);
    else if(str == "L1")
        return (L1);
    else if(str == "L2")
        return (L2);
    else if(str == "L3")
        return (L3);
    else
        throw(QString("Bad cell type %1").arg(str));
}

enum CellZone {
    UndefZone,
    CZ,
    PZ,
    AD,
    AB,
    MD
};
static CellZone stringToCellZone(const QString& str) {
    if(str == "UndefZone")
        return (UndefZone);
    else if(str == "CZ")
        return (CZ);
    else if(str == "PZ")
        return (PZ);
    else if(str == "AD")
        return (AD);
    else if(str == "AB")
        return (AB);
    else if(str == "MD")
        return (MD);
    else
        throw(QString("Bad cell type %1").arg(str));
}

enum CellPosition {
    Nothing,
    max,
    nine,
    eight,
    seven,
    six,
    five,
    four,
    three,
    two,
    one,
    pointfive,
    pointone
};

static CellPosition stringToCellPosition(const QString& str) {
    if(str == "Nothing")
        return (Nothing);
    else if(str == "1")
        return (max);
    else if(str == "0.9")
        return (nine);
    else if(str == "0.8")
        return (eight);
    else if(str == "0.7")
        return (seven);
    else if(str == "0.6")
        return (six);
    else if(str == "0.5")
        return (five);
    else if(str == "0.4")
        return (four);
    else if(str == "0.3")
        return (three);
    else if(str == "0.2")
        return (two);
    else if(str == "0.1")
        return (one);
    else if(str == "0.05")
        return (pointfive);
    else if(str == "0.01")
        return (pointone);
    else
        throw(QString("Bad cell position %1").arg(str));
}

enum InitialAuxin {
    hundredX,
    eightyX,
    fiftyX,
    twentyfiveX,
    tenX,
    fiveX,
    twoX,
    oneX
};

static InitialAuxin stringToInitialAuxin(const QString& str) {
    if(str == "100")
        return (hundredX);
    else if(str == "80")
        return (eightyX);
    else if(str == "50")
        return (fiftyX);
    else if(str == "25")
        return (twentyfiveX);
    else if(str == "10")
        return (tenX);
    else if(str == "5")
        return (fiveX);
    else if(str == "2")
        return (twoX);
    else if(str == "1")
        return (oneX);
    else
        throw(QString("Bad initial auxin type %1").arg(str));
}

struct FaceData: Data
{

  //VISUAL
    uint cellSide = UndefSide;
    uint cellLayer = UndefLayer;
    uint cellZone = UndefZone;
    uint cellPosition = Nothing;
    uint initialAuxin = oneX;
    bool OC = false;
    Point3d newcentroid, oldcentroid;
    Point3d force;
    Point3d stress;
    double  stiffnessXX;
    Point3d initStrain;
    Point3d strain;
    Point3d a1, a2;
    Point3d a1_inc;
    double prevArea;

  //Connect meshes
  std::vector<CCIndex> cellFaces;

  //Triangulation
  std::vector<CCIndex> crossVertex;
  std::vector<CCIndexVec> triangles;

  //Chemical
  double auxin=0; // Auxin
  double auxinProdRate=0;
  double auxinDecayRate=0;
  double a = 0; // activator
  int b = 0; //inhibitor
  double pin1=0; //Pin1 cytoplasmic
  Point3d auxinFluxVector=Point3d(0,0,0);
  Point3d auxinGradientVector=Point3d(0,0,0);

  //Mechanical
  double pressure;
  double strainStiffnessXX;
  std::vector<Point3d> restX0; //ShapePBD
  std::vector<CCIndex> cellVertices; //std::vector<CCIndex>* cellVertices =0; as a pointer
  Point3d restCm;
  Matrix3d invRestMat; //for shape
  bool shapeInit;
  std::vector<double> angles;
  Matrix2d strainInvRestMat; //for strain

  FaceData() {dim = 2;}

  bool read(const QByteArray& ba, size_t& pos) {
      readPOD<uint>(cellSide, ba, pos);
      readPOD<uint>(cellLayer, ba, pos);
      readPOD<uint>(cellZone, ba, pos);
      readPOD<bool>(OC, ba, pos);
      readPOD<uint>(cellPosition, ba, pos);
      readPOD<uint>(initialAuxin, ba, pos);
      readPOD<double>(a, ba, pos);
      readPOD<int>(b, ba, pos);
      readPOD<Point3d>(restCm, ba, pos);
      readPOD<bool>(shapeInit, ba, pos);
      return true;
  }

  bool write(QByteArray& ba) {
      writePOD<uint>(cellSide, ba);
      writePOD<uint>(cellLayer, ba);
      writePOD<uint>(cellZone, ba);
      writePOD<bool>(OC, ba);
      writePOD<uint>(cellPosition, ba);
      writePOD<uint>(initialAuxin, ba);
      writePOD<double>(a, ba);
      writePOD<int>(b, ba);
      writePOD<Point3d>(restCm, ba);
      writePOD<bool>(shapeInit, ba);
      return true;
  }

  bool operator==(const FaceData &other) const
  {
    return false;
  }
};

////----------------EDGE DATA----------------------//
struct EdgeData: Data
{
    int a=0;
    bool border =false;
    CCIndex f, g; // Attached faces, for fast lookup

    //Connect meshes
    bool wall=false;
    CCIndex sameEdge;

    //Mechanical
    double restLength=0;
    double stiffness=0;

    //Chemical
    double intercellularAuxin=0;
    Point3d midPoint=Point3d(0,0,0);
    std::map<int, double> pin1,traffickedPin1;
    std::map<int, double> edge_export, edge_export_prev;
    std::map<int, double> auxinFluxImpact, MFImpact;
    std::map<int, double> pin1Sensitivity, pin1SensitivityRaw;

    EdgeData() {dim=1;}

    bool read(const QByteArray& ba, size_t& pos) {
       readPOD<int>(a, ba, pos);
       return true;
    }

    bool write(QByteArray& ba) {
       writePOD<int>(a, ba);
       return true;
    }

    bool operator==(const EdgeData &other) const
    {
       return false;
    }
};

////----------------VERTEX DATA----------------------//
struct VerticeData : Data
{
    //
    int a=0;
    bool bottom=false;
    bool stiffRefVertex=false;

    //Mechanical
    Point3d velocity = Point3d(0,0,0);//PBD
    Point3d prevVelocity = Point3d(0,0,0);
    double invMass =1.0; //mass equal to 1;
    Point3d prevPos, lastPos;
    Point3d force;
    std::map<const char*, Point3d> corrections;

    VerticeData() {dim = 0;}

    bool read(const QByteArray& ba, size_t& pos) {

        readPOD<int>(a, ba, pos);
        readPOD<bool>(bottom, ba, pos);
        readPOD<double>(invMass, ba, pos);
        readPOD<bool>(stiffRefVertex, ba, pos);
        return true;
    }

    bool write(QByteArray& ba) {

        writePOD<int>(a, ba);
        writePOD<bool>(bottom, ba);
        writePOD<double>(invMass, ba);
        writePOD<bool>(stiffRefVertex, ba);
        return true;
    }

    bool operator==(const VerticeData& other) const { //
        return false;
    }

};


///////////////////////////////////////////////////////////////////////////////

// Read/write Cell Data
bool inline readAttr(FaceData& m, const QByteArray& ba, size_t& pos)
{
  return mdx::readChar((char*)&m, sizeof(FaceData), ba, pos);
}
bool inline writeAttr(const FaceData& m, QByteArray& ba)
{
  return mdx::writeChar((char*)&m, sizeof(FaceData), ba);
}
bool inline readAttr(EdgeData& m, const QByteArray& ba, size_t& pos)
{
  return mdx::readChar((char*)&m, sizeof(EdgeData), ba, pos);
}
bool inline writeAttr(const EdgeData& m, QByteArray& ba)
{
  return mdx::writeChar((char*)&m, sizeof(EdgeData), ba);
}
bool inline readAttr(VerticeData& m, const QByteArray& ba, size_t& pos)
{
  return mdx::readChar((char*)&m, sizeof(VerticeData), ba, pos);
}
bool inline writeAttr(const VerticeData& m, QByteArray& ba)
{
  return mdx::writeChar((char*)&m, sizeof(VerticeData), ba);
}

//////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////


class Chemics : public Process {
public:
    Chemics(const Process& process)
        : Process(process) {
        setName("Model/Folder/30 Chemics");

        //Auxin
        addParm("Auxin", "Auxin", "");
        addParm("Auxin Basal Production Rate", "", "0.1");
        addParm("Auxin Decay Rate", "Auxin Decay Rate", "0.01");
        addParm("Auxin Permeability", "Auxin Permeability", "0.1");
        addParm("Pin1-auxin export rate", "Pin1-auxin export rate", "0.1");
        //PIN1
        addParm("Pin1", "Pin1", "");
        addParm("Pin1 Basal Production Rate", "Pin1 Basal Production Rate", "0.2");
        addParm("Pin1 Max Auxin-induced Expression", "Pin1 Max Auxin-induced Expression", "50");
        addParm("Pin1 Half-max Auxin-induced K", "Pin1 Half-max Auxin-induced K", "0.01");
        addParm("Pin1 Cytoplasmic Max Concentration", "Pin1 Cytoplasmic Max Concentration (auxin per nm squared)", "2");
        addParm("Pin1 Cytoplasmic Decay Rate", "Pin1 Cytoplasmic Decay Rate", "0.08");
        addParm("Pin1 Max Trafficking Rate", "Pin1 Max Trafficking Rate", "0.75");
        addParm("Pin1 Membrane Max Concentration", "Pin1 Membrane Max Concentration(auxin per nm)", "20");
        addParm("Pin1 Membrane Max Decay Rate", "Pin1 Membrane Max Decay Rate (in satatic decay rate)", "0.08");
        //Source
        addParm("b0 source", "Source Auxin Bottom", "2");
        addParm("b1 source", "Source Auxin Bottom", "30");
        addParm("b2 source", "Source Auxin Bottom", "3.8");
        addParm("b3 source", "Source Auxin Bottom", "3");
        addParm("t0 source", "Source Auxin Top", "4");
        addParm("t1 source", "Source Auxin Top", "50");
        addParm("t2 source", "Source Auxin Top", "5");
        addParm("t3 source", "Source Auxin Top", "5");

    }
    void calcDerivsCell(const CCStructure& cs,
                           CCIndexDataAttr& indexAttr,
                           FaceDataAttr& faceAttr,
                           EdgeDataAttr& edgeAttr,
                           CCIndex f, double Dt,
                           double& prodAD,double& prodAB,double& prodADL1,double& prodABL1,
                           double& pinAD,double& pinAB,double& pinADL1,double& pinABL1,
                           double& finalAD,double& finalAB,double& finalADL1,double& finalABL1,
                           double& excessAD,double& excessAB,double& excessMD,
                           double& auxBasalAD, double& auxDecayAD, double& auxDiffAD, double& auxTransAD,  double& auxBasalAB, double& auxDecayAB, double& auxDiffAB, double& auxTransAB,
                           int stepCount,
                           double& topAuxin, double& bottomAuxin,
                           double& prodMD, double& prodMDL1, double& finalMD, double& finalMDL1, double& pinMD, double& pinMDL1, double& auxBasalMD,
                           double& auxDecayMD, double& auxDiffMD, double& auxTransMD);

    bool initialize(QWidget* parent);


    //AUXIN
    double permeability;
    double Kbase;
    double Kdecay;
    double auxinMaxCell;
    double auxinMaxEdge;
    double Kpin1;
    //PIN1
    double pin1BasalProduction;
    double pin1MaxInd;
    double pin1Kaux;
    double pin1MaxCyt;
    double pin1CytDecay;
    double pin1MaxTraffic;
    double pin1MaxEdge;
    double pin1MemDecay;
    //Source
    double b0, b1, b2, b3;
    double t0, t1, t2, t3;

private:
    Mesh* mesh = 0;
};

  class ChemicalSolver : public Process {
  public:
      ChemicalSolver(const Process& process)
          : Process(process) {
          setName("Model/Folder/31 Chemical Solver");
          addParm("Dt", "", "0.01");
          addParm( "Strain Measure Method",  "Strain Measure Method", "Cell Axis",
                                                                        QStringList() << "Given"
                                                                                      << "Cell Axis");

          addParm("Pin1", "Pin1Sensivity", "");//PIN1
          addParm("Auxin-Flux Impact Half-max", "Auxin-Flux Impact Half-max", "0.1");
          addParm("MF Impact Half-max","MF Impact Half-max", "0.5");
          addParm("Auxin Max Concentration Cell", "Auxin Max Amount (auxin per nm squared)", "3"); //should be the same as in Chemics
          addParm("Pin1 Sensitivity Average Method", "Pin1 Sensitivity Average Method", "Arithmetic Average",
                                                                        QStringList() << "Soft-max"
                                                                                      << "Arithmetic Average"
                                                                                      << "None");
          addParm("Chemics Process","Name of Chemicals derivatives process","Model/Folder/30 Chemics");
          addParm("Debug","Debug","True",QStringList() << "False"
                                                       << "True");

      }

      bool initialize(QWidget* parent,std::string debFileName);
      bool step(bool timeToDebug, int stepCount,double KPin1auxinflux, double KPin1MF, double KPin1MFauxinflux);


      double Dt;
      //PIN1Sensivity
      double KauxinFlux;
      double Kmf;

  private:
      Mesh* mesh = 0;
      Chemics* chemicsProcess = 0;
      //debugging
      double time = 0;
      std::ofstream output_file;
      string snapshotDir;
      int screenShotCount=0;

  };

  ////////////////////////////////////////////////////////////////////////////////////////////



  class Mechanics : public Process{
  public:
      Mechanics(const Process& process)
          : Process(process) {
          setName("Model/Folder/20 Mechanics");

          addParm("Pressure", "", "0.5");
      }

      void calcForces(CCIndex v, const CCStructure& cs, CCIndexDataAttr& indexAttr,
                       FaceDataAttr& faceAttr, EdgeDataAttr& edgeAttr, VerticeDataAttr& verticeAttr);

      bool initialize(QWidget* parent);

      double Pressure;

  private:
      Mesh* mesh = 0;
  };


  class MechanicalSolver : public Process {
  public:
      MechanicalSolver(const Process& process)
          : Process(process) {
          setName("Model/Folder/21 Mechanical Solver");
          addParm("Mechanics Process","", "Model/Folder/20 Mechanics");
          addParm("Dt", "Time step in hours", "0.01"); //same as in restUpdate
          addParm("Velocity Static Damping", "Damping", "0");
          addParm("Velocity Update","Solver Order", "First Order", QStringList() << "First Order"
                                                                                 << "Second Order");
          addParm("Max Velocity", "", "50");

          addParm("n Stiffness", "Distance", "4");
          addParm("K Stiffness", "Distance", "40");
          addParm("CZPZDiscretize", "spring stiffness", "True", QStringList() << "True"
                  << "False");
          addParm("CZStiffness", "Distance", "1");
          addParm("PZStiffness", "Distance", "1");


          addParm("PBD iterations", "PBD iterations", "5");

          addParm("-------------D", "DistanceConstraint", "69");
          addParm("PBD springs stiffness", "stiffness springs", "1");
          addParm("PBD spring active", "spring stiffness", "True", QStringList() << "True"
                                                                                  << "False");
          addParm("Wall eStifness", "Extension Stiffness of the wall springs", "0.1");
          addParm("Auxin-induced wall relaxation", "Auxin-induced wall relaxation", "0.01");
          addParm("Interior eStifness", "Extension Stiffness of the inner springs", "0.01");
          addParm("cStifness", "Compression Stiffness of springs", "1");

          addParm("-------------S", "ShapeConstraint", "69");
          addParm("PBD shape stiffness", "PBD stiffness shape", "0.5");
          addParm("PBD shape allow stretching", "Allow stretching", "True", QStringList() << "False"
                                                                                          << "True");

          addParm("-------------Sz", "BendingConstraint", "69");
          addParm("PBD bending stiffness", "PBD stiffness bending", "0.5");

          addParm("-------------St", "StrainConstraint", "69");
          addParm("PBD strain stiffness", "strain", "0.1");
          addParm("PBD strain XY stiffness", "XY", "0");
          addParm("PBD strain XX stiffness", "XX", "1");
          addParm("PBD strain YY stiffness", "YY", "1");
          addParm("PBD shear normalize", "normalize", "False", QStringList() << "True"
                                                                              << "False");
          addParm("PBD stretch normalize", "normalize", "False", QStringList() << "True"
                                                                              << "False");


          addParm("--------------P","Polarity Update","");
          addParm("Polarity Method", "Polarity Method", "Green Strain Tensor",
                                       QStringList() << "Cell Axis"
                                                     << "Green Strain Tensor");
          addParm("MF Strain Threshold for Growth", "MF Strain Threshold for Growth", "0.01"); //same as spring treshold in update //Marco no usa el mismo
          addParm("MF reorientation rate", "MF Reorientation Rate", "0.04");
          addParm("MF Degradation", "MF Degradation", "0.01");
          addParm("Momentum Preservation K", "", "0");
      }

      void momentumPreservation(const CCStructure& cs, CCIndexDataAttr& indexAttr, VerticeDataAttr& verticeAttr, Point3d P_r, Point3d L_r, double& meanCorrection);

      void solveEuler1stOrder(CCIndex v, const CCStructure& cs, CCIndexDataAttr& indexAttr, VerticeDataAttr& verticeAttr,
                              Point3d& P_r, Point3d& L_r, Point3d centerOfMasses);

      void solveConstraints(double& dStrainADIN, double& dStrainADL1, double& dStrainABIN, double& dStrainABL1, double& dStrainMDIN, double& dStrainMDL1, int& stepCount);

      bool initialize(QWidget* parent, std::string debFileName);
      bool step(bool timeToDebug, int stepCount);

      double Dt;
      double time;
      double staticDamping;
      QString velocityUpdate;
      double maxVelocity;

      int PBDiterations;

      double PBDspringsStiffness;
      bool springActive;
      double WallExtStiff;
      double InExtStiff;
      double ComStiffness;

      double PBDshapeStiffness;
      bool PBDshapeAllowStretching;

      double PBDbendingStiffness;

      double PBDstrainStiffness;
      double PBDstrainXYStiffness;
      double PBDstrainXXStiffness;
      double PBDstrainYYStiffness;
      bool PBDshearNormalize;
      bool PBDstretchNormalize;
      //polarity
      double MFspringTresh;
      double MFRORate;
      double MFdegradation;

  private:
      Mesh* mesh = 0;
      Mechanics* mechanicsProcess = 0;
  };


  class restUpdate : public Process
  {
  public:
    restUpdate(const Process &process): Process(process)
    {
      setName("Model/Folder/22 Update Rest");
      addParm("Dt", "Time step in hours", "0.01"); //same as mechnicalSolver
      addParm("Walls Growth", "Walls Growth", "");
      addParm("Strain Threshold for Growth", "Strain Threshold for Growth", "0.0001"); //same as in MF treshold for reorientation, only reorientate if the edge grows
      addParm("Cell Threshold for Growth", "Cell Threshold for Growth", "0.0001");
      addParm("Walls Growth Rate", "Walls Growth Rate", "10");
      addParm("AUXIN", "AUXIN", "");
      addParm("kL1", "kL1", "0.02");
      addParm("n1L1", "n1L1", "2");
      addParm("kmL1", "kmL1", "25");
      addParm("n2L1", "n2L1", "4");
      addParm("bL1", "bL1", "1");
      addParm("kIN", "kIN", "0.1");
      addParm("n1IN", "n1IN", "1");
      addParm("kmIN", "kmIN", "10");
      addParm("n2IN", "n2IN", "16");
      addParm("bIN", "bIN", "0.9");
      }

    bool initialize(QWidget *parent, std::string debFileName);
    bool step(bool timeToDebug, int stepCount);

    double Dt;
    double  springThresh;
    double  growthRateThresh;
    double  wallsGrowthRate;
    double divL1;
    double n1L1;
    double kL1;
    double n2L1;
    double iniL1;
    double divIN;
    double n1IN;
    double kIN;
    double n2IN;
    double iniIN;

    private:
      Mesh *mesh = 0;
  };


  class MySubdivide : public Process
  {
  public:
    MySubdivide (const Process &process) : Process(process)
    {
      setName("Model/Folder/100 Subdivide");

    }

    bool initialize(QWidget *parent);
    bool step(Dimension dim, CCStructure &cs, CCStructure::SplitStruct &ss,
              CCIndex otherP, CCIndex otherN, double interpPos, bool triangulate, std::vector<CCIndexVec> edgeDivisions, std::vector<edgeParData> edatas,
              double KPin1auxinflux, double KPin1MF, double KPin1MFauxinflux, std::vector<edgeStruct>  edgeData);



  private:
    Mesh *mesh = 0;

  };



  class MyTriangulate : public Process {
  public:
      MyTriangulate(const Process& process): Process(process) {
          setName("Model/Folder/06 Triangulate");
          addParm("Max Area", "Max area for triangles", "100");
          //addParm("First Triangulation", "", "False", QStringList() << "False"
          //                                                          << "True");
          addParm("Subdivide Process", "", "Model/Folder/100 Subdivide");
      }


      bool initialize(QWidget* parent);
      bool step(bool  firstTriangulation);

      double maxArea;
      bool woSplit;
      bool firstTriangulation;


private:
     Mesh *mesh = 0;
     MySubdivide *subdivideProcess = 0;

};



  class MyDivide : public Process {
  public:
    MyDivide(const Process &process) : Process(process) {
        setName("Model/Folder/05 Divide");
        setDesc("");

        addParm("Cell Max Area", "Maximum cell area before division", "50.0");
        addParm("AB Div Size", "AB Div Sizen", "2");
        addParm("Cell Wall Min", "Minimum distance during division to existing vertex", "0.1");
        addParm("Cell Wall Sample", "Cell wall sampling distance", "0.05");
        addParm("Cell Pinch", "Amount to pinch walls", "0.2");
        addParm("Cell Max Pinch", "Maximum absolute amount to pinch walls", "1.0");

        addParm("Subdivide Process", "", "Model/Folder/100 Subdivide");
        addParm("Manual Division", "Activate Manual Division", "False", QStringList() << "True"
                                                                                       << "False");

    }

    bool initialize(QWidget* parent);

    bool step(double KPin1auxinflux, double KPin1MF, double KPin1MFauxinflux);


  private:
    Mesh *mesh = 0;
    double CellMaxArea;
    double ABdiv;
    double CellWallMin;
    double CellWallSample;
    double CellPinch;
    double CellMaxPinch;
    char divAlg;
    bool Manual;

    MySubdivide *subdivideProcess = 0;
    //MyTriangulate *triangulateProcess = 0;

  };



  class SetCellType : public Process {
  public:
      SetCellType(const Process& process): Process(process) {
          setName("Model/Folder/11 Set Cell Type");
          addParm("Set Side", "", "False", QStringList() << "True"
                                                       << "False");
          addParm("Cell Side", "", "UndefSide", QStringList() << "UndefSide"
                                                               << "Left"
                                                               << "Right");

          addParm("Set Layer", "", "False", QStringList() << "True"
                                                       << "False");
          addParm("Cell Layer", "", "UndefLayer", QStringList() << "UndefLayer"
                                                               << "L1"
                                                               << "L2"
                                                               << "L3");

          addParm("Set Zone", "", "False", QStringList() << "True"
                                                       << "False");
          addParm("Cell Zone", "", "UndefZone", QStringList() << "UndefZone"
                                                              << "CZ"
                                                              << "PZ"
                                                              <<"AD"
                                                              <<"AB"
                                                              <<"MD");

          addParm("Set OC", "", "False", QStringList() << "True"
                                                       << "False");
          addParm("Set Position", "", "False", QStringList() << "True"
                                                       << "False");
          addParm("Cell Position", "", "Nothing", QStringList() << "Nothing"
                  <<"1"
                  <<"0.9"
                  <<"0.8"
                  <<"0.7"
                  <<"0.6"
                  <<"0.5"
                  <<"0.4"
                  <<"0.3"
                  <<"0.2"
                  <<"0.1"
                  <<"0.05"
                  <<"0.01");

          addParm("Set Initial Auxin", "", "False", QStringList() << "True"
                                                       << "False");
          addParm("Initial Auxin", "", "1", QStringList()
                  <<"100"
                  <<"80"
                  <<"50"
                  <<"25"
                  <<"10"
                  <<"5"
                  <<"2"
                  <<"1");

          addParm("Index", "","0");
      }
      bool step();
      bool rewind(QWidget *parent);

   private:
      Mesh *mesh = 0;
      bool SetSide;
      bool SetLayer;
      bool SetZone;
      bool SetOC;
      bool SetPosition;
      bool SetInitialAuxin;
      uint CellSide;
      uint CellLayer;
      uint CellZone;
      uint CellPosition;
      uint InitialAuxin;
      int i;
  };

  class SetAttr : public Process {
  public:
      SetAttr(const Process& process): Process(process) {
          setName("Model/Folder/12 Set Attributes");
          addParm("Delete", "", "False", QStringList() << "True"
                                                           << "False");
          addParm("Label", "", "True", QStringList() << "True"
                                                           << "False");
          addParm("StiffRef", "", "True", QStringList() << "True"
                                                           << "False");
          addParm("Bottom", "", "False", QStringList() << "True"
                                                           << "False");
          addParm("HighMass", "", "False", QStringList() << "True"
                                                           << "False");
          addParm("ZPOS", "", "False", QStringList() << "True"
                                                           << "False");
          addParm("EdgeLength", "", "False", QStringList() << "True"
                                                           << "False");
          addParm("Index", "","0");
      }
      bool step();

   private:
      Mesh *mesh = 0;
      int i;
      bool Delete;
      bool Label;
      bool StiffRef;
      bool Bottom;
      bool HighMass;
      bool ZPOS;
      bool EdgeLength;
  };



  class General : public Process
  {
  public:
    General(const Process &process) : Process(process)
    {
      setName("Model/Folder/00 General");
      setDesc("");
      setIcon(QIcon(":/images/CellDisk.png"));

      addParm("SetAttr Process", "", "Model/Folder/12 Set Attributes");
      addParm("MechanicalSolver Process", "", "Model/Folder/21 Mechanical Solver");
      addParm("ChemicalSolver Process", "", "Model/Folder/31 Chemical Solver");
      addParm("restUpdate Process", "", "Model/Folder/22 Update Rest");
      addParm("Divide Process", "", "Model/Folder/05 Divide");
      addParm("Remove Bottom", "", "Model/Folder/07 Remove Bottom");
      addParm("Triangulate Process", "", "Model/Folder/06 Triangulate");
      addParm("Chemical steps", "Chemical steps per one mechanical", "1");//20
      addParm("Max Step Interval", "Max Step Interval", "30");
      addParm("Take Snapshots", "Take snapshots of the simulation","False", QStringList() << "True"
                                                                                          << "False");
      addParm("Save View", "Save Final View Mesh","True", QStringList() << "True"
                                                                         << "False");
      addParm("SaveMesh", "SaveMesh", "Mesh/System/Save");
      addParm("Debug","Debug","True",QStringList() << "False"
                                                   << "True");
      addParm("Simulation", "Simulation Number", "1");
      addParm("Execution Steps", "Execution Steps", "1000");
      addParm("Execution Time", "Execution Time(s)", "3600");

      addParm("Wall eStifness", "Extension Stiffness of the wall springs", "0.1");

      //VisualTissue
      addParm("Draw Strain","Draw Strain","30");
      addParm("Draw Auxin-Flow","Draw Auxin-Flow","20");
      addParm("Draw PINs","Draw PINs","4");
      addParm("Draw Stress","Draw Stress","20");
      addParm("Draw StiffXX","Draw StiffXX","4");

      //Chemical
      addParm("Pin1 Sensitivity Auxin-flux K", "Pin1 Sensitivity Auxin-flux K", "5");
      addParm("Pin1 Sensitivity MF K", "Pin1 Sensitivity MF K", "1");
      addParm("Pin1 Sensitivity MF+Auxin-flux K", "Pin1 Sensitivity MF+Auxin-flux K", "20");

    }

    bool initialize(QWidget *parent);
    bool step();
    bool rewind(QWidget *parent);


    int chemicalSteps;
    int maxStepInterval;
    int executionSteps;
    double executionTime;
    bool timeToDebug=false;
    double drawStrain=0;
    double drawAuxinFlow=0;
    double drawPIN=0;
    double drawStress=0;
    double drawStiffXX=0;

    //chemical
    double KPin1auxinflux=0;
    double KPin1MF=0;
    double KPin1MFauxinflux=0;

  private:
    Mesh *mesh = 0;
    string snapshotDir;
    int screenShotCount=0;
    int stepCount = 0;
    int stepInterval =0;
    double time = 0;
    double realTime = 0;
    std::ofstream output_file;

    SetAttr *setAttrProcess = 0;
    MechanicalSolver *mechanicalProcess = 0;
    ChemicalSolver *chemicalProcess = 0;
    restUpdate *restUpdateProcess = 0;
    MyDivide *divideProcess = 0;
    MyTriangulate *triangulateProcess = 0;
    SaveViewFile* saveViewProcess = 0;
    MeshSave* saveMeshProcess = 0;


  };

}
#endif

