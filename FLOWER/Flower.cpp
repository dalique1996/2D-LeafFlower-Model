#include <Flower.hpp>
#include <CCUtils.hpp>

namespace Flower
{

// Initialize the main solver process
bool General::initialize(QWidget *parent)
{
    system("mkdir \"OUTPUT\""); //Folder for output information

    mesh = currentMesh();
    if(!mesh)
      throw(QString("General::initialize No current mesh"));
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");
    CCStructure& lm = mesh->ccStructure("Limits");

//SETUP PROCESSES

    //Initialize SetAttr
    if(!getProcess(parm("SetAttr Process"), setAttrProcess))
      throw(QString("General::initialize Cannot make SetAttr Process"));
    setAttrProcess->step();

    //Initialize Triangulate
    if(!getProcess(parm("Triangulate Process"), triangulateProcess))
      throw(QString("General::initialize Cannot make triangulate process"));
    triangulateProcess->initialize(parent);
    triangulateProcess->step(true);

//PLOT INITIAL LIMITS MORPHOLOGY (box)

    std::vector<Point3d> box; //rectangle which englobes the exterior of the cell
    Point3d axisMin, axisMax; //shortest and longest axis of the rectangle
    std::vector<Point3d> points;
    for(CCIndex v:cs.vertices())
        points.push_back(indexAttr[v].pos);
    MBR(points, box);
    axisMax = box[0] - box[1];
    axisMin = box[1] - box[2];
    if(norm(axisMax) < norm(axisMin)) {
        Point3d tmp = axisMin;
        axisMin = axisMax;
        axisMax = tmp;
    }

    //plot box
    CCIndex v0 = CCIndexFactory.getIndex();
    CCIndex v1 = CCIndexFactory.getIndex();
    CCIndex v2 = CCIndexFactory.getIndex();
    CCIndex v3 = CCIndexFactory.getIndex();
    CCIndex box_e0 = CCIndexFactory.getIndex();
    CCIndex box_e1 = CCIndexFactory.getIndex();
    CCIndex box_e2 = CCIndexFactory.getIndex();
    CCIndex box_e3 = CCIndexFactory.getIndex();

    lm.addCell(v0);
    lm.addCell(v1);
    lm.addCell(v2);
    lm.addCell(v3);
    lm.addCell(box_e0, +v0 - v1);
    lm.addCell(box_e1, +v1 - v2);
    lm.addCell(box_e2, +v2 - v3);
    lm.addCell(box_e3, +v3 - v0);

    indexAttr[v0].pos = box[0];
    indexAttr[v1].pos = box[1];
    indexAttr[v2].pos = box[2];
    indexAttr[v3].pos = box[3];

//DEBUGGING NAME
    std::string debFileName=parm("Simulation").toStdString();

//INITIALIZE PROCESSES

    //Initialize chemical solver
    if(!getProcess(parm("ChemicalSolver Process"), chemicalProcess))
      throw(QString("General::initialize Cannot make ChemicalSolver Process"));
    chemicalProcess->initialize(parent, debFileName);

    //Initialize mechanical solver
    if(!getProcess(parm("MechanicalSolver Process"), mechanicalProcess))
      throw(QString("General::initialize Cannot make MechanicalSolver Process"));
    mechanicalProcess->initialize(parent, debFileName);

    //Initialize dividing process
    if(!getProcess(parm("Divide Process"), divideProcess))
      throw(QString("General::initialize Cannot make divide process"));
    divideProcess->initialize(parent);

    //Initialize update reststate
    if(!getProcess(parm("restUpdate Process"), restUpdateProcess))
      throw(QString("General::initialize Cannot make restUpdate process"));
    restUpdateProcess->initialize(parent, debFileName);

//SAVE SCREENSHOTS OF THE SIMULATION

    if(!getProcess(parm("SaveMesh"), saveMeshProcess))
      throw(QString("Initialize Cannot make SaveMesh") + parm("SaveMesh"));

    if(parm("Take Snapshots") == "True") {
       system("mkdir \"snapshots\"");
       snapshotDir = "snapshots/" + currentDateTime() + "/";
       if(system((string("mkdir ") + snapshotDir).c_str()) == -1)
         throw("Error creating snapshot directory");
     }

//CREATE HEADER FOR GENERAL DEBUGGING

    if(parm("Debug") == "True") {
       if(!output_file.is_open()) {
         string debugFile = "OUTPUT/G" + debFileName + ".csv";
         output_file.open(debugFile);
         output_file << "Steps" << ","
                     << "Dt"  << ","
                     << "Velocity" << ","
                     << "minArea" << ","
                     << "meanArea" << ","
                     << "maxArea" << ","
                     << "minLength" << ","
                     << "meanLength" << ","
                     << "maxLength" << ","
                     << "AD_AREA" << ","
                     << "LEFT_AREA" << ","
                     << "AB_AREA" << ","
                     << "RIGHT_AREA" << ","
                     << "ADgrowth" << ","
                     << "ABgrowth" << ","
                     << "MDgrowth" << ","
                     << "ADL1growth" << ","
                     << "ABL1growth" << ","
                     << "MDL1growth" << ","
                     << "TOTgrowth" << ","
                     << "ADaniso" << ","
                     << "ABaniso" << ","
                     << "MDaniso" << ","
                     << "ADL1aniso" << ","
                     << "ABL1aniso" << ","
                     << "MDL1aniso" << ","
                     << "TOTaniso" << ","
                     << "ADnumber" << ","
                     << "ADL1number" << ","
                     << "ABnumber" << ","
                     << "ABL1number" << ","
                     << "MDnumber" << ","
                     << "MDL1number" << ","
                     << "LEFTL1_num" << ","
                     << "LEFT_num" << ","
                     << "RIGHTL1_num" << ","
                     << "RIGHT_num" << ","
                     << "TOTnumber" << ","
                     << "DistanceCorr" << ","
                     << "ShapeCorr" << ","
                     << "StrainCorr" << ","
                     << "BendingCorr" << ","
                     << "TotalCorr" << endl
                     << flush;
         }
       }

//READ PARAMETERS

     //Time Parameters
     chemicalSteps=parm("Chemical steps").toInt();
     maxStepInterval = parm("Max Step Interval").toInt();
     executionSteps =parm("Execution Steps").toInt();
     executionTime =parm("Execution Time").toDouble();

     //Draw Parameters
     drawStrain = parm("Draw Strain").toDouble(); //CMTs model predicted orientation
     drawAuxinFlow = parm("Draw Auxin-Flow").toDouble();
     drawPIN = parm("Draw PINs").toDouble();
     drawStress = parm("Draw Stress").toDouble();
     drawStiffXX = parm("Draw StiffXX").toDouble(); //CMTs predefined pattern

     //Chemical Parameters
     KPin1auxinflux = parm("Pin1 Sensitivity Auxin-flux K").toDouble();
     KPin1MF = parm("Pin1 Sensitivity MF K").toDouble();
     KPin1MFauxinflux = parm("Pin1 Sensitivity MF+Auxin-flux K").toDouble();

  return true;
}


bool General::step()
{

//DEBUGGING
    clock_t time_req;
    time_req = clock();

//INITIALIZE
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
    EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");
    VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");
    CCStructure& cs = mesh->ccStructure("Tissue");
    CCStructure& tt = mesh->ccStructure("Triangulate");
    CCStructure& cv = mesh->ccStructure("TissueVisual");
    CCStructure& lm = mesh->ccStructure("Limits");


//TISSUE VISUAL: Update old centroid positionn
    std::vector<int> labelsOld;
    for (CCIndex f:cs.faces()){
        labelsOld.push_back(indexAttr[f].label);
    }

//CLEAR OLD ATTRIBUTES: IMPORTANT: you have to check all the cell complexes
    //Edges
    std::set<CCIndex> old_cells;
    for(auto i : edgeAttr)
        if(!cs.hasCell(i.first) && !tt.hasCell(i.first) && !cv.hasCell(i.first) && !lm.hasCell(i.first))
            old_cells.insert(i.first);
    for(CCIndex e : old_cells)
            edgeAttr.erase(e);
    old_cells.clear();
    //Vertices
    for(auto i : verticeAttr)
        if(!cs.hasCell(i.first) && !tt.hasCell(i.first) && !cv.hasCell(i.first) && !lm.hasCell(i.first))
            old_cells.insert(i.first);
    for(CCIndex v : old_cells)
            verticeAttr.erase(v);
    old_cells.clear();
    //Faces
    for(auto i : faceAttr)
        if(!cs.hasCell(i.first) && !tt.hasCell(i.first) && !cv.hasCell(i.first) && !lm.hasCell(i.first))
            old_cells.insert(i.first);
    for(CCIndex f : old_cells)
            faceAttr.erase(f);
    old_cells.clear();
    //Indexes
    for(auto i : indexAttr)
        if(!cs.hasCell(i.first) && !tt.hasCell(i.first) && !cv.hasCell(i.first) && !lm.hasCell(i.first))
            old_cells.insert(i.first);
    for(CCIndex i : old_cells)
            indexAttr.erase(i);

//CLEAR SELECTION (in division for triangulation)
    for(CCIndex f: cs.faces())
        indexAttr[f].selected = false;
    for(CCIndex e: cs.edges())
        indexAttr[e].selected = false;
    for(CCIndex v: cs.vertices())
        indexAttr[v].selected = false;
    for(CCIndex f: tt.faces())
        indexAttr[f].selected = false;
    for(CCIndex e: tt.edges())
        indexAttr[e].selected = false;
    for(CCIndex v: tt.vertices())
        indexAttr[v].selected = false;


//CALCULATE TIME_TO_DEBUG
    if(stepInterval==maxStepInterval || stepCount==0){
        timeToDebug=true;
        stepInterval=0; //reset
    }
    else{
       timeToDebug=false;
       stepInterval ++;
    }
    mdxInfo<<"####STEP "<<stepCount<<endl;

//CALL MECHANICAL SOLVER
    mechanicalProcess->step(timeToDebug, stepCount);

//DEBUG PARAMETERS: BEFORE DIVISION
    double avg_velocity=0;
    double min_area=999999999, mean_area=0, max_area=0, AD_area=0, AB_area=0, MDleft_area=0, MDright_area=0;
    int face_count=0;
    double min_length=99999999, mean_length=0, max_length=0;
    int edge_count=0;
    double AD_mean_anisotropy=0, ADL1_mean_anisotropy=0;//anisotropy
    double AB_mean_anisotropy=0, ABL1_mean_anisotropy=0;
    double MD_mean_anisotropy=0, MDL1_mean_anisotropy=0;
    double total_mean_anisotropy=0;
    double AD_mean_growthRate=0, ADL1_mean_growthRate=0;//growthRate
    double AB_mean_growthRate=0, ABL1_mean_growthRate=0;
    double MD_mean_growthRate=0, MDL1_mean_growthRate=0;
    double total_mean_growthRate=0;
    int ADcount=0, ABcount=0, MDcount=0, LEFTcount=0, RIGHTcount=0; //cellNumber before division
    int ADL1count=0, ABL1count=0, MDL1count=0, LEFTL1count=0, RIGHTL1count=0; //cellNumber in L1
    double avg_totalCorrection=0, avg_distanceCorrection=0, avg_shapeCorrection=0, avg_strainCorrection=0, avg_bendingCorrection=0;

//AFTER DIVISION
    int ADnumber=0, ABnumber=0, MDnumber=0; //cellNumber after division (output number)
    int ADL1number=0, ABL1number=0, MDL1number=0;


//DEBUG BEFORE DIVISION
    if(timeToDebug){
        for (CCIndex v:tt.vertices())
            avg_velocity+=norm(verticeAttr[v].velocity);
        avg_velocity/=tt.vertices().size();

        //EdgeLength
        for(CCIndex e:cs.edges()){
            edge_count ++;
            std::pair<CCIndex, CCIndex> endpoints=cs.edgeBounds(e);
            CCIndexData& v1D=indexAttr[endpoints.first];
            CCIndexData& v2D=indexAttr[endpoints.second];
            double edgeLength=0;
            edgeLength=sqrt(pow(v2D.pos[0]-v1D.pos[0],2)+pow(v2D.pos[1]-v1D.pos[1],2));
            mean_length+=edgeLength;
            if(min_length>edgeLength)
               min_length=edgeLength;
            if(max_length<edgeLength)
               max_length=edgeLength;
        }

        //FaceArea
        for(CCIndex f:cs.faces()){
            FaceData& fD=faceAttr[f];

            if(fD.cellZone==AD || fD.cellZone==AB || fD.cellZone==MD){
              face_count ++;

              mean_area +=indexAttr[f].measure;
              if(min_area>indexAttr[f].measure)
                 min_area=indexAttr[f].measure;
              if(max_area<indexAttr[f].measure)
                 max_area=indexAttr[f].measure;

              //Growth; Anisotropy; Cell Number
              double cell_anisotropy=norm(fD.a1)/norm(fD.a2);
              double cell_growthRate=indexAttr[f].measure/fD.prevArea;
              total_mean_anisotropy +=cell_anisotropy;
              total_mean_growthRate +=cell_growthRate;

              if(fD.cellZone==AD){   //Adaxial
                AD_area +=indexAttr[f].measure;
                ADcount ++;
                LEFTcount ++;
                AD_mean_growthRate +=cell_growthRate;
                AD_mean_anisotropy +=cell_anisotropy;
                if(fD.cellLayer==L1){
                   ADL1count ++;
                   LEFTL1count ++;
                   ADL1_mean_growthRate +=cell_growthRate;
                   ADL1_mean_anisotropy +=cell_anisotropy;
                }
                } else if (fD.cellZone==AB){    //Abaxial
                   AB_area +=indexAttr[f].measure;
                   ABcount ++;
                   RIGHTcount ++;
                   AB_mean_growthRate +=cell_growthRate;
                   AB_mean_anisotropy +=cell_anisotropy;
                   if(fD.cellLayer==L1){
                      ABL1count ++;
                      RIGHTL1count ++;
                      ABL1_mean_growthRate +=cell_growthRate;
                      ABL1_mean_anisotropy +=cell_anisotropy;
                    }
                } else if (fD.cellZone==MD){     //Middle
                   if(fD.cellSide==left){ //Area and cell number per left/right zone
                      MDleft_area +=indexAttr[f].measure;
                      LEFTcount ++;
                      if (fD.cellLayer==L1){
                         LEFTL1count ++;
                      }
                   } else if(fD.cellSide==right){
                      MDright_area +=indexAttr[f].measure;
                      RIGHTcount ++;
                      if (fD.cellLayer==L1){
                         RIGHTL1count ++;
                      }
                   }
                   MDcount ++;
                   MD_mean_growthRate +=cell_growthRate;
                   MD_mean_anisotropy +=cell_anisotropy;
                   if(fD.cellLayer==L1){
                     MDL1count ++;
                     MDL1_mean_growthRate +=cell_growthRate;
                     MDL1_mean_anisotropy +=cell_anisotropy;
                   }
                }
            }
        }
        AD_mean_anisotropy/=ADcount;    ADL1_mean_anisotropy/=ADL1count;
        AB_mean_anisotropy/=ABcount;    ABL1_mean_anisotropy/=ABL1count;
        MD_mean_anisotropy/=MDcount;    MDL1_mean_anisotropy/=MDL1count;
        total_mean_anisotropy/=face_count;
        AD_mean_growthRate/=ADcount;   ADL1_mean_growthRate/=ADL1count;
        AB_mean_growthRate/=ABcount;   ABL1_mean_growthRate/=ABL1count;
        MD_mean_growthRate/=MDcount;   MDL1_mean_growthRate/=MDL1count;
        total_mean_growthRate/=face_count;


        //Constraints Corrections Info
        for(CCIndex v:cs.vertices()){
           for(auto deltap : verticeAttr[v].corrections) {
              const char* constraint = deltap.first;
              if(QString(constraint) == "distance"){
                Point3d corr = deltap.second;
                double ncorr= norm(corr);
                avg_distanceCorrection += ncorr;
                avg_totalCorrection += ncorr;
              } else if(QString(constraint) == "shape"){
                Point3d corr = deltap.second;
                double ncorr= norm(corr);
                avg_shapeCorrection += ncorr;
                avg_totalCorrection += ncorr;
              } else if(QString(constraint) == "strain"){
                Point3d corr = deltap.second;
                double ncorr= norm(corr);
                avg_strainCorrection += ncorr;
                avg_totalCorrection += ncorr;
              } else if(QString(constraint) == "bending"){
                Point3d corr = deltap.second;
                double ncorr= norm(corr);
                avg_bendingCorrection += ncorr;
                avg_totalCorrection += ncorr;
               } else
                mdxInfo<<"DEBUG: constraint unknown"<<endl;
             }
        }

        avg_distanceCorrection/=cs.vertices().size();
        avg_shapeCorrection/=cs.vertices().size();
        avg_strainCorrection/=cs.vertices().size();
        avg_bendingCorrection/=cs.vertices().size();
        avg_totalCorrection/=(cs.vertices().size()*3);
    }


//CALL DIVIDE PROCESS
    divideProcess->step(KPin1auxinflux, KPin1MF, KPin1MFauxinflux);

//CALL TRIANGULATION PROCES
    for(CCIndex f:cs.faces()){
       if(indexAttr[f].selected){ //required only if division happend (face selected in division)
         triangulateProcess->step(false);
         break;
       }
    }

//CALL UPDATE RESTSTATE (for constraints)
    restUpdateProcess->step(timeToDebug, stepCount); //Update restState --> distance: restLength; shape: restCm

//UPDATE CENTROID (FACE) AND MIDPOINT (EDGE) -->after division required for chemical solver
    for (CCIndex f: cs.faces()){
        FaceData& fD=faceAttr[f];
        std::vector<CCIndex> indices=faceVertices(cs,f);
        int size=indices.size();
        std::vector<Point3d> vertices;
        for (CCIndex v: indices){
            vertices.push_back(indexAttr[v].pos);
        }
        Point3d centroid=Centroid(vertices,size);
        fD.newcentroid=centroid;
    }

    for (CCIndex e:cs.edges()){
        EdgeData& eD=edgeAttr[e];
        auto eb=cs.edgeBounds(e);
        eD.midPoint= indexAttr[eb.first].pos+(indexAttr[eb.second].pos-indexAttr[eb.first].pos)/2.;
    }

//CALL CHEMICAL SOLVER
    for(int i=0; i<chemicalSteps; i++)
       chemicalProcess->step(timeToDebug, stepCount, KPin1auxinflux, KPin1MF, KPin1MFauxinflux);

    time += mechanicalProcess->Dt; //Dt for debugging

//DEBUG AFTER DIVISION & SIGNAL MAPS
    if(timeToDebug){
      cv.clear();
      //Signal Maps
      CCIndexIntAttr& celllayer = mesh->signalAttr<int>("CellLayer");
      celllayer.clear();
      CCIndexIntAttr& cellzone = mesh->signalAttr<int>("CellZone");
      cellzone.clear();
      CCIndexIntAttr& cellposition = mesh->signalAttr<int>("CellPosition"); //Identifier for initial cells to set different attr
      cellposition.clear();
      CCIndexDoubleAttr& growthRate = mesh->signalAttr<double>("Growth Rate"); //Change in area between steps
      growthRate.clear();
      CCIndexDoubleAttr& anisotropy = mesh->signalAttr<double>("Anisotropy"); //MaximalGrowth/MinimalGrowth
      anisotropy.clear();
      CCIndexDoubleAttr& cellstiffness = mesh->signalAttr<double>("CellStiffness"); //Mean stiffness of the edges of a cell
      cellstiffness.clear();

      for (CCIndex f: cs.faces()){
          FaceData& fD=faceAttr[f];
          cellzone[f] = fD.cellZone; //CellZone

          if(fD.cellZone ==AD || fD.cellZone ==MD || fD.cellZone ==AB){
            CCIndexData& iD =indexAttr[f];
            double cell_stiffness=0;
            uint i=0;
            for(CCIndex e: cs.incidentCells(f,1)){
               i++;
               CCIndex et =edgeAttr[e].sameEdge;
               cell_stiffness +=edgeAttr[et].stiffness;
            }
            cell_stiffness/=i;
            cellstiffness[f] = cell_stiffness; //CellStiffness
            celllayer[f] = fD.cellLayer; //CellLayer
            cellposition[f] = fD.cellPosition; //CellPosition
            growthRate[f] = iD.measure/fD.prevArea; //Growth Rate

            fD.prevArea= iD.measure; //UPDATE PREV AREA

            anisotropy[f]=norm(fD.a1)/norm(fD.a2); //Anisotropy

            //centroid vertex
            CCIndex v2 = CCIndexFactory.getIndex();
            cv.addCell(v2);
            indexAttr[v2].pos=fD.newcentroid;

            //DRAW PRINCIPAL GROWTH DIRECTIONS
            if(drawStrain > 0) {
              CCIndex v100 = CCIndexFactory.getIndex();
              CCIndex v101 = CCIndexFactory.getIndex();
              CCIndex v102 = CCIndexFactory.getIndex();
              CCIndex v103 = CCIndexFactory.getIndex();
              CCIndex v110 = CCIndexFactory.getIndex();
              CCIndex v111 = CCIndexFactory.getIndex();
              CCIndex v112 = CCIndexFactory.getIndex();
              CCIndex v113 = CCIndexFactory.getIndex();
              cv.addCell(v100);
              cv.addCell(v101);
              cv.addCell(v102);
              cv.addCell(v103);
              cv.addCell(v110);
              cv.addCell(v111);
              cv.addCell(v112);
              cv.addCell(v113);
              indexAttr[v100].pos=indexAttr[v2].pos+(Point3d(-fD.a1[1],fD.a1[0],0)*(drawStrain/16))+(fD.a1*drawStrain)/2;
              indexAttr[v101].pos=indexAttr[v2].pos+(Point3d(fD.a1[1],-fD.a1[0],0)*(drawStrain/16))+(fD.a1*drawStrain)/2;
              indexAttr[v110].pos=indexAttr[v2].pos+(Point3d(-fD.a2[1],fD.a2[0],0)*(drawStrain/16))+(fD.a2*drawStrain)/2;
              indexAttr[v111].pos=indexAttr[v2].pos+(Point3d(fD.a2[1],-fD.a2[0],0)*(drawStrain/16))+(fD.a2*drawStrain)/2;
              indexAttr[v103].pos=indexAttr[v100].pos-(fD.a1*drawStrain);
              indexAttr[v102].pos=indexAttr[v101].pos-(fD.a1*drawStrain);
              indexAttr[v113].pos=indexAttr[v110].pos-(fD.a2*drawStrain);
              indexAttr[v112].pos=indexAttr[v111].pos-(fD.a2*drawStrain);
              std::vector<CCIndex> dertices;
              dertices.push_back(v100);
              dertices.push_back(v103);
              dertices.push_back(v102);
              dertices.push_back(v101);
              std::vector<CCIndex> dertices2;
              dertices2.push_back(v110);
              dertices2.push_back(v113);
              dertices2.push_back(v112);
              dertices2.push_back(v111);
              CCIndex d1 = CCIndexFactory.getIndex();
              CCIndex d2 = CCIndexFactory.getIndex();
              addFace(cv, d1,  dertices);
              addFace(cv, d2,  dertices2);
              indexAttr[d1].label = 12;
              indexAttr[d2].label = 12;
            }

            //DRAW STRESS (only for L1)
            if(drawStress > 0) {
              CCIndex v7 = CCIndexFactory.getIndex();
              cv.addCell(v7);
              indexAttr[v7].pos=indexAttr[v2].pos+((fD.stress/norm(fD.stress))*drawStress);
              CCIndex e7 = CCIndexFactory.getIndex();
              cv.addCell(e7 , -v2 +v7);
            }

            //DRAW STRAIN CONSTRAINT DIRECTIONS (predicted by the model or predefined)
            if(drawStiffXX > 0) {
              CCIndex v8 = CCIndexFactory.getIndex();
              cv.addCell(v8);
              double stiffXX=0;
              for(CCIndex t:fD.cellFaces){
                  stiffXX=faceAttr[t].stiffnessXX;
                  break;
              }
              indexAttr[v8].pos=indexAttr[v2].pos+((Rotate90(fD.stress)/norm(fD.stress))*stiffXX*drawStiffXX);
              CCIndex e8 = CCIndexFactory.getIndex();
              cv.addCell(e8 , -v2 +v8);
            }

            //DRAW MAXIMAL AUXIN FLOW DIRECTION
            if(drawAuxinFlow > 0) {
              CCIndex v5 = CCIndexFactory.getIndex();
              CCIndex v6 = CCIndexFactory.getIndex();
              CCIndex v7 = CCIndexFactory.getIndex();
              CCIndex v8 = CCIndexFactory.getIndex();
              CCIndex v9 = CCIndexFactory.getIndex();
              CCIndex v99 = CCIndexFactory.getIndex();
              CCIndex v98 = CCIndexFactory.getIndex();
              CCIndex v97 = CCIndexFactory.getIndex();
              cv.addCell(v5);
              cv.addCell(v6);
              cv.addCell(v7);
              cv.addCell(v8);
              cv.addCell(v9);
              cv.addCell(v99);
              cv.addCell(v98);
              cv.addCell(v97);
              indexAttr[v9].pos=indexAttr[v2].pos+(Point3d(-fD.auxinFluxVector[1],fD.auxinFluxVector[0],0)*(drawAuxinFlow/16));
              indexAttr[v99].pos=indexAttr[v2].pos+(Point3d(fD.auxinFluxVector[1],-fD.auxinFluxVector[0],0)*(drawAuxinFlow/16));
              indexAttr[v97].pos=indexAttr[v9].pos+(fD.auxinFluxVector*drawAuxinFlow); //4 times more than the arrow connections
              indexAttr[v98].pos=indexAttr[v99].pos+(fD.auxinFluxVector*drawAuxinFlow);
              indexAttr[v5].pos=indexAttr[v2].pos+(fD.auxinFluxVector*drawAuxinFlow);
              indexAttr[v6].pos=indexAttr[v2].pos+(fD.auxinFluxVector*(drawAuxinFlow+(drawAuxinFlow/4)));
              indexAttr[v7].pos=indexAttr[v98].pos+(Point3d(-fD.auxinFluxVector[1],fD.auxinFluxVector[0],0)*(drawAuxinFlow/4));
              indexAttr[v8].pos=indexAttr[v5].pos+(Point3d(fD.auxinFluxVector[1],-fD.auxinFluxVector[0],0)*(drawAuxinFlow/4));
              std::vector<CCIndex> vertices;
              vertices.push_back(v8);
              vertices.push_back(v6);
              vertices.push_back(v7);
              vertices.push_back(v5);
              CCIndex h = CCIndexFactory.getIndex();
              addFace(cv, h,  vertices);
              indexAttr[h].label = 12; //label to color
              std::vector<CCIndex> vertices2;
              vertices2.push_back(v9);
              vertices2.push_back(v99);
              vertices2.push_back(v98);
              vertices2.push_back(v97);
              CCIndex h2 = CCIndexFactory.getIndex();
              addFace(cv, h2,  vertices2);
              indexAttr[h2].label = 12;
            }
          }
      }


      //DRAW PIN1 LOCALIZED IN THE MEMBRANE
      if(parm("Draw PINs").toDouble() > 0) {
      double shift = 0.1;
      double extension = parm("Draw PINs").toDouble();

       for (CCIndex e: cs.edges()){
           EdgeData& eD=edgeAttr[e];
           auto eb = cs.edgeBounds(e);
           CCIndex v1 = eb.first;
           CCIndex v2 = eb.second;
           Point3d versor = indexAttr[v2].pos - indexAttr[v1].pos;
           Point3d direction = Rotate90(versor);
           versor /= norm(versor);

           std::vector<CCIndex> fs;
           for(auto f : cs.incidentCells(e, 2))
               fs.push_back(f);
           Point3d insideFace0=faceAttr[fs[0]].newcentroid-eD.midPoint;
           double angle0 =acos((insideFace0*direction)/(norm(insideFace0)*norm(direction)));
           if (angle0<(M_PI/2))// always pointing inside the cell
              direction *=-1;
           Point3d c = direction/norm(direction);
           int label = indexAttr[fs[0]].label;
           double pin = eD.pin1[label] + EPS;
           Point3d p1 = indexAttr[v1].pos + -c * shift + versor * shift;
           Point3d p2 = indexAttr[v2].pos + -c * shift - versor * shift;
           Point3d p4 = p1 + (-c * extension * pin) + versor * shift*3;
           Point3d p3 = p2 + (-c * extension  * pin) - versor * shift*3;

           CCIndex PIN10 = CCIndexFactory.getIndex();
           CCIndex PIN11 = CCIndexFactory.getIndex();
           CCIndex PIN12 = CCIndexFactory.getIndex();
           CCIndex PIN13 = CCIndexFactory.getIndex();
           cv.addCell(PIN10);
           cv.addCell(PIN11);
           cv.addCell(PIN12);
           cv.addCell(PIN13);
           indexAttr[PIN10].pos = p1+ Point3d(0,0,0.04);
           indexAttr[PIN11].pos = p2+ Point3d(0,0,0.04);
           indexAttr[PIN12].pos = p3+ Point3d(0,0,0.04);
           indexAttr[PIN13].pos = p4+ Point3d(0,0,0.04);
           CCIndex PIN30 =CCIndexFactory.getIndex();
           CCIndex PIN31 =CCIndexFactory.getIndex();
           CCIndex PIN32 =CCIndexFactory.getIndex();
           CCIndex PIN33 =CCIndexFactory.getIndex();
           cv.addCell(PIN30,-PIN10 +PIN11);
           cv.addCell(PIN31,-PIN11 +PIN12);
           cv.addCell(PIN32,-PIN12 +PIN13);
           cv.addCell(PIN33,-PIN13 +PIN10);
           CCIndex PIN90 = CCIndexFactory.getIndex();
           cv.addCell(PIN90, -PIN30 -PIN31 -PIN32 -PIN33);
           indexAttr[PIN90].label = 13; //label to color

           if(fs.size() > 1) { //for the second faced delimited by the same edge
             Point3d insideFace1=faceAttr[fs[1]].newcentroid-eD.midPoint;
             double angle1 =acos((insideFace1*direction)/(norm(insideFace1)*norm(direction)));
             if (angle1<(M_PI/2))// always pointing inside the cell
                direction *=-1;
             c = direction/norm(direction);
             label =  indexAttr[fs[1]].label;
             pin = eD.pin1[label] + EPS;
             p1 = indexAttr[v1].pos + -c * shift + versor * shift;
             p2 = indexAttr[v2].pos + -c * shift - versor * shift;
             p4 = p1 + (-c * extension * pin) + versor * shift*3;
             p3 = p2 + (-c * extension * pin) - versor * shift*3;

             CCIndex PIN20 = CCIndexFactory.getIndex();
             CCIndex PIN21 = CCIndexFactory.getIndex();
             CCIndex PIN22 = CCIndexFactory.getIndex();
             CCIndex PIN23 = CCIndexFactory.getIndex();
             cv.addCell(PIN20);
             cv.addCell(PIN21);
             cv.addCell(PIN22);
             cv.addCell(PIN23);
             indexAttr[PIN20].pos = p1+ Point3d(0,0,0.04);
             indexAttr[PIN21].pos = p2+ Point3d(0,0,0.04);
             indexAttr[PIN22].pos = p3+ Point3d(0,0,0.04);
             indexAttr[PIN23].pos = p4+ Point3d(0,0,0.04);
             CCIndex PIN40 =CCIndexFactory.getIndex();
             CCIndex PIN41 =CCIndexFactory.getIndex();
             CCIndex PIN42 =CCIndexFactory.getIndex();
             CCIndex PIN43 =CCIndexFactory.getIndex();
             cv.addCell(PIN40,-PIN20 +PIN21);
             cv.addCell(PIN41,-PIN21 +PIN22);
             cv.addCell(PIN42,-PIN22 +PIN23);
             cv.addCell(PIN43,-PIN23 +PIN20);
             CCIndex PIN91 = CCIndexFactory.getIndex();
             cv.addCell(PIN91, +PIN40 +PIN41 +PIN42 +PIN43);
             indexAttr[PIN91].label = 13; //label to color
          }
       }
      }


    //SAVE A SNAPSHOT
      if(parm("Take Snapshots") == "True"){
        QString fileName = QString::fromStdString(snapshotDir) + QString("Root-%1.png").arg(screenShotCount++, 4, 10, QChar('0'));
        takeSnapshot(fileName, 1, 1280*2, 720*2, 100, true);
      }

    //DEBUG CELL NUMBER AFTER DIVISION
      for(CCIndex f:cs.faces()){
         FaceData& fD=faceAttr[f];
         if(fD.cellZone==AD){
           ADnumber ++;
           if(fD.cellLayer==L1)
              ADL1number ++;
          } else if (fD.cellZone==AB){
              ABnumber ++;
              if(fD.cellLayer==L1)
                ABL1number ++;
          } else if (fD.cellZone==MD){
              MDnumber ++;
              if(fD.cellLayer==L1)
                MDL1number ++;
          }
      }

    //DEBUG PRINTED
      if(parm("Debug") == "True"){
        output_file << stepCount << ","
                    << time << ","
                    << avg_velocity << ","
                    << min_area<< ","
                    << mean_area/face_count<< ","
                    << max_area<< ","
                    << min_length<< ","
                    << mean_length/edge_count<< ","
                    << max_length<< ","
                    << AD_area<<","
                    << AD_area+MDleft_area<<","
                    << AB_area<<","
                    << AB_area+MDright_area<<","
                    << AD_mean_growthRate << ","
                    << AB_mean_growthRate << ","
                    << MD_mean_growthRate << ","
                    << ADL1_mean_growthRate << ","
                    << ABL1_mean_growthRate << ","
                    << MDL1_mean_growthRate << ","
                    << total_mean_growthRate << ","
                    << AD_mean_anisotropy << ","
                    << AB_mean_anisotropy << ","
                    << MD_mean_anisotropy << ","
                    << ADL1_mean_anisotropy << ","
                    << ABL1_mean_anisotropy << ","
                    << MDL1_mean_anisotropy << ","
                    << total_mean_anisotropy << ","
                    << ADnumber << ","
                    << ADL1number << ","
                    << ABnumber << ","
                    << ABL1number << ","
                    << MDnumber << ","
                    << MDL1number << ","
                    << LEFTL1count <<","
                    << LEFTcount <<","
                    << RIGHTL1count << ","
                    << RIGHTcount << ","
                    << face_count << ","
                    << avg_distanceCorrection << ","
                    << avg_shapeCorrection << ","
                    << avg_strainCorrection << ","
                    << avg_bendingCorrection << ","
                    << avg_totalCorrection << endl
                    << flush;

        }

//UPDATE MESH VISUALIZATION
    mesh->updateAll();
    }

//EXIT RUN
    if(stepCount>executionSteps || realTime>executionTime){ //Max number of steps or Max execution time reached
       if(parm("Save View") == "True") {
          string meshString = "OUTPUT/" +  parm("Simulation").toStdString() + ".mdxm";
          QString meshFile = QString::fromStdString(meshString);
          saveMeshProcess->run(mesh, meshFile, false);
       }
       exit(0);
    }
    stepCount ++;
    time_req = clock() - time_req;
    realTime += (float)time_req/CLOCKS_PER_SEC;

  return true;
}

bool General::rewind(QWidget *parent){ //To rewind, we'll reload the mesh
   mesh = currentMesh();
   if(!mesh or mesh->file().isEmpty())
     throw(QString("General::rewind No current mesh, cannot rewind"));
   screenShotCount = 0;
   stepCount = 0;
   stepInterval=0;
   time = 0;
   realTime= 0;
   MeshLoad meshLoad(*this);
   meshLoad.setParm("File Name", mesh->file());
   return meshLoad.run();
}

////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////CHEMICAL SOLVER//////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////

bool Chemics::initialize(QWidget* parent) {
     mesh = getMesh("Mesh 1");
     if(!mesh)
       throw(QString("Chemics::step No current mesh"));

     CCStructure& cs = mesh->ccStructure("Tissue");
     CCIndexDataAttr& indexAttr = mesh->indexAttr();
     FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");

//INITIALIZE PARAMETERS
     //Auxin
      permeability = parm("Auxin Permeability").toDouble();
      Kpin1 = parm("Pin1-auxin export rate").toDouble();
      Kbase = parm("Auxin Basal Production Rate").toDouble();
      Kdecay = parm("Auxin Decay Rate").toDouble();

      //PIN1
      pin1BasalProduction = parm("Pin1 Basal Production Rate").toDouble();
      pin1MaxInd = parm("Pin1 Max Auxin-induced Expression").toDouble();
      pin1Kaux = parm("Pin1 Half-max Auxin-induced K").toDouble();
      pin1MaxCyt = parm("Pin1 Cytoplasmic Max Concentration").toDouble();
      pin1CytDecay = parm("Pin1 Cytoplasmic Decay Rate").toDouble();
      pin1MaxTraffic = parm("Pin1 Max Trafficking Rate").toDouble();
      pin1MaxEdge = parm("Pin1 Membrane Max Concentration").toDouble();
      pin1MemDecay = parm("Pin1 Membrane Max Decay Rate").toDouble();

      //Source
      b0=parm("b0 source").toDouble();
      b1=parm("b1 source").toDouble();
      b2=parm("b2 source").toDouble();
      b3=parm("b3 source").toDouble();
      t0=parm("t0 source").toDouble();
      t1=parm("t1 source").toDouble();
      t2=parm("t2 source").toDouble();
      t3=parm("t3 source").toDouble();

//INITIALIZE AUXIN CONCENTRATION
      for (CCIndex f:cs.faces()){
          FaceData& fD=faceAttr[f];
          double value;
          if(fD.initialAuxin== hundredX)
              value=100;
          else if(fD.initialAuxin==eightyX)
              value=80;
          else if(fD.initialAuxin==fiftyX)
              value=50;
          else if(fD.initialAuxin==twentyfiveX)
              value=25;
          else if(fD.initialAuxin==tenX)
              value=10;
          else if(fD.initialAuxin==fiveX)
              value=5;
          else if(fD.initialAuxin==twoX)
              value=2;
          else if(fD.initialAuxin==oneX)
              value=1;
          else
               throw(QString("Wrong Initial Auxin Value"));
          fD.auxin=(value/140)*indexAttr[f].measure; //multiply per area to use concentration instead of total amount
      }

      return true;
}

////-----------------------------------CHEMICS_DERIVATIVES-----------------------------//

void Chemics::calcDerivsCell(const CCStructure& cs,
                             CCIndexDataAttr& indexAttr,
                             FaceDataAttr& faceAttr,
                             EdgeDataAttr& edgeAttr,
                             CCIndex f, double Dt,
                             double& prodAD,double& prodAB,double& prodADL1,double& prodABL1,
                             double& pinAD,double& pinAB,double& pinADL1,double& pinABL1,
                             double& finalAD,double& finalAB,double& finalADL1,double& finalABL1,
                             double& excessAD, double& excessAB,double& excessMD,
                             double& auxBasalAD, double& auxDecayAD, double& auxDiffAD, double& auxTransAD,  double& auxBasalAB, double& auxDecayAB, double& auxDiffAB, double& auxTransAB,
                             int stepCount,
                             double& topAuxin, double& bottomAuxin,
                             double& prodMD, double& prodMDL1, double& finalMD, double& finalMDL1, double& pinMD, double& pinMDL1, double& auxBasalMD,
                             double& auxDecayMD, double& auxDiffMD, double& auxTransMD) {

//---------------------------------------------------AUXIN-----------------------------//
      FaceData&  fD = faceAttr[f];
      CCIndexData&  ifD = indexAttr[f];
      double diffusion = 0, activeTransport = 0;
      double pin1TrafficPerCell=0;
      int numEdgeIncome=0;


//SOURCE
      if(fD.cellLayer==L1){
        for(CCIndex fn:cs.neighbors(f)){
           if(faceAttr[fn].cellZone==CZ){

             fD.auxinProdRate=t0*(pow(stepCount/(t1*1000),t2))+t3;
             topAuxin=fD.auxinProdRate;
             for(const Flip& flip : cs.matchV(CCIndex::Q, f, fn, CCIndex::TOP))
                numEdgeIncome++;
             break;
           } else if (faceAttr[fn].cellZone==PZ){

             fD.auxinProdRate=b0*(pow(stepCount/(b1*1000),b2))+b3;
             bottomAuxin=fD.auxinProdRate;
             for(const Flip& flip : cs.matchV(CCIndex::Q, f, fn, CCIndex::TOP))
                numEdgeIncome++;
             break;
           }
       }
      } else
        fD.auxinProdRate=0;

//AUXIN DIFFUSION AND TRANSPORT
     std::set<CCIndex> faceEdges =cs.incidentCells(f,1);
     for (CCIndex e:faceEdges){ //for all the edges in the boundary of the cell
         std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);

         if(edgeFaces.size()==1) //No auxin in the borders of the tissue (neither transport nor diffusion)
             continue;
         else if(edgeFaces.size()==2){
             EdgeData& eD = edgeAttr[e];
             CCIndexData&  ieD = indexAttr[e];
             std::vector<CCIndex> edgeFacesV(edgeFaces.begin(), edgeFaces.end());
             CCIndex n; //neigbour cell
             if(edgeFacesV[0]==f)
                n=edgeFacesV[1];
             else if (edgeFacesV[1]==f)
                n=edgeFacesV[0];
             int nlabel=indexAttr[n].label;
             double edge_diffusion;
             if(faceAttr[n].cellZone==PZ || faceAttr[n].cellZone==CZ) //No diffusion to PZ/CZ
               edge_diffusion=0;
             else
               edge_diffusion = permeability * (eD.intercellularAuxin - fD.auxin) * ieD.measure/ifD.measure; //eD.intercellularAuxin is Auxin in the membrane
               //diffusion positive through the interior of the cell                                         // we assure intercell space 1 nm thick

             double edge_import=0; //for now I don't consider input carriers
             edge_import=eD.edge_export_prev[nlabel];//import equal to the neighbour export t-1
             eD.edge_export[ifD.label] = Kpin1 * eD.pin1[ifD.label]*fD.auxin * ieD.measure; //eD.Pin1 is Pin1 in the membrane oriented from the cell to the neighbour
             double edge_transport = edge_import - eD.edge_export[ifD.label];//transport positive: import higher than export

             if(faceAttr[n].cellZone==PZ || faceAttr[n].cellZone==CZ) //No Auxin in the membrane PZ/CZ to simulate a sink (this auxin is lost)
               eD.intercellularAuxin=0;
             else
               eD.intercellularAuxin += (-edge_diffusion -edge_transport -(Kdecay+fD.auxinDecayRate) * eD.intercellularAuxin) * Dt;

             fD.auxinFluxVector += ((fD.newcentroid - eD.midPoint)/norm(fD.newcentroid - eD.midPoint)) * (edge_transport+edge_diffusion) * Dt; //outwards if export higher than import

             if (fD.auxinProdRate>0){//Add the auxin from the source for the auxinFluxVector
                if((faceAttr[n].cellZone==PZ || faceAttr[n].cellZone==CZ) && faceAttr[n].cellLayer==L1)
                  fD.auxinFluxVector +=((fD.newcentroid - eD.midPoint)/norm(fD.newcentroid - eD.midPoint)) * (fD.auxinProdRate/numEdgeIncome) * Dt;
             }

             diffusion += edge_diffusion;
             activeTransport += edge_transport;
      }

      else
          throw("edge with more than two faces attached");
     }

//AUXIN UPDATE
       double basalProduction = Kbase; //production
       double decay = Kdecay * fD.auxin; //decay

       fD.auxin += (basalProduction +fD.auxinProdRate - decay + diffusion + activeTransport) * Dt;

       //Debug
       if(fD.cellZone==AD){
          auxBasalAD +=basalProduction;
          auxDecayAD +=decay;
          auxDiffAD +=diffusion;
          auxTransAD +=activeTransport;
       } else if(fD.cellZone==AB){
          auxBasalAB +=basalProduction;
          auxDecayAB +=decay;
          auxDiffAB +=diffusion;
          auxTransAB +=activeTransport;
       } else if(fD.cellZone==MD){
           auxBasalMD +=basalProduction;
           auxDecayMD +=decay;
           auxDiffMD +=diffusion;
           auxTransMD +=activeTransport;
        }

       if(fD.auxin < 0 || fD.auxin > 1000)
           mdxInfo<<"WARNING: wierd auxin level "<<ifD.label<<endl;

//---------------------------------------------------PIN1-----------------------------//

//CYTOPLASMIC
       double inducedExpression = pin1MaxInd * //by Auxin
                                       (
                                        (pow(fD.auxin / ifD.measure, 2) / (pow(pin1Kaux, 2) + pow(fD.auxin / ifD.measure, 2))) *
                                        (pow(pin1MaxCyt, 8) / (pow(pin1MaxCyt, 8) + pow(fD.pin1 / ifD.measure, 8)))
                                       );

       fD.pin1 += (pin1BasalProduction + inducedExpression - pin1CytDecay*fD.pin1) * Dt;

       //Debug
       if (fD.cellZone==AD){
           prodAD +=fD.pin1;
           if(fD.cellLayer==L1)
              prodADL1 +=fD.pin1;
       } else if (fD.cellZone==AB){
           prodAB +=fD.pin1;
           if(fD.cellLayer==L1)
              prodABL1 +=fD.pin1;
       } else if (fD.cellZone==MD){
           prodMD +=fD.pin1;
           if(fD.cellLayer==L1)
              prodMDL1 +=fD.pin1;
       }


//MEMBRANE
     double pin1Tr = pin1MaxTraffic * (fD.auxin / (fD.auxin + 1)); //depends on Auxin(Higher-->Higer max traficking)

     for (CCIndex e:faceEdges){ //for all the edges in the boundary of the cell
         std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);
         if(edgeFaces.size()==1)
             continue;
         else if(edgeFaces.size()==2){
             EdgeData& eD = edgeAttr[e];
             CCIndexData&  ieD = indexAttr[e];
             std::vector<CCIndex> edgeFacesV(edgeFaces.begin(), edgeFaces.end());
             CCIndex n; //neigbour cell
             if(edgeFacesV[0]==f)
                n=edgeFacesV[1];
             else if (edgeFacesV[1]==f)
                n=edgeFacesV[0];

         eD.traffickedPin1[ifD.label] = fD.pin1* pin1Tr * eD.pin1Sensitivity[ifD.label] *
                                        (pow(pin1MaxEdge,10) / (pow(pin1MaxEdge,10) + pow(eD.pin1[ifD.label]/ieD.measure,10)));

         pin1TrafficPerCell+= eD.traffickedPin1[ifD.label]; //total per cell

         //pin1 update edge
         eD.pin1[ifD.label] += (eD.traffickedPin1[ifD.label] - eD.pin1[ifD.label]*pin1MemDecay) * Dt; //static decay

           //Debug
           if (fD.cellZone==AD){
                pinAD += eD.pin1[ifD.label];
                if(fD.cellLayer==L1)
                  pinADL1 += eD.pin1[ifD.label];
           } else if (fD.cellZone==AB){
                pinAB += eD.pin1[ifD.label];
                if(fD.cellLayer==L1)
                  pinABL1 += eD.pin1[ifD.label];
           } else if (fD.cellZone==MD){
               pinMD += eD.pin1[ifD.label];
               if(fD.cellLayer==L1)
                 pinMDL1 += eD.pin1[ifD.label];

           } else
               throw(QString("edge with more than two faces attached"));
     }
     }

     fD.pin1 -= pin1TrafficPerCell*Dt;

     if(fD.pin1<0){ //limit PIN1 traffic to PIN1 amount in cytoplasm
         double excess=0-fD.pin1;
       //Debug
           if(fD.cellZone==AD)
           excessAD +=excess;
           else if(fD.cellZone==AB)
           excessAB +=excess;
           else if(fD.cellZone==MD)
           excessMD +=excess;

         for (CCIndex e:faceEdges){
             EdgeData& eD = edgeAttr[e];
             fD.pin1=0;
             std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);
             if(edgeFaces.size()==1)
                continue;
             else if(edgeFaces.size()==2){
                double rest=(eD.traffickedPin1[ifD.label]*excess)/pin1TrafficPerCell;
                eD.pin1[ifD.label] -= rest;
                if(eD.pin1[ifD.label]<0)
                  eD.pin1[ifD.label]=0;

                  //Debug
                  if (fD.cellZone==AD){
                       pinAD -= rest;
                       if(fD.cellLayer==L1)
                         pinADL1 -=rest;
                  }
                  else if (fD.cellZone==AB){
                       pinAB -=rest;
                       if(fD.cellLayer==L1)
                          pinABL1 -=rest;
                  }
                  else if (fD.cellZone==MD){
                       pinMD -=rest;
                       if(fD.cellLayer==L1)
                          pinMDL1 -=rest;
                  }
             }
          }
     }

       //Debug
       if (fD.cellZone==AD){
           finalAD += fD.pin1;
           if(fD.cellLayer==L1)
             finalADL1 += fD.pin1;
       } else if (fD.cellZone==AB){
           finalAB += fD.pin1;
           if(fD.cellLayer==L1)
              finalABL1 += fD.pin1;
       } else if (fD.cellZone==MD){
           finalMD += fD.pin1;
           if(fD.cellLayer==L1)
              finalMDL1 += fD.pin1;
       }


}

////-----------------------------------CHEMICAL SOLVER-----------------------------//
bool ChemicalSolver::initialize(QWidget* parent,std::string debFileName) {
     mesh = getMesh("Mesh 1");
     if(!mesh)
        throw(QString("ChemicalSolver::initialize No current mesh"));
      if(!getProcess(parm("Chemics Process"), chemicsProcess))
          throw(QString("ChemicalSolver::initialize Cannot make Chemics Process:") +
                parm("Chemics Process"));
       Dt = parm("Dt").toDouble();

       //PIN1Sensivity
       KauxinFlux = parm("Auxin-Flux Impact Half-max").toDouble();
       Kmf = parm("MF Impact Half-max").toDouble();

       chemicsProcess->initialize(parent);

       if(parm("Take Snapshots") == "True") {
          system("mkdir \"snapshots\"");
          snapshotDir = "snapshots/" + currentDateTime() + "/";
          if(system((string("mkdir ") + snapshotDir).c_str()) == -1)
              throw("Error creating snapshot directory");
       }

//CHEMICAL DEBUG; CREATE THE HEADER
      if(parm("Debug") == "True") {
         if(!output_file.is_open()) {
            string debugFile = "OUTPUT/CH" + debFileName + ".csv";
            output_file.open(debugFile);
            output_file << "steps" << ","
                        << "TOP INCOME" <<","
                        << "BOTTOM INCOME" <<","
                        << "ADinner AUXIN" <<","
                        << "ABinner AUXIN" <<","
                        << "MDinner AUXIN" <<","
                        << "ADL1 AUXIN" <<","
                        << "ABL1 AUXIN" <<","
                        << "MDL1 AUXIN" <<","
                        << "AD PIN1 CYT BEF TRAFFIC" << ","
                        << "AB PIN1 CYT BEF TRAFFIC" << ","
                        << "MD PIN1 CYT BEF TRAFFIC" << ","
                        << "AD L1 PIN1 CYT BEF TRAFFIC" << ","
                        << "AB L1 PIN1 CYT BEF TRAFFIC" << ","
                        << "MD L1 PIN1 CYT BEF TRAFFIC" << ","
                        << "AD PIN1 EDGE" << ","
                        << "AB PIN1 EDGE" << ","
                        << "MD PIN1 EDGE" << ","
                        << "AD L1 PIN1 EDGE" << ","
                        << "AB L1 PIN1 EDGE" << ","
                        << "MD L1 PIN1 EDGE" << ","
                        << "AD PIN1 TRAFFIC" << ","
                        << "AB PIN1 TRAFFIC" << ","
                        << "MD PIN1 TRAFFIC" << ","
                        << "AD L1 TRAFFIC" << ","
                        << "AB L1 TRAFFIC" << ","
                        << "MD L1 TRAFFIC" << ","
                        << "AD PIN1 SENSIVITY" << ","
                        << "AB PIN1 SENSIVITY" << ","
                        << "MD PIN1 SENSIVITY" << ","
                        << "AD L1 PIN1 SENSIVITY" << ","
                        << "AB L1 PIN1 SENSIVITY" << ","
                        << "MD L1 PIN1 SENSIVITY" << ","
                        << "AD PIN1 EXCESS" << ","
                        << "AB PIN1 EXCESS" << ","
                        << "MD PIN1 EXCESS" << ","
                        << "AD AUX BASAL" << ","
                        << "AD AUX DECAY" << ","
                        << "AD AUX DIFF" << ","
                        << "AD AUX TRANS" << ","
                        << "AB AUX BASAL" << ","
                        << "AB AUX DECAY" << ","
                        << "AB AUX DIFF" << ","
                        << "AB AUX TRANS" << ","
                        << "MD AUX BASAL" << ","
                        << "MD AUX DECAY" << ","
                        << "MD AUX DIFF" << ","
                        << "MD AUX TRANS" << ","
                        << "TOTAL AUXIN CYT" << ","
                        << "TOTAL AUXIN INTER" << ","
                        << "TOTAL AUXIN" << endl
                        << flush;
          }
      }
      return true;
}


  bool ChemicalSolver::step(bool timeToDebug, int stepCount,double KPin1auxinflux, double KPin1MF, double KPin1MFauxinflux) {
       if(!mesh)
           throw(QString("ChemicalsSolver::step No current mesh"));

       CCStructure& cs = mesh->ccStructure("Tissue");
       CCStructure& tt = mesh->ccStructure("Triangulate");
       CCIndexDataAttr& indexAttr = mesh->indexAttr();
       FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
       EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");

       double prodAD=0, prodAB=0, prodMD=0;
       double prodADL1=0, prodABL1=0, prodMDL1=0;
       double pinAD=0, pinAB=0, pinMD=0;
       double pinADL1=0, pinABL1=0, pinMDL1=0;
       double finalAD=0, finalAB=0, finalMD=0;
       double finalADL1=0, finalABL1=0,finalMDL1=0;
       double sensAD=0, sensAB=0, sensMD=0;
       double sensADL1=0, sensABL1=0, sensMDL1=0;
       double excessAD=0, excessAB=0, excessMD=0;
       double auxBasalAD=0, auxBasalAB=0, auxBasalMD=0;
       double auxDecayAD=0, auxDecayAB=0, auxDecayMD=0;
       double auxDiffAD=0, auxDiffAB=0, auxDiffMD=0;
       double auxTransAD=0, auxTransAB=0, auxTransMD=0;
       double topAuxin=0;
       double bottomAuxin=0;


       // Calculate chemicals changes on cells
       for(CCIndex f:cs.faces()){
               if(faceAttr[f].cellZone==AD || faceAttr[f].cellZone==AB ||faceAttr[f].cellZone==MD){
                  chemicsProcess->calcDerivsCell(cs, indexAttr, faceAttr, edgeAttr, f, Dt, prodAD, prodAB, prodADL1, prodABL1,
                                                 pinAD, pinAB,pinADL1, pinABL1,
                                                 finalAD,finalAB, finalADL1,finalABL1,
                                                 excessAD, excessAB,excessMD,
                                                 auxBasalAD,auxDecayAD,auxDiffAD, auxTransAD,auxBasalAB,auxDecayAB,auxDiffAB,auxTransAB,
                                                 stepCount,
                                                 topAuxin, bottomAuxin,
                                                 prodMD, prodMDL1, finalMD, finalMDL1, pinMD, pinMDL1, auxBasalMD, auxDecayMD, auxDiffMD, auxTransMD);
               }
       }

       //Debug
       double fprodAD =prodAD;
       double fprodAB=prodAB;
       double fprodMD=prodMD;
       double fprodADL1=prodADL1;
       double fprodABL1=prodABL1;
       double fprodMDL1=prodMDL1;
       double ffinalAD =finalAD;
       double ffinalAB =finalAB;
       double ffinalMD =finalMD;
       double ffinalADL1 =finalADL1;
       double ffinalABL1 =finalABL1;
       double ffinalMDL1 =finalMDL1;
       double fpinAD =pinAD;
       double fpinAB=pinAB;
       double fpinMD=pinMD;
       double fpinADL1 =pinADL1;
       double fpinABL1 =pinABL1;
       double fpinMDL1 =pinMDL1;
       double fexcessAD =excessAD;
       double fexcessAB =excessAB;
       double fexcessMD =excessMD;
       double fauxBasalAD=auxBasalAD;
       double fauxDecayAD=auxDecayAD;
       double fauxDiffAD=auxDiffAD;
       double fauxTransAD=auxTransAD;
       double fauxBasalAB=auxBasalAB;
       double fauxDecayAB=auxDecayAB;
       double fauxDiffAB=auxDiffAB;
       double fauxTransAB=auxTransAB;
       double fauxBasalMD=auxBasalMD;
       double fauxDecayMD=auxDecayMD;
       double fauxDiffMD=auxDiffMD;
       double fauxTransMD=auxTransMD;

       double normSum = 0;

       //Update chemicals
       for(CCIndex f:cs.faces()) { //neighbours cells
          FaceData& fD=faceAttr[f];

  //PIN1 SENSITIVITY
          //Take Strain direction for MF PIN1 polarization
          if(parm("Strain Measure Method") == "Given"){
            fD.strain=norm(fD.stiffnessXX)*(Rotate90(fD.stress)/(norm(fD.stress)+0.00001));
          }
          else if(parm("Strain Measure Method") == "Cell Axis"){
             fD.strain=fD.a1;
          }
       }

       //outwardNormal vector and aligment with stress/strain vector and auxin flux/gradient
       for (CCIndex e:cs.edges()){
           EdgeData& eD=edgeAttr[e];
           std::pair<CCIndex, CCIndex> eb=cs.edgeBounds(e);
           Point3d edge = indexAttr[eb.second].pos-indexAttr[eb.first].pos;
           std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);

           if(edgeFaces.size()==1) //No PIN1 in the borders
               continue;
           else if(edgeFaces.size()==2){
           for (CCIndex f:edgeFaces){   //incident cells
               FaceData& fD=faceAttr[f];
               int label =indexAttr[f].label;
               eD.edge_export_prev[label]=eD.edge_export[label]; //*1

               Point3d outwardNormal = Rotate90(edge);

               Point3d insideFace=fD.newcentroid-eD.midPoint; //Careful with orientation, from midPoint to centroid (head of the arrow)
               double angle =acos((insideFace*outwardNormal)/(norm(insideFace)*norm(outwardNormal))); //gives the smallest angle, <M_PI
               if (angle<(M_PI/2)) //change the orientation to be always out the cell
                   outwardNormal *=-1;

               // Calculate auxin flux impact
               eD.auxinFluxImpact[label] = 0;
               double auxinImpactValue=0;
               if(norm(fD.auxinFluxVector) > 0)
                   auxinImpactValue = ((fD.auxinFluxVector*outwardNormal)/(norm(fD.auxinFluxVector)*norm(outwardNormal)));
               if(auxinImpactValue < 0)
                   auxinImpactValue = 0;
               eD.auxinFluxImpact[label] =
                   norm(fD.auxinFluxVector)*(pow(auxinImpactValue, 4) / (pow(KauxinFlux, 4) + pow(auxinImpactValue, 4)));

               // Calculate MF impact (strain impact)
               double MFvalue = ((fD.strain*outwardNormal)/(norm(fD.strain)*norm(outwardNormal))); //direction
               if(norm(fD.strain)==0) //when a cell divides there is no strain
                   MFvalue=0;
               if(MFvalue < 0)
                   MFvalue *= -1;
               eD.MFImpact[label] =
                       norm(fD.strain)*(pow(MFvalue, 4) / (pow(Kmf, 4) + pow(MFvalue, 4))) ;

               // Calculate raw Pin sensitivity
               eD.pin1Sensitivity[label] =
                     (
                         (KPin1MF * eD.MFImpact[label])
                     +   (KPin1auxinflux * eD.auxinFluxImpact[label])
                     +   (KPin1MFauxinflux * eD.MFImpact[label] * eD.auxinFluxImpact[label])
                          );

               // Soft PIN1 sensitivity
               if(parm("Pin1 Sensitivity Average Method") == "Soft-max")
                   normSum += exp(eD.pin1Sensitivity[label]);
               else if(parm("Pin1 Sensitivity Average Method") == "Arithmetic Average")
                   normSum += eD.pin1Sensitivity[label];
             }
             }
       }
       for (CCIndex e:cs.edges()){ //iter over all although only required for edge between 2 faces in L1
           EdgeData& eD=edgeAttr[e];
           for(CCIndex f : cs.incidentCells(e, 2)){
              if(normSum > 0) {
                int label=indexAttr[f].label;
                if(parm("Pin1 Sensitivity Average Method") == "None")
                  break;
                else if(parm("Pin1 Sensitivity Average Method") == "Soft-max")
                  eD.pin1Sensitivity[label] = exp(eD.pin1Sensitivity[label]) / normSum;
                else if(parm("Pin1 Sensitivity Average Method") == "Arithmetic Average")
                  eD.pin1Sensitivity[label] = eD.pin1Sensitivity[label] / normSum;
              }
          }
      }

      //Debug(average sensivity)
       for (CCIndex e:cs.edges()){
           EdgeData& eD=edgeAttr[e];
           std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);
           if(edgeFaces.size()==1) //No PIN1 in the borders
               continue;
           else if(edgeFaces.size()==2){
           for (CCIndex f:edgeFaces){   //incident cells
               FaceData& fD=faceAttr[f];
               int label=indexAttr[f].label;
               if(fD.cellZone==AD){
                 sensAD +=eD.pin1Sensitivity[label];
                 if(fD.cellLayer==L1)
                     sensADL1 +=eD.pin1Sensitivity[label];
               }
               else if(fD.cellZone==AB){
                 sensAB +=eD.pin1Sensitivity[label];
                 if(fD.cellLayer==L1)
                     sensABL1 +=eD.pin1Sensitivity[label];
               }
               else if(fD.cellZone==MD){
                 sensMD +=eD.pin1Sensitivity[label];
                 if(fD.cellLayer==L1)
                     sensMDL1 +=eD.pin1Sensitivity[label];
               }
           }
           }
       }


  //UPDATE MAPS AND DEBUG
       //Cytoplasmic auxin
       if(timeToDebug){
       CCIndexDoubleAttr& auxinCyt = mesh->signalAttr<double>("Auxin Cyt");
       auxinCyt.clear();
       double totalCyt=0;
       for(CCIndex f:cs.faces()){
           auxinCyt[f] = (faceAttr[f].auxin/indexAttr[f].measure);
           totalCyt+=faceAttr[f].auxin;
       }

       //Intercelular auxin and PIN1 in walls
       double totalOut=0;
       for(CCIndex e:cs.edges()){
           EdgeData& eD=edgeAttr[e];
           totalOut +=eD.intercellularAuxin;
       }

       CCIndexDoubleAttr& auxinOut = mesh->signalAttr<double>("Auxin Out");
       CCIndexDoubleAttr& pin1 = mesh->signalAttr<double>("PIN1");
       auxinOut.clear();
       pin1.clear();
       for(CCIndex f:tt.faces()){
           double intercellularAuxin=0;
           double edgePin1=0;
           for(CCIndex e:tt.incidentCells(f,1)){
              CCIndex n=edgeAttr[e].sameEdge;
              intercellularAuxin+=edgeAttr[n].intercellularAuxin;
              edgePin1+=edgeAttr[n].pin1[indexAttr[f].label];
           }
           auxinOut[f]=intercellularAuxin;
            pin1[f] = edgePin1;
       }

       mesh->updateAll();

  //PRINT INFORMATION
       int nAD=0; int nAB=0;
       int nADL1=0; int nABL1=0;
       int nADedge=0; int nABedge=0;
       int nADL1edge=0; int nABL1edge=0;
       int nMD=0; int nMDL1=0; int nMDedge=0; int nMDL1edge=0;
       int nADinner=0, nABinner=0, nMDinner=0;
       double L1AD=0, innerAD=0, L1AB=0, innerAB=0, L1MD=0, innerMD=0;

       for(CCIndex f:cs.faces()){
           if(faceAttr[f].cellZone==AD){
               nAD ++;
               if(faceAttr[f].cellLayer==L1){
                   nADL1 ++;
                   L1AD += faceAttr[f].auxin/indexAttr[f].measure;
               } else{
                 nADinner ++;
                 innerAD += faceAttr[f].auxin/indexAttr[f].measure;
               }
           std::set<CCIndex> faceEdges =cs.incidentCells(f,1);
           for (CCIndex e:faceEdges){ //for all the edges in the boundary of the cell
                std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);
                if(edgeFaces.size()==2){
                    nADedge ++;
                    if(faceAttr[f].cellLayer==L1)
                        nADL1edge ++;
                }
           }
           }
           else if (faceAttr[f].cellZone==AB){
               nAB ++;
               if(faceAttr[f].cellLayer==L1){
                   nABL1 ++;
                   L1AB += faceAttr[f].auxin/indexAttr[f].measure;
               } else{
                   nABinner ++;
                   innerAB += faceAttr[f].auxin/indexAttr[f].measure;
               }
               std::set<CCIndex> faceEdges =cs.incidentCells(f,1);
               for (CCIndex e:faceEdges){ //for all the edges in the boundary of the cell
                    std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);
                    if(edgeFaces.size()==2){
                        nABedge ++;
                        if(faceAttr[f].cellLayer==L1)
                            nABL1edge ++;
                    }
               }
           }
           else if (faceAttr[f].cellZone==MD){
               nMD ++;
               if(faceAttr[f].cellLayer==L1){
                   nMDL1 ++;
                   L1MD += faceAttr[f].auxin/indexAttr[f].measure;
               } else{
                   nMDinner ++;
                   innerMD += faceAttr[f].auxin/indexAttr[f].measure;
               }
               std::set<CCIndex> faceEdges =cs.incidentCells(f,1);
               for (CCIndex e:faceEdges){ //for all the edges in the boundary of the cell
                    std::set<CCIndex> edgeFaces=cs.incidentCells(e,2);
                    if(edgeFaces.size()==2){
                        nMDedge ++;
                        if(faceAttr[f].cellLayer==L1)
                            nMDL1edge ++;
                    }
               }
           }
       }

       if(parm("Debug") == "True"){
       output_file << stepCount << ","
                   << topAuxin <<","
                   << bottomAuxin <<","
                   << innerAD/nADinner <<","
                   << innerAB/nABinner <<","
                   << innerMD/nMDinner <<","
                   << L1AD/nADL1 <<","
                   << L1AB/nABL1 <<","
                   << L1MD/nMDL1 <<","
                   << fprodAD/nAD  << ","
                   << fprodAB/nAB  << ","
                   << fprodMD/nMD  << ","
                   << fprodADL1/nADL1  << ","
                   << fprodABL1/nABL1  << ","
                   << fprodMDL1/nMDL1  << ","
                   << fpinAD/nADedge << ","
                   << fpinAB/nABedge << ","
                   << fpinMD/nMDedge << ","
                   << fpinADL1/nADL1edge << ","
                   << fpinABL1/nABL1edge << ","
                   << fpinMDL1/nMDL1edge << ","
                   << (ffinalAD-fprodAD)/nAD  << ","
                   << (ffinalAB-fprodAB)/nAB  << ","
                   << (ffinalMD-fprodMD)/nMD  << ","
                   << (ffinalADL1-fprodADL1)/nADL1  << ","
                   << (ffinalABL1-fprodABL1)/nABL1  << ","
                   << (ffinalMDL1-fprodMDL1)/nMDL1  << ","
                   << sensAD/nADedge  << ","
                   << sensAB/nABedge  << ","
                   << sensMD/nMDedge  << ","
                   << sensADL1/nADL1edge  << ","
                   << sensABL1/nABL1edge  << ","
                   << sensMDL1/nMDL1edge  << ","
                   << fexcessAD/nAD << ","
                   << fexcessAB/nAB  << ","
                   << fexcessMD/nMD  << ","
                   << fauxBasalAD/nAD << ","
                   << fauxDecayAD/nAD << ","
                   << fauxDiffAD/nAD << ","
                   << fauxTransAD/nAD << ","
                   << fauxBasalAB/nAB << ","
                   << fauxDecayAB/nAB << ","
                   << fauxDiffAB/nAB << ","
                   << fauxTransAB/nAB << ","
                   << fauxBasalMD/nMD << ","
                   << fauxDecayMD/nMD << ","
                   << fauxDiffMD/nMD << ","
                   << fauxTransMD/nMD << ","
                   <<totalCyt <<","
                   <<totalOut<< ","
                   << totalCyt+totalOut<< endl
                   << flush;
       }

    }

      return true;
  }

  ////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////MECHANICAL SOLVER////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////

  bool Mechanics::initialize(QWidget* parent) {
      mesh = currentMesh();
      if(!mesh)
        throw(QString("General::initialize No current mesh"));
      CCStructure& cs = mesh->ccStructure("Triangulate"); //TISSUE TO GROWTH
      CCIndexDataAttr& indexAttr = mesh->indexAttr();
      EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");

      Pressure=parm("Pressure").toDouble();

      //border (of the mesh) and wall(border of the cell) edges (update when triangulate)//if triangulate always first here is not required
      for(CCIndex e:cs.edges()) {
         EdgeData& eD = edgeAttr[e];
         CCIndexSet cb = cs.cobounds(e);
         std::vector<CCIndex> cbV(cb.begin(), cb.end());
           if(cb.size() == 1) { // edge in the border
             eD.border = true;
             eD.f = *cb.begin();
             eD.g = (CCIndex)0;
             eD.wall=true;
           } else { //edge not in the border
               eD.border = false;
               if (indexAttr[cbV[0]].label !=indexAttr[cbV[1]].label)
                   eD.wall = true;
               else
                   eD.wall=false;

           }
      }

        return true;
  }


  void Mechanics::calcForces(CCIndex v, const CCStructure& cs, CCIndexDataAttr& indexAttr,
                             FaceDataAttr& faceAttr, EdgeDataAttr& edgeAttr, VerticeDataAttr& verticeAttr) {

      VerticeData& vD = verticeAttr[v];
      vD.force=0;

      for(const Flip& flip : cs.matchV(CCIndex::BOTTOM, v, CCIndex::Q, CCIndex::Q)) { //Move for each vertex
         CCIndexData& vIdx = indexAttr[v];
         CCIndex e = flip.interior; //edge btw
         CCIndex n = flip.otherFacet(v); //other vertex
         CCIndexData& nIdx = indexAttr[n];
         EdgeData& eD = edgeAttr[e];
         VerticeData& vD = verticeAttr[v];
         if(eD.border) {
         //calculate the length of the edge
         Point3d dir =nIdx.pos - vIdx.pos;
         double length = norm(dir);
         if(length == 0) {
           continue;
         }
        dir /= length;
        // Pressure on the outside edge
        // Find orientation, translate to + and -
        int ro = cs.ro(eD.f, e) * cs.ro(e, v);
        ro = ro == 0 ? 1 : -1;
        vD.force += ro * Point3d(dir.y(), -dir.x(), 0) * Pressure * length * 0.5;
                  }
      }
  }

  ////----------------------------------PBD CONSTRAINTS-----------------------------//

  bool  solve_DihedralConstraint(
       const Point3d &p0, double invMass0,
       const Point3d &p1, double invMass1,
       const Point3d &p2, double invMass2,
       const Point3d &p3, double invMass3,
       const double restAngle,
       const double stiffness,
       Point3d &corr0, Point3d &corr1, Point3d &corr2, Point3d &corr3, bool verbose=false)
   {
       // derivatives from Bridson, Simulation of Clothing with Folds and Wrinkles
       // his modes correspond to the derivatives of the bending angle arccos(n1 dot n2) with correct scaling

       if (invMass0 == 0.0 && invMass1 == 0.0)
           return false;

       Point3d e = p3-p2;
       double  elen = e.norm();
       if (elen < EPS)
           return false;

       double invElen = 1.0 / elen;

       Point3d n1 = (p2-p0).cross(p3-p0); n1 /= n1.norm() * n1.norm();
       Point3d n2 = (p3 - p1).cross(p2 - p1); n2 /= n2.norm() * n2.norm();

       Point3d d0 = elen*n1;
       Point3d d1 = elen*n2;
       Point3d d2 = ((p0-p3) * e) * invElen * n1 + ((p1-p3) * e) * invElen * n2;
       Point3d d3 = ((p2-p0) * e) * invElen * n1 + ((p2-p1) * e) * invElen * n2;

       n1 /= n1.norm();
       n2 /= n2.norm();
       double dot = n1 * n2;

       if (dot < -1.0) dot = -1.0;
       if (dot >  1.0) dot =  1.0;
       double phi = acos(dot);

       if(verbose)
           mdxInfo << "Dihedral angle: " << phi* (180./M_PI) << " restAngle: " << restAngle* (180./M_PI)  <<  endl;

       double lambda =
           invMass0 * d0.norm() * d0.norm() +
           invMass1 * d1.norm() * d1.norm() +
           invMass2 * d2.norm() * d2.norm() +
           invMass3 * d3.norm() * d3.norm();

       if (lambda == 0.0)
           return false;

       lambda = (phi - restAngle) / lambda * stiffness;

       if (n1.cross(n2) * e > 0.0) ///////////// I changed < to >
           lambda = -lambda;

      corr0 = - invMass0 * lambda * d0;
      corr1 = - invMass1 * lambda * d1;
      corr2 = - invMass2 * lambda * d2;
      corr3 = - invMass3 * lambda * d3;

      return true;
  }

  ////-----------------------------------------------------------------------//

  bool init_StrainTriangleConstraint( //just calculate inverse matrix of the rest situation
      const Point3d &p0,
      const Point3d &p1,
      const Point3d &p2,
      Matrix2d &invRestMat)
  {
      double a = p1[0] - p0[0]; double b = p2[0] - p0[0];
      double c = p1[1] - p0[1]; double d = p2[1] - p0[1];

      // inverse
      double det = a*d - b*c; //si no es capaz de initializar es por esto
      if (fabs(det) < EPS)
          return false;

      double s = 1. / det;
      invRestMat[0][0] =  d*s;  invRestMat[0][1] = -b*s;
      invRestMat[1][0] = -c*s;  invRestMat[1][1] =  a*s;

      return true;
  }


  bool solve_StrainTriangleConstraint(
           const Point3d &p0, double invMass0,
           const Point3d &p1, double invMass1,
           const Point3d &p2, double invMass2,
           const Matrix2d &invRestMat,
           const double xxStiffness, //in the x direction; 0
           const double yyStiffness, //in the y direction; 0
           const double xyStiffness,
           const bool normalizeStretch, //normalize values before (true as default)
           const bool normalizeShear,
           Point3d &corr0, Point3d &corr1, Point3d &corr2)
   {
       Point3d c[2];
       c[0] = Point3d(invRestMat[0][0], invRestMat[1][0], 0.0);
       c[1] = Point3d(invRestMat[0][1], invRestMat[1][1], 0.0);

       Point3d r[3];

       for (int i = 0; i < 2; i++) {
           for (int j = 0; j <= i; j++) {


               r[0] = Point3d((p1[0] + corr1[0]) - (p0[0] + corr0[0]), (p2[0] + corr2[0]) - (p0[0] + corr0[0]), 0.0);     // Gauss - Seidel
               r[1] = Point3d((p1[1] + corr1[1]) - (p0[1] + corr0[1]), (p2[1] + corr2[1]) - (p0[1] + corr0[1]), 0.0);
               r[2] = Point3d((p1[2] + corr1[2]) - (p0[2] + corr0[2]), (p2[2] + corr2[2]) - (p0[2] + corr0[2]), 0.0);


               double Sij = 0.0;
               for (int k = 0; k < 3; k++)
                   Sij += (r[k] * (c[i])) * (r[k] * (c[j]));

               Point3d d[3];
               d[0] = Point3d(0.0, 0.0, 0.0);

               for (int k = 0; k < 2; k++) {
                   d[k+1]  = Point3d((r[0] * (c[j])), (r[1] * (c[j])), (r[2] * (c[j]))) * invRestMat(k, i);
                   d[k+1] += Point3d((r[0] * (c[i])), (r[1] * (c[i])), (r[2] * (c[i]))) * invRestMat(k, j);
                   d[0] -= d[k+1];
               }

               if (i != j && normalizeShear) {
                   double fi2 = 0.0;
                   double fj2 = 0.0;
                   for (int k = 0; k < 3; k++) {
                       fi2 += (r[k] * (c[i])) * (r[k] * (c[i]));
                       fj2 += (r[k] * (c[j])) * (r[k] * (c[j]));
                   }
                   double fi = sqrt(fi2);
                   double fj = sqrt(fj2);

                   d[0] = Point3d(0.0, 0.0, 0.0);
                   double s = Sij / (fi2*fi*fj2*fj);
                   for (int k = 0; k < 2; k++) {
                       d[k+1] /= fi * fj;
                       d[k+1] -= fj*fj * Point3d((r[0] * (c[i])), (r[1] * (c[i])), (r[2] * (c[i]))) * invRestMat[k][i] * s;
                       d[k+1] -= fi*fi * Point3d((r[0] * (c[j])), (r[1] * (c[j])), (r[2] * (c[j]))) * invRestMat[k][j] * s;
                       d[0] -= d[k+1];
                   }
                   Sij = Sij / (fi * fj);
               }

               double lambda =
                   invMass0 * d[0].norm() * d[0].norm()  +
                   invMass1 * d[1].norm() * d[1].norm()  +
                   invMass2 * d[2].norm() * d[2].norm();

               if (lambda == 0.0)
                   continue;

               if (i == 0 && j == 0) {
                   if (normalizeStretch) {
                       double s = sqrt(Sij);
                       lambda = static_cast<double>(2.0) * s * (s - static_cast<double>(1.0)) / lambda * xxStiffness;
                   }
                   else {
                       lambda = (Sij - static_cast<double>(1.0)) / lambda * xxStiffness;
                   }
               }
               else if (i == 1 && j == 1) {
                   if (normalizeStretch) {
                       double s = sqrt(Sij);
                       lambda = static_cast<double>(2.0) * s * (s - static_cast<double>(1.0)) / lambda * yyStiffness;
                   }
                   else {
                       lambda = (Sij - static_cast<double>(1.0)) / lambda * yyStiffness;
                   }
               }
               else {
                   lambda = Sij / lambda * xyStiffness;
               }

               corr0 -= lambda * invMass0 * d[0];
               corr1 -= lambda * invMass1 * d[1];
               corr2 -= lambda * invMass2 * d[2];
           }
       }
       return true;
  }

  ////-----------------------------------------------------------------------//
  bool init_ShapeMatchingConstraint(const std::vector<Point3d> x0, //rest configuration
                                    const std::vector<double> invMasses,
                                    int numPoints, //number of vertices
                                    Point3d& restCm, //return center of mass of the rest configuration = centroid no mass?¿
                                    Matrix3d& invRestMat) { //return a matrix required by the solver
      invRestMat.identity();

      // center of mass
      restCm = 0;
      double wsum = 0.0;
      for(int i = 0; i < numPoints; i++) {
          double wi = 1. / (invMasses[i] + EPS);
          restCm += x0[i] * wi;
          wsum += wi;
      }
      if(wsum == 0.0)
          return false;
      restCm /= wsum;

      // A
      Matrix3d A;
      A = 0;
      for(int i = 0; i < numPoints; i++) {
          const Point3d qi = x0[i] - restCm;
          double wi = 1.0 / (invMasses[i] + EPS);
          double x2 = wi * qi[0] * qi[0]; //each dimension
          double y2 = wi * qi[1] * qi[1];
          //double z2 = wi * qi[2] * qi[2];
          double xy = wi * qi[0] * qi[1];
          double xz = wi * qi[0] * qi[2];
          double yz = wi * qi[1] * qi[2];
          A[0] += Point3d(x2, xy, xz);
          A[1] += Point3d(xy, y2, yz);
          //A[2] += Point3d(xz, yz, z2);
          A[2] = Point3d(0, 0, 1); //to fill
          ;
      }
      double determinant = det(A);
      if(fabs(determinant) > EPS) {
          invRestMat = inverse(A);
          return true;
      }
      return false;
  }

  bool solve_ShapeMatchingConstraint(const std::vector<Point3d> x0, //rest configuration
                                     const std::vector<Point3d> x,  //current configuration
                                     const std::vector<double> invMasses,
                                     int numPoints,  //number of vertices
                                     const Point3d& restCm,
                                     const Matrix3d& invRestMat,
                                     bool allowStretch,  //allow stretching--TRUE
                                     std::vector<Point3d>& corr, //corrections
                                     Matrix3d* rot=NULL) {   //rotation matrix
      // center of mass
      Point3d cm(0.0, 0.0, 0.0);
      double wsum = 0.0;
      for(int i = 0; i < numPoints; i++) {
          double wi = 1.0 / (invMasses[i] + EPS);
          cm += x[i] * wi;
          wsum += wi;
      }
      if(wsum == 0.0)
          return false;
      cm /= wsum;

      // A
      Matrix3d mat;
      for(int i = 0; i < numPoints; i++) {
          Point3d q = x0[i] - restCm;
          Point3d p = x[i] - cm;
          double w = 1.0 / (invMasses[i] + EPS);
          p *= w;
          mat[0] += Point3d(p[0] * q[0], p[0] * q[1], p[0] * q[2]);
          mat[1] += Point3d(p[1] * q[0], p[1] * q[1], p[1] * q[2]);
          mat[2] += Point3d(p[2] * q[0], p[2] * q[1], p[2] * q[2]);
      }

      mat = mat * invRestMat;
      Matrix3d R, U;
      if (allowStretch)
        R = mat;
     // else
     //   PolarDecompX(mat, R, U);

      for(int i = 0; i < numPoints; i++) {
          Point3d goal = cm + R * (x0[i] - restCm);
          corr.push_back((goal - x[i]));
      }
      if(rot)
          *rot = R;

      return true;
  }

  ////-----------------------------------------------------------------------//
   bool solve_DistanceConstraint(
       const Point3d &p0, double invMass0,
       const Point3d &p1, double invMass1,
       const double restLength,
       const double compressionStiffness,
       const double extensionStiffness,
       Point3d &corr0,
       Point3d &corr1)
   {
       double wSum = invMass0 + invMass1;
       if (wSum == 0.0)
           return false;
       Point3d n = p1 - p0;
       double dis = n.norm(); //current length
       n /= dis; //direction

       Point3d corr;
       if (dis < restLength)
           corr = ((dis - restLength) * n  * compressionStiffness) / wSum; //high to not allow compression
       else
           corr = ((dis - restLength) * n  * extensionStiffness) / wSum; //low to allow extension

       corr0 = invMass0 * corr;
       corr1 = -invMass1 * corr;
       return true;
  }


   ////----------------------------------MECHANICAL SOLVER-----------------------------//

   bool MechanicalSolver::initialize(QWidget* parent, std::string debFileName) {

     mesh = getMesh("Mesh 1");
     if(!mesh)
        throw(QString("MechanicalSolver::initialize No current mesh"));

     if(!getProcess(parm("Mechanics Process"), mechanicsProcess))
        throw(QString("MechanicalSolver::initialize Cannot make Mechanics Process:") +
              parm("Mechanics Process"));

     CCIndexDataAttr& indexAttr = mesh->indexAttr();
     VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");
     EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");
     FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
     CCStructure& cs = mesh->ccStructure("Triangulate");
     CCStructure& ts = mesh->ccStructure("Tissue");

   //PARAMETERS
     //General
     Dt = parm("Dt").toDouble();
     staticDamping = parm("Velocity Static Damping").toDouble();
     velocityUpdate = parm("Velocity Update");
     maxVelocity = parm("Max Velocity").toDouble();
     PBDiterations = parm("PBD iterations").toInt();

     //Stiffness
     PBDspringsStiffness = parm("PBD springs stiffness").toDouble();
     springActive = stringToBool(parm("PBD spring active"));
     WallExtStiff = parm("Wall eStifness").toDouble();
     InExtStiff = parm("Interior eStifness").toDouble();
     ComStiffness = parm("cStifness").toDouble();

     PBDshapeStiffness = parm("PBD shape stiffness").toDouble();
     PBDshapeAllowStretching = parm("PBD shape allow stretching") == "True";

     PBDbendingStiffness = parm("PBD bending stiffness").toDouble();

     PBDstrainStiffness = parm("PBD strain stiffness").toDouble();
     PBDstrainXYStiffness = parm("PBD strain XY stiffness").toDouble();
     PBDstrainXXStiffness = parm("PBD strain XX stiffness").toDouble();
     PBDstrainYYStiffness = parm("PBD strain YY stiffness").toDouble();
     PBDshearNormalize = stringToBool(parm("PBD shear normalize"));
     PBDstretchNormalize = stringToBool(parm("PBD stretch normalize"));

     //Polarity
     for (CCIndex f:ts.faces()){
        if(norm(faceAttr[f].a1) == 0) //avoid errors in divisions
            faceAttr[f].a1[1] = EPS;
        if(norm(faceAttr[f].a2) == 0)
            faceAttr[f].a2[0] = EPS;
     }
     MFspringTresh = parm("MF Strain Threshold for Growth").toDouble();
     MFRORate = parm("MF reorientation rate").toDouble();
     MFdegradation = parm("MF Degradation").toDouble();

   //Velocity: initialize
     for (CCIndex v:cs.vertices())
        verticeAttr[v].velocity=Point3d(0,0,0);

   //PBDShape: initialize restCm
     for(CCIndex f:ts.faces()) {
        FaceData& fD = faceAttr[f];
        fD.restX0.clear();
        fD.cellVertices= faceVertices(ts, f);
        std::vector<double> invMasses;
        for(CCIndex v : fD.cellVertices) {
            fD.restX0.push_back(Point3d(indexAttr[v].pos));
            invMasses.push_back(verticeAttr[v].invMass);
        }
        init_ShapeMatchingConstraint(fD.restX0, invMasses, fD.restX0.size(),  fD.restCm, fD.invRestMat);
     }

   //PBDSprings: initialize restLength and stiffness
     for(CCIndex e:cs.edges()) {
        EdgeData& eD = edgeAttr[e];
        auto eb = cs.edgeBounds(e);
        double length = norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos);
        eD.restLength = length;
     }

   //CALL MECHANICS
       mechanicsProcess->initialize(parent);

       return true;
   }


   ////----------------------------------MOMENTUM-----------------------------//
   void MechanicalSolver::momentumPreservation(const CCStructure& cs, CCIndexDataAttr& indexAttr, VerticeDataAttr& verticeAttr, Point3d P_r, Point3d L_r, double& meanCorrection) {

        double totalMass=0;
        Point3d centerOfMasses;
        Point3d P_pbd; //linear momentum
        Point3d L_pbd; //angular momentum
        Matrix3d I; //inertia tensor
        for (CCIndex v:cs.vertices()){
            VerticeData& vD=verticeAttr[v];
            CCIndexData& iD=indexAttr[v];
            double vertexMass=1/vD.invMass;
            totalMass += vertexMass;
            centerOfMasses += iD.pos * vertexMass;
        }
        centerOfMasses /=totalMass;

        for(CCIndex v:cs.vertices()){
            VerticeData& vD=verticeAttr[v];
            CCIndexData& iD=indexAttr[v];
            double vertexMass=1/vD.invMass;
            Point3d r=iD.pos-centerOfMasses;
            Matrix3d r_star;
            r_star[0] = Point3d(0, -r[2], r[1]);
            r_star[1] = Point3d(r[2], 0, -r[0]);
            r_star[2] = Point3d(-r[1], r[0], 0);
            I +=vertexMass*r_star*transpose(r_star);
            P_pbd +=vertexMass*vD.velocity;
            L_pbd +=vertexMass*r.cross(vD.velocity);
        }
        Point3d v_cor=(P_r-P_pbd)/totalMass;
        Point3d omega_cor = inverse(I) * (L_r-L_pbd);

        for(CCIndex v:cs.vertices()){
            VerticeData& vD=verticeAttr[v];
            CCIndexData& iD=indexAttr[v];
            vD.velocity+=v_cor+omega_cor.cross(iD.pos - centerOfMasses);
            meanCorrection+= norm(v_cor+omega_cor.cross(iD.pos - centerOfMasses));
        }
        meanCorrection /=cs.vertices().size();
   }


   ////----------------------------------VELOCITY-----------------------------//

   void MechanicalSolver::solveEuler1stOrder(CCIndex v, const CCStructure& cs, CCIndexDataAttr& indexAttr, VerticeDataAttr& verticeAttr,
                                               Point3d& P_r, Point3d& L_r, Point3d centerOfMasses) {

       VerticeData& vD = verticeAttr[v];
       CCIndexData& iD = indexAttr[v];
       vD.prevVelocity = vD.velocity;
       vD.velocity += vD.force * vD.invMass *Dt;
       vD.velocity -= vD.velocity * staticDamping;

       vD.lastPos = vD.prevPos;
       vD.prevPos = iD.pos;
       iD.pos += vD.velocity * Dt;

       //momentumPreservation
       Point3d p=(1/vD.invMass)*(vD.velocity - vD.prevVelocity);//linear moemntum increase
       Point3d r=(iD.pos - centerOfMasses);
       P_r +=p;
       L_r +=r.cross(p);
   }

  //PBD SOLVE

   void MechanicalSolver::solveConstraints(double& dStrainADIN, double& dStrainADL1, double& dStrainABIN, double& dStrainABL1, double& dStrainMDIN, double& dStrainMDL1, int& stepCount) {

      CCIndexDataAttr& indexAttr = mesh->indexAttr();
      VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");
      EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");
      FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
      CCStructure& cs = mesh->ccStructure("Triangulate");
      CCStructure& ts = mesh->ccStructure("Tissue");

      time +=Dt;
      //Inner stiffness
      for(CCIndex e:cs.edges()){
         edgeAttr[e].stiffness=InExtStiff;
      }

      for(CCIndex f: ts.faces()){
        double stiffness=0;
        FaceData& fD = faceAttr[f];

        if(faceAttr[f].cellPosition==max)
              stiffness=1;
        else if(faceAttr[f].cellPosition==nine)
           stiffness=0.9;
        else if(faceAttr[f].cellPosition==eight)
           stiffness=0.8;
        else if(faceAttr[f].cellPosition==seven)
           stiffness=0.7;
        else if(faceAttr[f].cellPosition==six)
           stiffness=0.6;
        else if(faceAttr[f].cellPosition==five)
           stiffness=0.5;
        else if(faceAttr[f].cellPosition==four)
           stiffness=0.4;
        else if(faceAttr[f].cellPosition==three)
           stiffness=0.27;
        else if(faceAttr[f].cellPosition==two)
           stiffness=0.2;
        else if(faceAttr[f].cellPosition==one)
           stiffness=0.1;
        else if(faceAttr[f].cellPosition==pointfive)
           stiffness=0.01;
        else if(faceAttr[f].cellPosition==pointone)
           stiffness=0.01;
        else
           throw(QString("Wrong cellPosition defined"));

        if(faceAttr[f].cellLayer==L1){
            stiffness=stiffness+(stiffness*(stepCount/FinalStep)); //Increase with time
        }

        if(faceAttr[f].cellZone==AB || faceAttr[f].cellZone==AD)
           stiffness/=2.4; //scale

        for(CCIndex e:ts.incidentCells(f,1)){ //always get the highest stiffness
            CCIndex es=edgeAttr[e].sameEdge;
                if(stiffness>edgeAttr[es].stiffness)
                   edgeAttr[es].stiffness=stiffness;
                if(edgeAttr[es].stiffness>1) //max limit always 1
                    edgeAttr[es].stiffness=1;
                else if(edgeAttr[es].stiffness<InExtStiff)
                    edgeAttr[es].stiffness=InExtStiff; //min limit always the interior
        }
       }

      //Strain: Initialize Strain direction
      Point3d ref_pos=Point3d(0,-1000,0); //?¿?¿

      for(CCIndex v:cs.vertices()){
         if(verticeAttr[v].stiffRefVertex){
            ref_pos=indexAttr[v].pos;
            break;
         }
      }

      //Pass attr to traingulate faces
        for(CCIndex f:ts.faces()){
          for(CCIndex t:faceAttr[f].cellFaces){
              faceAttr[t].a1=faceAttr[f].a1;
              faceAttr[t].a2=faceAttr[f].a2;
              faceAttr[t].force=faceAttr[f].force;
              faceAttr[t].stress=faceAttr[f].stress;
              faceAttr[t].auxinFluxVector=faceAttr[f].auxinFluxVector;
              faceAttr[t].newcentroid=faceAttr[f].newcentroid;
          }
        }

     //Pass strain constraint direction and initialize
       for(CCIndex f: cs.faces()) {
           FaceData& fD = faceAttr[f];
           std::vector<CCIndex> fv = faceVertices(cs, f);
           if (fD.cellLayer==L1)
              fD.stress=fD.force;
           else
              fD.stress=Rotate90(fD.a1);

           double faceAngle = mdx::angle(Point3d(1,0,0), fD.stress);//xx direction of maximal stiffness

           Point3d x0 = RotateCounterClock(verticeAttr[fv[0]].prevPos, faceAngle);
           Point3d x1 = RotateCounterClock(verticeAttr[fv[1]].prevPos, faceAngle);
           Point3d x2 = RotateCounterClock(verticeAttr[fv[2]].prevPos, faceAngle);
           if(!init_StrainTriangleConstraint(x0, x1, x2, fD.strainInvRestMat)) {
               mdxInfo << "solveConstraints: Failed initialization of strain constraint for face: " << f << endl;
               continue;
           }

       }


      //APPLY CONSTRAINTS
            //Distance
            for(int inter = 0; inter < PBDiterations; inter++) {
            if(PBDspringsStiffness > 0) {
                for (CCIndex e:cs.edges()){
                    EdgeData& eD = edgeAttr[e];
                    auto eb = cs.edgeBounds(e);
                    Point3d p0 = indexAttr[eb.first].pos;
                    Point3d p1 = indexAttr[eb.second].pos;
                    double m0 = verticeAttr[eb.first].invMass;
                    double m1 = verticeAttr[eb.second].invMass;

                    double eStiffness;
                    if (springActive)
                       eStiffness=eD.stiffness;
                    else{
                        if(eD.wall==true)
                        eStiffness=WallExtStiff;
                        else
                        eStiffness=InExtStiff;
                    }
                    Point3d corr0, corr1;
                    solve_DistanceConstraint(p0, m0, p1,m1, eD.restLength, ComStiffness, eStiffness, corr0, corr1);
                    double stiffness = PBDspringsStiffness;
                    corr0 *= stiffness;
                    corr1 *= stiffness;
                    verticeAttr[eb.first].corrections["distance"] += corr0;
                    verticeAttr[eb.second].corrections["distance"] += corr1;
                    indexAttr[eb.first].pos += corr0;
                    indexAttr[eb.second].pos += corr1;
                }
            }

            //Shape
            if(PBDshapeStiffness> 0) {
                ////#pragma omp parallel for schedule(static)
                for(CCIndex f: ts.faces()) { //For tissue vertices
                    FaceData& fD = faceAttr[f];
                    std::vector<Point3d> X;
                    std::vector<double> invMasses;
                    for(CCIndex v : fD.cellVertices) { //take all vertices of a cell
                        X.push_back(indexAttr[v].pos);
                        invMasses.push_back(verticeAttr[v].invMass);
                    }
                    std::vector<Point3d> corr;
                    Matrix3d rot;
                    solve_ShapeMatchingConstraint(fD.restX0, X, invMasses, fD.cellVertices.size(), fD.restCm,
                                                  fD.invRestMat, PBDshapeAllowStretching, corr, &rot);
                    double stiffness;
                    if ((fD.cellZone==PZ || fD.cellZone==CZ) && fD.cellLayer==L1){
                      stiffness=0.8;
                    }
                    else
                      stiffness = PBDshapeStiffness;
                    int k=0;
                    for(CCIndex v : fD.cellVertices) {
                        Point3d corrv = corr[k] * stiffness;
                        verticeAttr[v].corrections["shape"] += corrv;
                        indexAttr[v].pos += corrv;
                        k++;
                    }
                }
            }

            //Bending
            if(PBDbendingStiffness > 0) { //only for tissue vertices
                for(CCIndex f: ts.faces()) {
                    FaceData& fD = faceAttr[f];
                    int n=fD.cellVertices.size(); //updated in shape initialize
                    double stiffness;
                    stiffness = PBDbendingStiffness;
                for(int i = 0; i < n; i++) {
                        CCIndex v = fD.cellVertices[i];
                        CCIndex prev = (i>0) ? fD.cellVertices[i-1] : fD.cellVertices[n-1];
                        CCIndex next = (i<n-1) ? fD.cellVertices[i+1] : fD.cellVertices[0];
                        CCIndex e0, e1;
                        for(const Flip& flip : cs.matchV(v, CCIndex::Q, CCIndex::Q, f)) {
                            e0=flip.facet[0];
                            e1=flip.facet[1];
                        }
                        Point3d a0 = indexAttr[prev].pos - indexAttr[v].pos;
                        Point3d a1 = indexAttr[next].pos - indexAttr[v].pos;
                        double angle = mdx::angle(a0, a1);
                        double restAngle = fmod(angle-M_PI, M_PI);
                        if(restAngle < 0)
                            restAngle *= -1;
                        Point3d p0 = indexAttr[v].pos;
                        Point3d p1 = indexAttr[v].pos + Point3d(0,0,-1);
                        Point3d p2 = indexAttr[prev].pos;
                        Point3d p3 = indexAttr[next].pos;
                        double m0 = verticeAttr[v].invMass;
                        double m1 = verticeAttr[v].invMass;
                        double m2 = verticeAttr[prev].invMass;
                        double m3 = verticeAttr[next].invMass;
                        Point3d corr0, corr1, corr2, corr3;

                        solve_DihedralConstraint(p2, m2, p3, m3, p0, m0, p1, m1, restAngle, stiffness, corr2, corr3, corr0, corr1, indexAttr[v].selected);
                        verticeAttr[v].corrections["bending"] += corr0;
                        verticeAttr[prev].corrections["bending"] += corr2;
                        verticeAttr[next].corrections["bending"] += corr3;
                        indexAttr[v].pos += corr0;
                        indexAttr[prev].pos += corr2;
                        indexAttr[next].pos += corr3;
                }
                }
            }

            //Strain
            if(PBDstrainStiffness > 0) {
                for(CCIndex f: cs.faces()) {
                    FaceData& fD = faceAttr[f];
                    std::vector<CCIndex> fv = faceVertices(cs, f);
                    Point3d p0 = indexAttr[fv[0]].pos;
                    Point3d p1 = indexAttr[fv[1]].pos;
                    Point3d p2 = indexAttr[fv[2]].pos;
                    double m0 = verticeAttr[fv[0]].invMass;
                    double m1 = verticeAttr[fv[1]].invMass;
                    double m2 = verticeAttr[fv[2]].invMass;

                    Point3d corr0, corr1, corr2;

                    double stiffnessXX=0;
                    double stiffnessYY=0;
                    double stiffnessXY=0;

                        if(fD.cellZone==AD){
                          if(fD.cellLayer==L1)
                            stiffnessXX = 1*(norm(fD.a1)-norm(fD.a2));
                          else
                            stiffnessXX= 0.5*(norm(fD.a1)-norm(fD.a2));
                        }
                        else if(fD.cellZone==AB){
                            if(fD.cellLayer==L1)
                                stiffnessXX= 1*(norm(fD.a1)-norm(fD.a2));
                            else
                                stiffnessXX=0.5*(norm(fD.a1)-norm(fD.a2));
                        }
                        else if(fD.cellZone==MD){
                            if(fD.cellLayer==L1)
                                stiffnessXX=1*(norm(fD.a1)-norm(fD.a2));
                            else
                                stiffnessXX= 0.5*(norm(fD.a1)-norm(fD.a2));
                        }
                        else if(fD.cellZone==CZ){
                            stiffnessXX= 0;
                            stiffnessYY = 0;
                            stiffnessXY = 0;
                        }

                        else if(fD.cellZone==PZ){
                            stiffnessXX = 0;
                            stiffnessYY = 0;
                            stiffnessXY = 0;
                        }
                        else{
                            stiffnessXX = (norm(fD.a1) - norm(fD.a2))*3;
                            stiffnessYY = PBDstrainYYStiffness;
                            stiffnessXY = PBDstrainXYStiffness;
                        }

                   if (stiffnessXX>1)
                       stiffnessXX=1;

                   fD.stiffnessXX=stiffnessXX;

                    solve_StrainTriangleConstraint(p0, m0, p1, m1, p2, m2, fD.strainInvRestMat, stiffnessXX, stiffnessYY, stiffnessXY,
                                                   PBDstretchNormalize, PBDshearNormalize, corr0, corr1, corr2);

                    double  stiffness = PBDstrainStiffness;
                    corr0 *= stiffness;
                    corr1 *= stiffness;
                    corr2 *= stiffness;
                    verticeAttr[fv[0]].corrections["strain"] += corr0;
                    verticeAttr[fv[1]].corrections["strain"] += corr1;
                    verticeAttr[fv[2]].corrections["strain"] += corr2;
                    indexAttr[fv[0]].pos += corr0;
                    indexAttr[fv[1]].pos += corr1;
                    indexAttr[fv[2]].pos += corr2;
                }
            }
            }
   }







  ////---------------------------------MECHANICAL SOLVER-----------------------------//

  bool MechanicalSolver::step(bool timeToDebug, int stepCount) {

      CCStructure& cs = mesh->ccStructure("Triangulate"); //TISSUE TO GROWTH
      CCStructure& ts = mesh->ccStructure("Tissue");
      CCIndexDataAttr& indexAttr = mesh->indexAttr();
      VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");
      EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");
      FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");

  // Update Forces
      for(CCIndex v:cs.vertices())
         mechanicsProcess->calcForces(v, cs, indexAttr, faceAttr, edgeAttr, verticeAttr);

  //Calculate the stress (fD.force; sum of forces over the cell) and pass to triangulate;
      for (CCIndex f:ts.faces()){
          FaceData& fD=faceAttr[f];
          fD.cellVertices=faceVertices(ts,f);
          fD.force=0;
          for (uint i=0;i<fD.cellVertices.size();i++){ //calculate force
              CCIndex v = fD.cellVertices[i];
              fD.force +=verticeAttr[v].force;
          }
          for(CCIndex t:fD.cellFaces)//pass to triangulate
             faceAttr[t].force=fD.force;
      }

  // Update positions && clear corrections debug
      Point3d centerOfMasses;//momentumPreservation
      Point3d P_r;//linear
      Point3d L_r;//angular
      for(CCIndex v:cs.vertices()){
          centerOfMasses += indexAttr[v].pos * (1/verticeAttr[v].invMass);//center of masses before update positions
      }
      for(CCIndex v:cs.vertices()){
          solveEuler1stOrder(v, cs, indexAttr, verticeAttr, P_r, L_r,centerOfMasses);
          verticeAttr[v].corrections.clear();
      }


  //CALL SOLVE CONSTRAINTS
      //Debug
      double dStrainADIN=0; double dStrainADL1=0; double dStrainABIN=0; double dStrainABL1=0; double dStrainMDIN=0; double dStrainMDL1=0;
      solveConstraints(dStrainADIN,dStrainADL1,dStrainABIN,dStrainABL1,dStrainMDIN,dStrainMDL1,stepCount);

  //UPDATE VELOCITY
      for(CCIndex v:cs.vertices()) {
          VerticeData& vD = verticeAttr[v];
          if(velocityUpdate == "First Order")
              vD.velocity = (1.0 / Dt) * (indexAttr[v].pos - vD.prevPos) ;
          else if (velocityUpdate == "Second Order")
              vD.velocity = (1.0 / Dt) * (1.5 * indexAttr[v].pos - 2.0 * vD.prevPos + 0.5 * vD.lastPos);
          else
              throw(QString("Wrong velocity order"));
      }

  //CALL MOMENTUM PRESERVATION (in theory the constraints are build to conserve it)
      double meanCorrection=0;
      if(parm("Momentum Preservation K").toDouble() > 0)
          momentumPreservation(cs, indexAttr, verticeAttr, P_r, L_r, meanCorrection);
        //Debug
        Point3d linearMomentum, angularMomentum;
        for(uint i = 0; i < cs.vertices().size(); i++) {
            CCIndex v = cs.vertices()[i];
            VerticeData& vD = verticeAttr[v];
            for(auto deltap : vD.corrections) {
                Point3d corr = deltap.second;
                linearMomentum += 1./vD.invMass * corr;
                angularMomentum += indexAttr[v].pos ^ (1./vD.invMass * corr);
            }
        }

        //Polarity-->avoid error in division
        for (CCIndex f:ts.faces()){
            FaceData& fD=faceAttr[f];
            if(norm(fD.a1) == 0)
                fD.a1[1] = EPS;
            if(norm(fD.a2) == 0)
                fD.a2[0] = EPS;
            Point3d a1_inc, a2_inc;

  //CELL POLARITY (CMTs orientation model)
        if(parm("Polarity Method") == "Green Strain Tensor") {
            int i = 0;
            Matrix3d G;
            Matrix3d F;
            double gMax=0;
            double gMin=0;
            double gAngle=0;

            for (CCIndex ft:fD.cellFaces){ //do directly for perimeter vertices?¿
                std::vector<CCIndex> vs = faceVertices(cs, ft);
                Point3d x_p[3] = {verticeAttr[vs[0]].prevPos,verticeAttr[vs[1]].prevPos,verticeAttr[vs[2]].prevPos}; ////// NW
                Point3d X_p[3] = {indexAttr[vs[0]].pos, indexAttr[vs[1]].pos, indexAttr[vs[2]].pos};
                F = DefGradient(x_p, X_p);
                F[2][2] = 1;
                G += 0.5 * (transpose(F) * F - F.identity()); //only if fD.type==membrane?¿
                i++;
            }
            G /= i;
            if((G[0][0] - G[1][1]) != 0) {
                gAngle = (0.5 * atan(2 * G[0][1] / (G[0][0] - G[1][1]) )) ;
                double left = (G[0][0] + G[1][1]) / 2;
                double right = sqrt(pow((G[0][0] - G[1][1]) / 2, 2) + pow(G[0][1]/2, 2));
                gMax = left + right;
                gMin = left - right;
            }
            if(G[0][0] < G[1][1]){
                if(gAngle > -(M_PI/4) && gAngle < (M_PI/4))
                   gAngle += (M_PI/2);
            }

            if(gMax > MFspringTresh) //treshold
                a1_inc = RotateCounterClock(Point3d(gMax,0,0), gAngle); //value gMax which is rotated according to the gAngle

            if(gMin > MFspringTresh)
                a2_inc = RotateCounterClock(Point3d(gMin,0,0), gAngle+(M_PI/2));

        }
        else if (parm("Polarity Method") == "Cell Axis") {
           std::vector<Point3d> box; //rectangle which englobes the exterior of the cell
           Point3d axisMin, axisMax; //shortest and longest axis of the rectangle
           std::vector<Point3d> points;
           for(CCIndex v:faceVertices(ts,f))
               points.push_back(indexAttr[v].pos);
           MBR(points, box);
           axisMax = box[0] - box[1];
           axisMin = box[1] - box[2];
           if(norm(axisMax) < norm(axisMin)) {
               Point3d tmp = axisMin;
               axisMin = axisMax;
               axisMax = tmp;
           }

          if(MFRORate > 0) {
             if(findClosestLineToLine(axisMax,fD.a1, fD.a2) == fD.a1) { //if fD.a1 is closer to axisMAx then takes its direction// *norm(fD.a1)--previous value
                fD.a1 = axisMax/norm(axisMax) * norm(fD.a1); //fD.a1 and fD.a2 = 0 when initialize -->no polarity
                fD.a2 = axisMin/norm(axisMin) * norm(fD.a2);
            } else {
                fD.a2 = axisMax/norm(axisMax) * norm(fD.a2);
                fD.a1 = axisMin/norm(axisMin) * norm(fD.a1);
            }

            for(CCIndex e : ts.incidentCells(f,1)) {
                EdgeData& eD = edgeAttr[e];
                auto eb = ts.edgeBounds(e);
                Point3d edge = indexAttr[eb.first].pos - indexAttr[eb.second].pos;
                Point3d dv = verticeAttr[eb.first].velocity - verticeAttr[eb.second].velocity; //velocity after update
                double prevLength=norm(verticeAttr[eb.first].prevPos-verticeAttr[eb.second].prevPos);
                double strainRate = (edge * dv) / prevLength;
                double length=norm(edge);
                double strain = (length - eD.restLength) / eD.restLength;
                if(strain< MFspringTresh) //if the edge hasn't deform enough is not considered--in fact is not going to elongate--plastic effect
                    continue;
                edge /= norm(edge);
                a1_inc += fD.a1/norm(fD.a1) * abs(fD.a1/norm(fD.a1) * edge * strainRate);
                a2_inc += fD.a2/norm(fD.a2) * abs(fD.a2/norm(fD.a2) * edge * strainRate);
            }
        }

      }
       else
           throw(QString("Wrong Polarity Method"));

        //Update Polarity (Common)
        Point3d a1U = norm(fD.a1) > 0 ? fD.a1 / norm(fD.a1) : Point3d(0, 0, 0);
        Point3d a2U = norm(fD.a2) > 0 ? fD.a2 / norm(fD.a2) : Point3d(0, 0, 0);
        Point3d degr_a1 =
            (MFRORate > 0) ? MFdegradation * a1U * norm(fD.a1): Point3d(0, 0, 0);
        Point3d degr_a2 =
            (MFRORate > 0) ? MFdegradation * a2U * norm(fD.a2) : Point3d(0, 0, 0);

        fD.a1 += (a1_inc * MFRORate - degr_a1) * Dt;
        fD.a1_inc=a1_inc;
        fD.a2 += (a2_inc * MFRORate - degr_a2) * Dt;
        if(norm(fD.a1) > 1)
            fD.a1 = shortenLength(fD.a1, norm(fD.a1) - 1);
        if(norm(fD.a2) > 1)
            fD.a2 = shortenLength(fD.a2, norm(fD.a2) - 1);
        if(norm(fD.a1) < norm(fD.a2)) {
            Point3d tmp = fD.a1;
            fD.a1 = fD.a2;
            fD.a2 = tmp;
        }
        if(norm(fD.a1) == 0) //again to avoid errors in divisions
            fD.a1[1] = EPS;
        if(norm(fD.a2) == 0)
            fD.a2[0] = EPS;

      }

      return true;

  }


////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////UPDATE RESTSTATE/////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////

bool restUpdate::initialize(QWidget *parent, std::string debFileName)
{
  mesh = currentMesh();
  if(!mesh)
    throw(QString("CellDiskGrowth::initialize No current mesh"));

//PARAMETERS
  Dt = parm("Dt").toDouble();
  springThresh = parm("Strain Threshold for Growth").toDouble(); //same as in MF treshold for reorientation, only reorientate if the edge grows
  wallsGrowthRate = parm("Walls Growth Rate").toDouble();
  divL1 = parm("kL1").toDouble();
  n1L1 = parm("n1L1").toDouble();
  kL1 = parm("kmL1").toDouble();
  n2L1 = parm("n2L1").toDouble();
  iniL1 = parm("bL1").toDouble();
  divIN = parm("kIN").toDouble();
  n1IN = parm("n1IN").toDouble();
  kIN = parm("kmIN").toDouble();
  n2IN = parm("n2IN").toDouble();
  iniIN = parm("bIN").toDouble();

  return true;
}

bool restUpdate::step(bool timeToDebug,int stepCount)
{
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");
    EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");
    FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
    CCStructure& cs = mesh->ccStructure("Triangulate"); //TISSUE TO GROWTH
    CCStructure& ts = mesh->ccStructure("Tissue");

//DISTANCE CONSTRAINT UPDATE (restLength) -->AUXIN EFFECT
    //elastic-plastic behaivour, the rest length is only updated if the strain is over a treshold (consider edge length)
    for(CCIndex e:cs.edges()) {
        EdgeData& eD = edgeAttr[e];
        auto eb=cs.edgeBounds(e);
        double length=norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos);
            if (eD.restLength==0){ //for new edges that comes from retriangulation after division
                eD.restLength = length;
            }

            double scalefact=225;

            if(eD.wall){
                double auxinEffect;
                CCIndex es=edgeAttr[e].sameEdge;
                auto ecb=ts.cobounds(es);
                if(ecb.size()==1){
                  CCIndex cb=*ecb.begin();
                  double auxinConc=scalefact*(faceAttr[cb].auxin/indexAttr[cb].measure);
                  auxinEffect=(pow(auxinConc*divL1,n1L1)*(pow(auxinConc,n2L1)/(pow(kL1,n2L1)+pow(auxinConc,n2L1))))+iniL1;
                } else if(ecb.size()==2){
                   std::vector<CCIndex> cb(ecb.begin(),ecb.end());
                   double auxinConc_1, auxinConc_2;
                   auxinConc_1=faceAttr[cb[0]].auxin/indexAttr[cb[0]].measure;
                   auxinConc_2=faceAttr[cb[1]].auxin/indexAttr[cb[1]].measure;
                   double auxinConc=scalefact*0.5*(auxinConc_1+auxinConc_2); //half effect for edges attached to two cells

                   if(faceAttr[cb[0]].cellLayer==L1 && faceAttr[cb[1]].cellLayer==L1){
                      auxinEffect=(pow(auxinConc*divL1,n1L1)*(pow(auxinConc,n2L1)/(pow(kL1,n2L1)+pow(auxinConc,n2L1))))+iniL1;
                   }
                     else{
                       auxinEffect=(pow(auxinConc*divIN,n1IN)*(pow(auxinConc,n2IN)/(pow(kIN,n2IN)+pow(auxinConc,n2IN))))+iniIN;
                   }

                } else
                    throw(QString("RestUpdate: more than 2 faces for 1 edge"));

              double strain = (length - eD.restLength) / eD.restLength;

              if(strain >= springThresh) {
                 double strainDiff = strain - springThresh;

                 eD.restLength += wallsGrowthRate * auxinEffect * strainDiff * Dt;

                 if(eD.restLength > length){
                    mdxInfo << "WARNING: restLength=length" << endl;
                    eD.restLength = length;
                   }
              }
            }
            else
              eD.restLength = length;
            }

//SHAPE CONSTRAINT UPDATE (restCm)
    for(CCIndex f:ts.faces()) {
        FaceData& fD = faceAttr[f];
        fD.restX0.clear();
        fD.cellVertices= faceVertices(ts, f);
        std::vector<double> invMasses;
        for(CCIndex v : fD.cellVertices) {
            fD.restX0.push_back(Point3d(indexAttr[v].pos));
            invMasses.push_back(verticeAttr[v].invMass);
        }
        init_ShapeMatchingConstraint(fD.restX0, invMasses, fD.restX0.size(),  fD.restCm, fD.invRestMat);
    }


  return true;
}


////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////DIVISION/////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////

////---------------------------------SUBDIVIDE----------------------------------------//
bool MySubdivide::initialize(QWidget *parent)
  {
      mesh = currentMesh();
      if(!mesh)
        throw(QString("General::initialize No current mesh"));

    return true;
}

bool MySubdivide::step(Dimension dim, CCStructure &cs, CCStructure::SplitStruct &ss,
                         CCIndex otherP, CCIndex otherN, double interpPos, bool triangulate, std::vector<CCIndexVec> edgeDivisions, std::vector<edgeParData> edatas,
                         double KPin1auxinflux, double KPin1MF, double KPin1MFauxinflux,std::vector<edgeStruct>  edgeData)
  {
      CCIndexDataAttr& indexAttr = mesh->indexAttr();
      FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
      EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");
      VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");

//DIMENSION 1
      if(dim ==1) {
        //otherP, otherN bound vertices
        //ss.membrane new vertex
        //ss.childP, ss.childN new edges
        //ss.parent parent edge
        CCIndexData &iP = indexAttr[ss.childP], &iN = indexAttr[ss.childN]/*, &iPar = indexAttr[ss.parent]*/;
        EdgeData &eP = edgeAttr[ss.childP], &eN = edgeAttr[ss.childN], &ePar = edgeAttr[ss.parent];
        VerticeData &vOP = verticeAttr[otherP], &vON = verticeAttr[otherN], &vM = verticeAttr[ss.membrane];

        //Length
        auto endP=cs.edgeBounds(ss.childP);
        iP.measure =norm(indexAttr[endP.first].pos - indexAttr[endP.second].pos);
        auto endN=cs.edgeBounds(ss.childN);
        iN.measure =norm(indexAttr[endN.first].pos - indexAttr[endN.second].pos);

        double totalLength=iP.measure+iN.measure;

        //mass and velocity
        if(vOP.invMass<1 && vON.invMass<1)
            vM.invMass=1/BIG_VAL;
        vM.velocity=(vOP.velocity+vON.velocity)*0.5;

        if(vOP.bottom && vON.bottom)
            vM.bottom=true;

        //to split in triangulate
        iP.selected=true;
        iN.selected=true;

        //EDGES
        eP.intercellularAuxin=ePar.intercellularAuxin*(iP.measure/totalLength);
        eN.intercellularAuxin=ePar.intercellularAuxin*(iN.measure/totalLength);

        for(uint i = 0 ; i < edatas.size() ; i++) {
           eP.pin1[edatas[i].label]=edatas[i].pin*(iP.measure/totalLength);
           eN.pin1[edatas[i].label]=edatas[i].pin*(iN.measure/totalLength);
           eP.edge_export_prev[edatas[i].label]=edatas[i].edge_export_prev*(iP.measure/totalLength);
           eN.edge_export_prev[edatas[i].label]=edatas[i].edge_export_prev*(iN.measure/totalLength);
        }
      }

//DIMENSION 2
      if(dim ==2) {
        //ss.membrane new edge
        //ss.childP, ss.childN new faces
        //ss.parent parent face

        CCIndexData &iP = indexAttr[ss.childP], &iN = indexAttr[ss.childN], &iPar = indexAttr[ss.parent], &iM = indexAttr[ss.membrane];
        FaceData &fP = faceAttr[ss.childP], &fN = faceAttr[ss.childN], &fPar = faceAttr[ss.parent];
        EdgeData &eM = edgeAttr[ss.membrane];

        //LABEL, SELECTED, AREA
         iP.label=LABEL;
         LABEL++;
         iN.label=LABEL;
         LABEL++;
         iP.selected = iN.selected = iPar.selected;
         iP.measure = (1. - interpPos) * iPar.measure;
         iN.measure = interpPos * iPar.measure;

         //CELL TYPE
         fP.cellSide=fN.cellSide=fPar.cellSide;
         fP.cellLayer=fN.cellLayer=fPar.cellLayer;
         fP.cellZone=fN.cellZone=fPar.cellZone;
         fP.cellPosition=fN.cellPosition=fPar.cellPosition;
         fP.prevArea=iPar.measure/2;
         fN.prevArea=iPar.measure/2;

         //CHEMICAL
         //face data
         double totalArea=iP.measure+iN.measure;
         fP.auxin=fPar.auxin*(iP.measure/totalArea);
         fP.pin1=fPar.pin1*(iP.measure/totalArea);
         fN.auxin=fPar.auxin*(iN.measure/totalArea);
         fN.pin1=fPar.pin1*(iN.measure/totalArea);
         fP.auxinFluxVector=fPar.auxinFluxVector;
         fN.auxinFluxVector=fPar.auxinFluxVector;
         fP.strain=fPar.strain;
         fN.strain=fPar.strain;

         //edge data
         double maxPin=0;
         for(CCIndex e:cs.incidentCells(ss.childP,1)){
             EdgeData& eD = edgeAttr[e];
             CCIndex fpos, fneg;
             std::set<CCIndex> edgeFace =cs.incidentCells(e,2);
             if(edgeFace.size()==1){
                fpos=*edgeFace.begin();
             }
              else if(edgeFace.size()==2){
                  std::vector<CCIndex> cb(edgeFace.begin(),edgeFace.end());
                  if(indexAttr[cb[0]].label==iP.label){
                     fpos=cb[0];
                     fneg=cb[1];
                  }
                  else if(indexAttr[cb[1]].label==iP.label){
                      fpos=cb[1];
                      fneg=cb[0];
                  } else
                      mdxInfo<<"WRONG ORIENTATION SUBDV"<<endl;
           }

           auto eb=cs.edgeBounds(e);
           Point3d p1=indexAttr[eb.first].pos;
           Point3d p2=indexAttr[eb.second].pos;
           for (uint k = 0; k < edgeData.size(); k++){
               if((p1==edgeData[k].fBound && p2==edgeData[k].sBound) || (p2==edgeData[k].fBound && p1==edgeData[k].sBound)){
                   eD.intercellularAuxin=edgeData[k].intercellularAuxin;
                   eD.pin1[indexAttr[fpos].label]=edgeData[k].posPin1;
                   eD.pin1[indexAttr[fneg].label]=edgeData[k].negPin1;
                   eD.edge_export_prev[indexAttr[fpos].label]=edgeData[k].posPrevExport;
                   eD.edge_export_prev[indexAttr[fneg].label]=edgeData[k].negPrevExport;

                   if( eD.pin1[indexAttr[fpos].label]>maxPin)
                       maxPin=eD.pin1[indexAttr[fpos].label];
                   if( eD.pin1[indexAttr[fneg].label]>maxPin)
                       maxPin=eD.pin1[indexAttr[fneg].label];
                   break;
                 }
           }
         }

         for(CCIndex e:cs.incidentCells(ss.childN,1)){
             EdgeData& eD = edgeAttr[e];
             CCIndex fpos, fneg;
             std::set<CCIndex> edgeFace =cs.incidentCells(e,2);
             if(edgeFace.size()==1){
                fpos=*edgeFace.begin();
             }
              else if(edgeFace.size()==2){
                  std::vector<CCIndex> cb(edgeFace.begin(),edgeFace.end());
                  if(indexAttr[cb[0]].label==iN.label){
                     fpos=cb[0];
                     fneg=cb[1];
                  }
                  else if(indexAttr[cb[1]].label==iN.label){
                      fpos=cb[1];
                      fneg=cb[0];
                  } else
                      mdxInfo<<"WRONG ORIENTATION SUBDV"<<endl;

             }

             auto eb=cs.edgeBounds(e);
             Point3d p1=indexAttr[eb.first].pos;
             Point3d p2=indexAttr[eb.second].pos;
             for (uint k = 0; k < edgeData.size(); k++){
                 if((p1==edgeData[k].fBound && p2==edgeData[k].sBound) || (p2==edgeData[k].fBound && p1==edgeData[k].sBound)){
                     eD.intercellularAuxin=edgeData[k].intercellularAuxin;
                     eD.pin1[indexAttr[fpos].label]=edgeData[k].posPin1;
                     eD.pin1[indexAttr[fneg].label]=edgeData[k].negPin1;
                     eD.edge_export_prev[indexAttr[fpos].label]=edgeData[k].posPrevExport;
                     eD.edge_export_prev[indexAttr[fneg].label]=edgeData[k].negPrevExport;

                     if( eD.pin1[indexAttr[fpos].label]>maxPin)
                         maxPin=eD.pin1[indexAttr[fpos].label];
                     if( eD.pin1[indexAttr[fneg].label]>maxPin)
                         maxPin=eD.pin1[indexAttr[fneg].label];
                     break;
                 }
            }
         }


         //CHEMICAL FOR NEW MEMBRANE
         auto eb = cs.edgeBounds(ss.membrane);
         Point3d edge = indexAttr[eb.second].pos-indexAttr[eb.first].pos;
         Point3d outwardNormal = Rotate90(edge);
         Point3d midPoint= indexAttr[eb.first].pos+(edge)/2.;

         std::vector<CCIndex> childFaces;
         childFaces.push_back(ss.childP);
         childFaces.push_back(ss.childN);

         for (CCIndex f:childFaces){
             FaceData& fD=faceAttr[f];
             CCIndexData &iD=indexAttr[f];
             int label =iD.label;

             //Centroid
             std::vector<CCIndex> indices=faceVertices(cs,f);
             int size=indices.size();
             std::vector<Point3d> vertices;
             for (CCIndex v: indices){
                  vertices.push_back(indexAttr[v].pos);
              }
              Point3d centroid=Centroid(vertices,size);
              fD.newcentroid=centroid;

             Point3d insideFace=fD.newcentroid-midPoint;
             double angle =acos((insideFace*outwardNormal)/(norm(insideFace)*norm(outwardNormal)));
             if (angle<(M_PI/2))
                outwardNormal *=-1;

              //AUXIN
              double auxinImpactValue=0;
              if(norm(fD.auxinFluxVector) > 0){
                  auxinImpactValue = ((fD.auxinFluxVector*outwardNormal)/(norm(fD.auxinFluxVector)*norm(outwardNormal)));
                  fD.auxinFluxVector=0;
              }
              if(auxinImpactValue < 0)
                  auxinImpactValue = 0;

              //MF
              double MFvalue = ((fD.strain*outwardNormal)/(norm(fD.strain)*norm(outwardNormal)));
              if(norm(fD.strain)==0)
                  MFvalue=0;
              if(MFvalue < 0)
                  MFvalue *= -1;

              double pinInit=((auxinImpactValue*KPin1auxinflux)+(MFvalue*KPin1MF)+(auxinImpactValue*MFvalue*KPin1MFauxinflux))/
                      (KPin1auxinflux+KPin1MF+KPin1MFauxinflux);

              eM.pin1[label] +=pinInit*maxPin;

              fD.auxinFluxVector += ((fD.newcentroid-midPoint)/norm(fD.newcentroid-midPoint)) * -1*eM.pin1[label] * 0.01; //?¿

              fD.pin1 -=eM.pin1[label];
              if(fD.pin1<0)
                 fD.pin1=0;
         }
         iM.selected=true; //to split in triangulate
     }


      return true;
}


////---------------------------------DIVIDE------------------------------------//

bool MyDivide::initialize(QWidget* parent)
{
    mesh = currentMesh();
    if(!mesh)
      throw(QString("MyDivide::initialize No current mesh"));

    if(!getProcess(parm("Subdivide Process"), subdivideProcess))
      throw(QString("Divide::initialize Cannot make Subdivide Process"));
    subdivideProcess->initialize(parent);

//PARAMETERS
    CellMaxArea = parm("Cell Max Area").toDouble();
    ABdiv = parm("AB Div Size").toDouble();
    CellWallMin = parm("Cell Wall Min").toDouble();
    CellWallSample = parm("Cell Wall Sample").toDouble();
    CellPinch = parm("Cell Pinch").toDouble();
    CellMaxPinch = parm("Cell Max Pinch").toDouble();
    Manual= stringToBool(parm("Manual Division"));

  return true;
}

bool MyDivide::step(double KPin1auxinflux, double KPin1MF, double KPin1MFauxinflux)
{

    CCStructure& cs = mesh->ccStructure("Tissue");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
    EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");

//FIND CELLS TO DIVIDE
    CCIndexVec toDivide;
    if(Manual){
      for(CCIndex f:cs.faces()) {
         if(indexAttr[f].selected){
            toDivide.push_back(f);
      }
      }
    }
    else{
    for(CCIndex f:cs.faces()) {
      if((indexAttr[f].measure > CellMaxArea && faceAttr[f].cellZone != AB) || indexAttr[f].measure > ABdiv*CellMaxArea){
        toDivide.push_back(f);
        indexAttr[f].selected=true; //select the faces to retriangulate (the one than didvides and the neighbours)
        for(CCIndex e:cs.incidentCells(f,1)){
            indexAttr[e].selected=true;
            for(CCIndex fi:cs.incidentCells(e,2))
                indexAttr[fi].selected=true;
        }
      } else
         faceAttr[f].auxinFluxVector=0;
    }
    }
    if(toDivide.empty())
      return false;

    // Structure representing a division wall in 2D
    struct DivWall2d
    {
    struct {
    Point3d pos;    // position of endpoint
    CCIndex vA, vB; // vertices of wall to be split
    double sfrac;   // fraction of position along wall from A to B
    }
    endpoints[2];
    };

    DivWall2d divWall;
    for(const CCIndex &cell : toDivide) {

//FIND WALLS FOR DIVISION
    for(uint i = 0 ; i < 2 ; i++)
      divWall.endpoints[i].vA = divWall.endpoints[i].vB = CCIndex::UNDEF;

    // Check dimension
    Dimension dimension = cs.dimensionOf(cell);
    if(dimension != 2)
      throw(QString("findCellDiv2d: Can only split a two-dimensional cell!"));

    // get the vertex position in the cell to divide
    CCStructure::CellTuple ct(cs,cell);
    std::vector<Point3d> pos;
    CCIndex firstV = ct[0];
    do {
      pos.push_back(indexAttr[ct[0]].pos);
      ct.flip(0,1);
    } while(ct[0] != firstV);
    uint numVertices = pos.size(); //number  of vertex

    // get centroid and normal //cdata.centroid, cdata.normal
      auto cdata = polygonCentroidData(pos);

//L1 Stress Throught Centroid
      if (faceAttr[cell].cellLayer==L1) { //Stress Throught Centroid
          //division plane
          Point3d divVector=cdata.centroid + faceAttr[cell].force;
          CCIndexVec edgesToDivide;
          CCIndexVec boundVertices;
          Point3dVec intersectPoints;

          for(CCIndex e:cs.incidentCells(cell,1)){
             auto eb=cs.edgeBounds(e);
              Point3d intersect {0.0,0.0,0.0};
              if(lineSegmentIntersection(cdata.centroid,divVector,indexAttr[eb.first].pos,indexAttr[eb.second].pos, intersect)){
                 edgesToDivide.push_back(e);
                 boundVertices.push_back(eb.first);
                 boundVertices.push_back(eb.second);
                 intersectPoints.push_back(intersect);
              }
          }

          if(edgesToDivide.size()==2){
            divWall.endpoints[0] = { intersectPoints[0],  boundVertices[0],  boundVertices[1],
                                    norm(intersectPoints[0] - indexAttr[boundVertices[0]].pos) / norm(indexAttr[boundVertices[1]].pos - indexAttr[boundVertices[0]].pos) };
            divWall.endpoints[1] = { intersectPoints[1], boundVertices[2], boundVertices[3],
                                     norm(intersectPoints[1] - indexAttr[boundVertices[2]].pos) / norm(indexAttr[boundVertices[3]].pos - indexAttr[boundVertices[2]].pos) };
          }

          else{
              mdxInfo<<"STRESS DIVISION FAILED"<<endl;

              double shortestDist = HUGE_VAL;
              for(uint i = 0 ; i < numVertices ; i++) {
                uint i1 = (i+1) % numVertices;
                Point3d u1 = pos[i], u2 = pos[i1], u1u2 = u2 - u1;
                double u1u2l = norm(u1u2);
                Point3d u1u2d = u1u2 / u1u2l;
                double ulEnd = u1u2l - CellWallMin;
                for(double ul = CellWallMin ; ul <= ulEnd ; ul += CellWallSample) {
                  Point3d u = u1 + ul * u1u2d, uc = cdata.centroid - u;
                  Point3d nu = cdata.normal ^ uc;
                  double centroidDist = norm(uc);

                  // Now we check the upcoming walls for an intersection
                  CCStructure::CellTuple ct1(ct);
                  ct1.flip(0,1);
                  for(uint j = i+1; j < numVertices; j++) {
                    uint j1 = (j+1) % numVertices;
                    Point3d v;
                    double s;
                    if(planeLineIntersect(u, nu, pos[j], pos[j1], s, v) and (s >= 0. && s <= 1.)) {

                      double dist = norm(u - v);

                      // Check that the line actually goes through the centroid (needed if the cell is not convex)
                      if(dist < centroidDist) continue;


                        Point3d v1v2 = pos[j1] - pos[j];
                        double v1v2l = norm(v1v2);
                        v1v2 /= v1v2l;

                        // len is the distance along the edge from v to v1
                        double len = (v - pos[j]) * v1v2, len1 = len;
                        if(v1v2l <= CellWallMin * 2.0)
                          // If the edge is too short, we put v at its midpoint
                          len1 = (v1v2l/2.0);
                        else if(len < CellWallMin)
                          // If v is too close to v1, we put it at the minimum distance
                          len1 = CellWallMin;
                        else if(len > v1v2l - CellWallMin)
                          // If v is too close to v2, we put it at the minimum distance
                          len1 = v1v2l - CellWallMin;

                        v = pos[j] + v1v2 * len1;
                        dist += fabs(len - len1);

                      if(dist < shortestDist) {
                        shortestDist = dist;
                        divWall.endpoints[0] = { u,  ct[0],  ct.other(0), norm(u - pos[i]) / norm(pos[i1] - pos[i]) };
                        divWall.endpoints[1] = { v, ct1[0], ct1.other(0), norm(v - pos[j]) / norm(pos[j1] - pos[j]) };
                      }
                    }
                    ct1.flip(0,1);
                  }
                }
                ct.flip(0,1);
              }
            }

    }
//Inner Shortest Path
      else {
      double shortestDist = HUGE_VAL;
      for(uint i = 0 ; i < numVertices ; i++) {
        uint i1 = (i+1) % numVertices;
        Point3d u1 = pos[i], u2 = pos[i1], u1u2 = u2 - u1;
        double u1u2l = norm(u1u2);
        Point3d u1u2d = u1u2 / u1u2l;
        double ulEnd = u1u2l - CellWallMin;
        for(double ul = CellWallMin ; ul <= ulEnd ; ul += CellWallSample) {
          Point3d u = u1 + ul * u1u2d, uc = cdata.centroid - u;
          Point3d nu = cdata.normal ^ uc;
          double centroidDist = norm(uc);

          // Now we check the upcoming walls for an intersection
          CCStructure::CellTuple ct1(ct);
          ct1.flip(0,1);
          for(uint j = i+1; j < numVertices; j++) {
            uint j1 = (j+1) % numVertices;
            Point3d v;
            double s;
            if(planeLineIntersect(u, nu, pos[j], pos[j1], s, v) and (s >= 0. && s <= 1.)) {

              double dist = norm(u - v);

              // Check that the line actually goes through the centroid
              // (needed if the cell is not convex)
              if(dist < centroidDist) continue;

                    Point3d v1v2 = pos[j1] - pos[j];
                    double v1v2l = norm(v1v2);
                    v1v2 /= v1v2l;

                    // len is the distance along the edge from v to v1
                    double len = (v - pos[j]) * v1v2, len1 = len;
                    if(v1v2l <= CellWallMin * 2.0)
                      // If the edge is too short, we put v at its midpoint
                      len1 = (v1v2l/2.0);
                    else if(len < CellWallMin)
                      // If v is too close to v1, we put it at the minimum distance
                      len1 = CellWallMin;
                    else if(len > v1v2l - CellWallMin)
                      // If v is too close to v2, we put it at the minimum distance
                      len1 = v1v2l - CellWallMin;

                    v = pos[j] + v1v2 * len1;
                    dist += fabs(len - len1);

                  if(dist < shortestDist) {
                    shortestDist = dist;
                    divWall.endpoints[0] = { u,  ct[0],  ct.other(0), norm(u - pos[i]) / norm(pos[i1] - pos[i]) };
                    divWall.endpoints[1] = { v, ct1[0], ct1.other(0), norm(v - pos[j]) / norm(pos[j1] - pos[j]) };
                  }
                }
                ct1.flip(0,1);
              }
            }
            ct.flip(0,1);
          }
        }

        if(divWall.endpoints[0].vA.isPseudocell() or divWall.endpoints[0].vB.isPseudocell() or
           divWall.endpoints[1].vA.isPseudocell() or divWall.endpoints[1].vB.isPseudocell()) {
          mdxInfo << "CellTissue::findCellDiv Unable to find divide wall for cell:"
                  << to_string(cell).c_str() << ", label:" << indexAttr[cell].label
                  << ", CellWallMin or CellWallSample too high?" <<  endl;
          return false;
        }

// Now we pinch the walls.
        for(uint i = 0 ; i < 2 ; i++) {
          auto &ep = divWall.endpoints[i];
          CCIndex edge = cs.join(ep.vA,ep.vB);
          // Don't pinch a given cell wall if it is external
          if(!cs.onBorder(edge)) {
            Point3d centroidDir = cdata.centroid - ep.pos;
            double cdl = norm(centroidDir);
            centroidDir /= cdl;

            Point3d posA = indexAttr[ep.vA].pos, posB = indexAttr[ep.vB].pos;
            double distA = norm(ep.pos - posA), distB = norm(ep.pos - posB);
            double minD = std::min(distA , distB);
            double pinchAmt = std::min(minD , cdl) * CellPinch;
            pinchAmt = std::min(pinchAmt , CellMaxPinch);

            ep.pos += centroidDir * pinchAmt;
          }
        }

//DIVIDE WALLS
      CCIndex ep[2];
       CCIndexVec edgeDivision;
       std::vector<CCIndexVec> edgeDivisions;

      for(int i = 0 ; i < 2 ; i++) {
        CCIndex edge = edgeBetween(cs, divWall.endpoints[i].vA, divWall.endpoints[i].vB);

        std::vector<edgeParData> edatas;
        for(CCIndex f:cs.incidentCells(edge,2)){
            edgeParData edata;
            CCIndexData &iD=indexAttr[f];
            edata.label=iD.label;
            edata.pin=edgeAttr[edge].pin1[iD.label];
            edata.edge_export_prev=edgeAttr[edge].edge_export_prev[iD.label];
            edatas.push_back(edata);
        }

        SplitStruct ss(edge);
        CCIndexFactory.fillSplitStruct(ss);

        ep[i] = ss.membrane;

        bool neg = cs.ro(cell, edge) == ccf::NEG;
        cs.splitCell(ss);
        // Update the position ourselves since the subdivider may need it
        indexAttr[ep[i]].pos = divWall.endpoints[i].pos;


        double empty2=0;
        std::vector<edgeStruct>  empty3;

         subdivideProcess->step(1, cs, ss, divWall.endpoints[i].vA, divWall.endpoints[i].vB, //call twice
                                neg ? divWall.endpoints[i].sfrac : (1.0 - divWall.endpoints[i].sfrac),false, edgeDivisions,edatas, empty2, empty2,empty2,empty3);

         edgeDivision.push_back(ss.childP);
         edgeDivision.push_back(ss.childN);
         edgeDivision.push_back(ss.parent);
         edgeDivisions.push_back(edgeDivision);
         edgeDivision.clear();
      }

      std::vector<edgeStruct> initTriEdges;

      for(CCIndex e:cs.incidentCells(cell,1)){
          edgeStruct initTriEdge;
          EdgeData& eD=edgeAttr[e];
          auto eb=cs.edgeBounds(e);
          Point3d p1=indexAttr[eb.first].pos;
          Point3d p2=indexAttr[eb.second].pos;
          initTriEdge.intercellularAuxin=eD.intercellularAuxin;
          initTriEdge.fBound=p1;
          initTriEdge.sBound=p2;

          CCIndex fpos, fneg;
          std::set<CCIndex> edgeFace =cs.incidentCells(e,2);
          if(edgeFace.size()==1){
             fpos=*edgeFace.begin();
          }
           else if(edgeFace.size()==2){
               std::vector<CCIndex> cb(edgeFace.begin(),edgeFace.end());
               if(indexAttr[cb[0]].label==indexAttr[cell].label){
                  fpos=cb[0];
                  fneg=cb[1];
               }
               else if(indexAttr[cb[1]].label==indexAttr[cell].label){
                   fpos=cb[1];
                   fneg=cb[0];
               } else
                   mdxInfo<<"MYDIVIDE: WRONG ORIENTATION TRIANGULATES"<<endl;
          }
          initTriEdge.posPin1=eD.pin1[indexAttr[fpos].label];
          initTriEdge.negPin1=eD.pin1[indexAttr[fneg].label];
          initTriEdge.posPrevExport=eD.edge_export_prev[indexAttr[fpos].label];
          initTriEdge.posPrevExport=eD.edge_export_prev[indexAttr[fneg].label];

          initTriEdges.push_back(initTriEdge);
      }

//DIVIDE CELL
      CCStructure::SplitStruct ss(cell);
      CCIndexFactory.fillSplitStruct(ss);
      cs.splitCell(ss, +ep[0] -ep[1]);
      updateFaceGeometry(cs, indexAttr, ss.childP);
      updateFaceGeometry(cs, indexAttr, ss.childN);
      double s = indexAttr[ss.childP].measure / (indexAttr[ss.childP].measure + indexAttr[ss.childN].measure);

      std::vector<edgeParData> empty;
      subdivideProcess->step(2, cs, ss, CCIndex(), CCIndex(), s, false, edgeDivisions, empty, KPin1auxinflux, KPin1MF, KPin1MFauxinflux,initTriEdges);
      edgeDivisions.clear();
  }
    return true;
}

  ////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////PREVIOUS PROCESSES///////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////


  ////---------------------------------CELL TYPE----------------------------------------//
  bool SetCellType::step() {
        mesh = currentMesh();
        if(!mesh or mesh->file().isEmpty())
            throw(QString("SetCellType::step No current mesh, cannot rewind"));

        CCStructure& cs = mesh->ccStructure("Tissue");
        CCIndexDataAttr& indexAttr = mesh->indexAttr();
        FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");

  //INITIALIZE MAPS
        CCIndexIntAttr& cellside = mesh->signalAttr<int>("CellSide");
        CCIndexIntAttr& celllayer = mesh->signalAttr<int>("CellLayer");
        CCIndexIntAttr& cellzone = mesh->signalAttr<int>("CellZone");
        CCIndexIntAttr& cellposition = mesh->signalAttr<int>("CellPosition");
        CCIndexIntAttr& initialauxin = mesh->signalAttr<int>("InitialAuxin");

  //PARAMETERS
        SetSide = stringToBool(parm("Set Side"));
        SetLayer = stringToBool(parm("Set Layer"));
        SetZone = stringToBool(parm("Set Zone"));
        SetPosition = stringToBool(parm("Set Position"));
        SetInitialAuxin = stringToBool(parm("Set Initial Auxin"));
        CellSide = stringToCellSide(parm("Cell Side"));
        CellLayer = stringToCellLayer(parm("Cell Layer"));
        CellZone = stringToCellZone(parm("Cell Zone"));
        CellPosition = stringToCellPosition(parm("Cell Position"));
        InitialAuxin = stringToInitialAuxin(parm("Initial Auxin"));
        i= parm("Index").toInt();

        if(SetSide){
        for (CCIndex f:cs.faces()){
            FaceData& fD = faceAttr[f];
            auto fv=faceVertices(cs,f);
                 for (uint j=0;j<fv.size();j++){
                     if(indexAttr[fv[j]].selected){
                         fv.erase(fv.begin()+j);
                         for(CCIndex k:fv){
                             if(indexAttr[k].selected){
                                 fD.cellSide = CellSide;
                                 break;
                             }
                         }
                     }
                 }
        }
        }

        if(SetLayer){
        for (CCIndex f:cs.faces()){
            FaceData& fD = faceAttr[f];
            auto fv=faceVertices(cs,f);
                 for (uint j=0;j<fv.size();j++){
                     if(indexAttr[fv[j]].selected){
                         fv.erase(fv.begin()+j);
                         for(CCIndex k:fv){
                             if(indexAttr[k].selected){
                                 fD.cellLayer = CellLayer;
                                 break;
                             }
                         }
                     }
                 }
        }
        }

        if(SetZone){
        for (CCIndex f:cs.faces()){
            FaceData& fD = faceAttr[f];
            auto fv=faceVertices(cs,f);
                 for (uint j=0;j<fv.size();j++){
                     if(indexAttr[fv[j]].selected){
                         fv.erase(fv.begin()+j);
                         for(CCIndex k:fv){
                             if(indexAttr[k].selected){
                                 fD.cellZone = CellZone;
                                 break;
                             }
                         }
                     }
                 }
        }
        }

        if(SetPosition){
        for (CCIndex f:cs.faces()){
            FaceData& fD = faceAttr[f];
            auto fv=faceVertices(cs,f);
                 for (uint j=0;j<fv.size();j++){
                     if(indexAttr[fv[j]].selected){
                         fv.erase(fv.begin()+j);
                         for(CCIndex k:fv){
                             if(indexAttr[k].selected){
                                 fD.cellPosition = CellPosition;
                                 break;
                             }
                         }
                     }
                 }
        }
        }

        if(SetInitialAuxin){
        for (CCIndex f:cs.faces()){
            FaceData& fD = faceAttr[f];
            auto fv=faceVertices(cs,f);
                 for (uint j=0;j<fv.size();j++){
                     if(indexAttr[fv[j]].selected){
                         fv.erase(fv.begin()+j);
                         for(CCIndex k:fv){
                             if(indexAttr[k].selected){
                                 fD.initialAuxin = InitialAuxin;
                                 break;
                             }
                         }
                     }
                 }
        }
        }

        for(CCIndex f:cs.faces()){
            FaceData& fD = faceAttr[f];
            cellside[f] = fD.cellSide;
            celllayer[f] = fD.cellLayer;
            cellzone[f] = fD.cellZone;
            cellposition[f] = fD.cellPosition;
            initialauxin[f] = fD.initialAuxin;
        }

        mesh->updateAll();

        return true;
  }

  bool SetCellType::rewind(QWidget *parent)
  {
    // To rewind, we'll reload the mesh
    mesh = currentMesh();
    if(!mesh or mesh->file().isEmpty())
      throw(QString("General::rewind No current mesh, cannot rewind"));
    MeshLoad meshLoad(*this);
    meshLoad.setParm("File Name", mesh->file());
    return meshLoad.run();
  }

  ////---------------------------------CELL ATTR----------------------------------------//
  bool SetAttr::step() {
      mesh = getMesh("Mesh 1");
      if(!mesh or mesh->file().isEmpty())
          throw(QString("SetAttr::step No current mesh, cannot rewind"));

      CCStructure& cs = mesh->ccStructure("Tissue");
      CCIndexDataAttr& indexAttr = mesh->indexAttr();
      VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");

      Label =stringToBool(parm("Label"));
      if(Label){
          for(CCIndex f:cs.faces()) {
              CCIndexData& iD = indexAttr[f];
              iD.label=LABEL;
              LABEL++;
          }
      }

       StiffRef =stringToBool(parm("StiffRef"));
       if(StiffRef){
           for (CCIndex v:cs.vertices()){
              if(indexAttr[v].selected)
                   verticeAttr[v].stiffRefVertex=true;
              else
                  verticeAttr[v].stiffRefVertex=false;

           }
       }

       Bottom = stringToBool(parm("Bottom"));
       HighMass = stringToBool(parm("HighMass"));

        if(Bottom){
          for (CCIndex v:cs.vertices()){
              if(indexAttr[v].selected){
                 verticeAttr[v].bottom=true;
                 if(HighMass)
                    verticeAttr[v].invMass=1/BIG_VAL;
              }
              else{
                  verticeAttr[v].bottom=false;
                  verticeAttr[v].invMass=1;
              }
          }
        }

         CCIndexVec faces=cs.faces();
         i= parm("Index").toInt(); //NOW IS PER LABEL
         for(CCIndex f:cs.faces()) {
             mdxInfo<<indexAttr[f].label<<endl;
             if(indexAttr[f].label==i)
                indexAttr[f].selected=true;
         }

         Delete = stringToBool(parm("Delete"));
             if(Delete){
                 cs.deleteCell(faces[i]);
                 for (CCIndex e: cs.edges()){
                     CCIndexSet cb = cs.cobounds(e);
                     if(cb.size() == 0)
                         cs.deleteCell(e);
                 }
         }

         ZPOS= stringToBool(parm("ZPOS"));
         if (ZPOS){
          for (CCIndex v:cs.vertices()){
              indexAttr[v].pos[2]=0;
             }
         }

         EdgeLength= stringToBool(parm("EdgeLength"));
         if (EdgeLength){
          for (CCIndex e:cs.edges()){
              std::pair<CCIndex, CCIndex> endpoints=cs.edgeBounds(e);
              CCIndexData& v1D=indexAttr[endpoints.first];
              CCIndexData& v2D=indexAttr[endpoints.second];
              indexAttr[e].measure=sqrt(pow(v2D.pos[0]-v1D.pos[0],2)+pow(v2D.pos[1]-v1D.pos[1],2));
             }
         }

         mesh->updateAll("Tissue");

         return true;
  }



  ////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////TRIANGULATION////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////

  bool MyTriangulate::initialize(QWidget* parent){

        mesh = currentMesh();
        if(!mesh)
          throw(QString("Triangulate::initialize No current mesh"));

        if(!getProcess(parm("Subdivide Process"), subdivideProcess))
          throw(QString("Triangulate::initialize Cannot make Subdivide Process"));
        subdivideProcess->initialize(parent);

        maxArea=parm("Max Area").toDouble();

    return true;
  }

  bool MyTriangulate::step(bool firstTriangulation) {

      CCStructure& cs = mesh->ccStructure("Tissue");
      CCStructure& tt = mesh->ccStructure("Triangulate");

      CCIndexDataAttr& indexAttr = mesh->indexAttr();
      FaceDataAttr& faceAttr = mesh->attributes().attrMap<CCIndex, FaceData>("FaceData");
      EdgeDataAttr& edgeAttr = mesh->attributes().attrMap<CCIndex, EdgeData>("EdgeData");
      VerticeDataAttr& verticeAttr = mesh->attributes().attrMap<CCIndex, VerticeData>("VerticeData");

  //if triangulates the mesh for first time select all the faces and edges for triangulation
      if(firstTriangulation){
          for (CCIndex e:cs.edges())
              indexAttr[e].selected=true;
          for (CCIndex f:cs.faces())
              indexAttr[f].selected=true;
       }

  //First split edges
       std::vector<CCIndexVec> empty;
       double empty2=0;
       std::vector<edgeStruct>  empty3;
       double maxLength = sqrt(2 * maxArea);
       do {
           CCIndexVec edges;
           for(CCIndex e : cs.edges()) {
           if(indexAttr[e].selected){ //if divided seleted and then split
               auto eb = cs.edgeBounds(e);
               if(norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos) > maxLength)
                    edges.push_back(e);
            }
            }
            if(edges.size() > 0){
                for(CCIndex e : edges) {
                    EdgeData &eD = edgeAttr[e];
                    auto eb = cs.edgeBounds(e);
                    std::vector<edgeParData> edatas;
                    for(CCIndex f:cs.incidentCells(e,2)){
                        edgeParData edata;
                        CCIndexData &iD=indexAttr[f];
                        edata.label=iD.label;
                        edata.pin=eD.pin1[iD.label];
                        edatas.push_back(edata);
                    }

                    SplitStruct ss(e);
                    CCIndexFactory.fillSplitStruct(ss);
                    if(!cs.splitCell(ss)) {
                        mdxInfo << "SplitEdges::run Problem splitting cell" << endl;
                        continue;
                    }
                    // Set the position
                    indexAttr[ss.membrane].pos = (indexAttr[eb.first].pos + indexAttr[eb.second].pos)/2.0;
                    // We'll only update the membrane data
                    subdivideProcess->step(1, cs, ss, eb.first, eb.second, 0.5,true, empty,edatas,empty2,empty2,empty2,empty3);
                }
            }
            else
                break;
        } while(true);

        bool boundary = false;
        bool inOrder = true;
        bool center = true;

        int dim = cs.maxDimension;
        MeshBuilder mb(indexAttr, cs.maxDimension);

  // triangulate using triangle
        for(CCIndex f : cs.faces()) {
            auto& fIdx = indexAttr[f];
            int label = indexAttr[f].label;

            if(fIdx.selected==false){
                std::vector<Point3dVec> triPos;
                for (CCIndex ft:tt.faces()){
                if (indexAttr[ft].label==label){
                    CCIndexVec fv=faceVertices(tt,ft);
                    Point3dVec points;
                    for (CCIndex v:fv){
                        points.push_back(indexAttr[v].pos);
                    }
                    triPos.push_back(points);
                }
                }

                if(dim == 2) {
                    for(auto& tri:triPos){
                        mb.addFace({tri[0], tri[1], tri[2]}, label);
                    }
                }

        }
            else{
            CCIndexVec fVertices = faceVertices(cs, f);
            Point3dVec poly(fVertices.size());
            for(uint i = 0; i < fVertices.size(); i++)
                poly[i] = indexAttr[fVertices[i]].pos;
            Point3d polyNormal;
            std::vector<Point3i> triList;
            Point3dVec ptList;

            // calculate the polygon plane using pca
            Point3d planePos, planeNrml;
            Point3dVec polyProj = poly;
            findPolygonPlane(polyProj, planePos, planeNrml);
            if(fIdx.nrml * planeNrml > 0)
                planeNrml = -planeNrml;

            // project all points onto the plane
            projectPointsOnPlane(polyProj, planePos, planeNrml);

            // triangulate the polygon
            if(!triangulatePolygon3D(
                   maxArea, planeNrml, polyProj, triList, ptList, boundary, inOrder, center))
                mdxInfo << QString("%1: Error calling triangulatePolygon3D").arg(name()) << endl;

            // move the border points back to their original position
            if(ptList.size() < poly.size()) {
                mdxInfo << QString("%1: Error in ptList: %2 poly: %3 triList %4")
                               .arg(name())
                               .arg(ptList.size())
                               .arg(triList.size())
                               .arg(poly.size())
                        << endl;
                continue;
            }
            for(uint i = 0; i < poly.size(); i++)
                ptList[i] = poly[i];

            // create the triangular faces
            if(dim == 2) {
                for(const Point3i& tri : triList)
                    mb.addFace({ptList[tri[0]], ptList[tri[1]], ptList[tri[2]]}, label);
            }
        }
      }


      //Pass velocity, for new vertices in the border (cs) is 1/2 (subdivide)
      std::vector<CCIndex> allIndex;
      std::vector<Point3d> allPos;
      for(CCIndex v:tt.vertices()) //Save pos and index to restore
          allIndex.push_back(v);

      for(CCIndex v:cs.vertices())
          allIndex.push_back(v);

      sort(allIndex.begin(), allIndex.end() );
      allIndex.erase( unique(allIndex.begin(), allIndex.end() ), allIndex.end());

      for(CCIndex v:allIndex)
          allPos.push_back(indexAttr[v].pos);

      //Take the position to find the new stiffRefVertex later

      //Delete Previous
      for (CCIndex f: tt.faces())
           tt.deleteCell(f);
      for (CCIndex e: tt.edges())
           tt.deleteCell(e);
      for (CCIndex v: tt.vertices())
           tt.deleteCell(v);

      QString ccNameOut = "Triangulate";          //Create a CCStructure
      CCStructure& csOut = mesh->ccStructure(ccNameOut);
      csOut = CCStructure(dim);
      mb.createMesh(csOut);

      for(CCIndex e:csOut.edges()){
          auto edgeFaces=csOut.incidentCells(e,2);
          if(edgeFaces.size()>2){
              for(CCIndex f:edgeFaces)
                  indexAttr[f].selected = true;
              throw(QString("WRONG"));
          }
      }

      //UPDATE TISSUE WITH SAME CCINDEX AS IN TRIANGULATED TO GROWTH BOTH AT THE SAME TIME
  //Pass Vertice Attr
      for(CCIndex v:csOut.vertices()){
          std::pair<bool, int> vertex = findInVector<Point3d>(allPos,indexAttr[v].pos);
          if(vertex.first){
             verticeAttr[v].velocity=verticeAttr[allIndex[vertex.second]].velocity; //pass velocity
             verticeAttr[v].invMass=verticeAttr[allIndex[vertex.second]].invMass; //pass invMass
             verticeAttr[v].bottom=verticeAttr[allIndex[vertex.second]].bottom; //pass bottom
             verticeAttr[v].prevPos=verticeAttr[allIndex[vertex.second]].prevPos; //pass prevPos for updating lastPos
             verticeAttr[v].stiffRefVertex=verticeAttr[allIndex[vertex.second]].stiffRefVertex; //stiff ref vertex as attr for the cluster simulation
          }
      }

  //Save EdgeAttr
      struct edgeTriStruct
        {
      double restLength;
      double intercellularAuxin;
      double posPin1;
      double negPin1;
      double posSensitivity;
      double negSensitivity;
      double posPrevExport;
      double negPrevExport;
      Point3d fBound; //position
      Point3d sBound;
         };
      std::vector<edgeTriStruct> initTriEdges;
      edgeTriStruct initTriEdge;

      for (CCIndex e:cs.edges()){ //link edges from tissue and triangulate
          EdgeData& eD=edgeAttr[e];
          CCIndex es=eD.sameEdge;
          auto eb=cs.edgeBounds(e);
          Point3d p1=indexAttr[eb.first].pos;
          Point3d p2=indexAttr[eb.second].pos;
          initTriEdge.restLength=edgeAttr[es].restLength; //For restLength I use triangulate because is stored directly here
          initTriEdge.intercellularAuxin=eD.intercellularAuxin;
          initTriEdge.fBound=p1;
          initTriEdge.sBound=p2;
          CCIndex fpos, fneg;

          std::set<CCIndex> edgeFace =cs.incidentCells(e,2);
          if(edgeFace.size()==1){
             fpos=*edgeFace.begin();
          }
           else if(edgeFace.size()==2){
               std::vector<CCIndex> cb(edgeFace.begin(),edgeFace.end());
               if(indexAttr[cb[0]].measure>indexAttr[cb[1]].measure){
                  fpos=cb[0];
                  fneg=cb[1];
               }
               else if(indexAttr[cb[1]].measure>indexAttr[cb[0]].measure){
                   fpos=cb[1];
                   fneg=cb[0];
               } else
                   mdxInfo<<"WRONG ORIENTATION TRIANGULATION"<<endl;

          }
          initTriEdge.posPin1=eD.pin1[indexAttr[fpos].label];
          initTriEdge.posSensitivity=eD.pin1Sensitivity[indexAttr[fpos].label];
          initTriEdge.posPrevExport=eD.edge_export_prev[indexAttr[fpos].label];
          initTriEdge.negPin1=eD.pin1[indexAttr[fneg].label];
          initTriEdge.negSensitivity=eD.pin1Sensitivity[indexAttr[fneg].label];
          initTriEdge.negPrevExport=eD.edge_export_prev[indexAttr[fneg].label];
          initTriEdges.push_back(initTriEdge);
      }


  //Save FaceAttr
       struct faceStruct
             {
               std::vector<CCIndex> borderVertex;
               int label;
               double measure;
               int cellside;
               int celllayer;
               int cellzone;
               int cellposition;
               int initialauxin;
               bool oc;
               Point3d oldcentroid;
               Point3d force;
               Point3d stress;
               double a;
               Point3d a1;
               Point3d a2;
               double prevarea;
               double auxin;
               Point3d fluxVector;
               double pin1;
             };

        std::vector<faceStruct> initFaces;
        faceStruct initFace;

       for (CCIndex f:cs.faces()){
           CCIndexData& ifD=indexAttr[f];
           initFace.borderVertex=faceVertices(cs,f); //counterclockwise order
           initFace.label=ifD.label;
           initFace.measure=ifD.measure;
           initFace.cellside=faceAttr[f].cellSide;
           initFace.celllayer=faceAttr[f].cellLayer;
           initFace.cellzone=faceAttr[f].cellZone;
           initFace.cellposition=faceAttr[f].cellPosition;
           initFace.initialauxin=faceAttr[f].initialAuxin;
           initFace.oldcentroid=faceAttr[f].oldcentroid;
           initFace.force=faceAttr[f].force;
           initFace.stress=faceAttr[f].stress;
           initFace.a=faceAttr[f].a;
           initFace.a1=faceAttr[f].a1;
           initFace.a2=faceAttr[f].a2;
           initFace.prevarea=faceAttr[f].prevArea;
           initFace.auxin=faceAttr[f].auxin;
           initFace.fluxVector=faceAttr[f].auxinFluxVector;
           initFace.pin1=faceAttr[f].pin1;
           initFaces.push_back(initFace);
        }
        //If vertex position in Tissue==Triangulate, take CCIndex from Triangulate to Tissue
        for (uint i = 0; i < initFaces.size(); i++){
             for (uint j = 0; j < initFaces[i].borderVertex.size(); j++){
                 CCIndex v1=initFaces[i].borderVertex[j];
                 CCIndexData& v1Idx = indexAttr[v1];
                 for(CCIndex v2 : csOut.vertices()) {
                     CCIndexData& v2Idx = indexAttr[v2];
                     if (v1Idx.pos==v2Idx.pos) {
                        initFaces[i].borderVertex[j]=v2;

                        if(v1==refVertex) //REFVERTEX FOR VISUALIZE
                            refVertex=v2;

                     }
                 }
             }
        }

  //Delete Original Tissue
        for (CCIndex f: cs.faces())
             cs.deleteCell(f);
        for (CCIndex e: cs.edges())
             cs.deleteCell(e);
        for (CCIndex v: cs.vertices())
             cs.deleteCell(v);
  //Add New Faces (and Edges)
           for (uint i = 0; i < initFaces.size(); i++){
           CCIndex ft = CCIndexFactory.getIndex();
           addFace(cs, ft, initFaces[i].borderVertex);
           indexAttr[ft].label = initFaces[i].label;
           indexAttr[ft].measure= initFaces[i].measure;
           faceAttr[ft].cellSide =initFaces[i].cellside;
           faceAttr[ft].cellLayer =initFaces[i].celllayer;
           faceAttr[ft].cellZone =initFaces[i].cellzone;
           faceAttr[ft].cellPosition =initFaces[i].cellposition;
           faceAttr[ft].initialAuxin =initFaces[i].initialauxin;
           faceAttr[ft].oldcentroid =initFaces[i].oldcentroid;
           faceAttr[ft].force =initFaces[i].force;
           faceAttr[ft].stress =initFaces[i].stress;
           faceAttr[ft].a1 =initFaces[i].a1;
           faceAttr[ft].a2 =initFaces[i].a2;
           faceAttr[ft].prevArea =initFaces[i].prevarea;
           faceAttr[ft].auxin =initFaces[i].auxin;
           faceAttr[ft].auxinFluxVector =initFaces[i].fluxVector;
           faceAttr[ft].pin1 =initFaces[i].pin1;
           }

  //Link faces from triangulate to tissue and Pass faceAttr
        for (CCIndex fs:cs.faces()){
            for (CCIndex ft:csOut.faces()){
                if(indexAttr[ft].label==indexAttr[fs].label){
                  faceAttr[fs].cellFaces.push_back(ft);

                  faceAttr[ft].cellSide=faceAttr[fs].cellSide;
                  faceAttr[ft].cellLayer=faceAttr[fs].cellLayer;
                  faceAttr[ft].cellZone=faceAttr[fs].cellZone;
                  faceAttr[ft].cellPosition=faceAttr[fs].cellPosition;
                }
            }
        }


  //Update edges length (measure)
        for(CCIndex e:csOut.edges()){
            auto eb=csOut.edgeBounds(e);
            indexAttr[e].measure=norm(indexAttr[eb.first].pos-indexAttr[eb.second].pos);
        }

        for(CCIndex e:cs.edges()){ //Update CS, if not we lose the info in the division
           auto eb=cs.edgeBounds(e);
           indexAttr[e].measure=norm(indexAttr[eb.first].pos-indexAttr[eb.second].pos);
        }

  //Assign border and wall(boundary) edge again to triangulate (the ones that grow)
        for(CCIndex e:csOut.edges()) {
           EdgeData& eD = edgeAttr[e];
           CCIndexSet cb = csOut.cobounds(e); //adjacent faces
           std::vector<CCIndex> cbV(cb.begin(), cb.end());
             if(cb.size() == 1) { // edge in the border
               eD.border = true;
               eD.f = *cb.begin(); //le asocia la face al edge del borde
               eD.g = (CCIndex)0; //le asocia nada?
               eD.wall = true;
             } else { //edge not in the border
                 eD.border = false;
                 if (indexAttr[cbV[0]].label != indexAttr[cbV[1]].label)
                     eD.wall = true;
                 else
                     eD.wall = false;
             }
        }

  //Pass EdgeAttr
         //restLLength (directly to triangulate)
        for (CCIndex e:csOut.edges()){
            EdgeData& eD = edgeAttr[e];
            auto eb=csOut.edgeBounds(e);
            if(eD.wall){
                Point3d p1=indexAttr[eb.first].pos;
                Point3d p2=indexAttr[eb.second].pos;
            for (uint k = 0; k < initTriEdges.size(); k++){
            if((p1==initTriEdges[k].fBound && p2==initTriEdges[k].sBound) || (p2==initTriEdges[k].fBound && p1==initTriEdges[k].sBound)){
                eD.restLength=initTriEdges[k].restLength;
                break;
            }
            else
                eD.restLength=0; //for the new walls after dividing a cell
            }
            }
            else
                eD.restLength=0; //for the inner edges of triangulate
        }
        //chemicalAttr
        for (CCIndex e:cs.edges()){
            EdgeData& eD = edgeAttr[e];
            CCIndex fpos, fneg;
            std::set<CCIndex> edgeFace =cs.incidentCells(e,2);
            if(edgeFace.size()==1){
               fpos=*edgeFace.begin();
            }
            else if(edgeFace.size()==2){
               std::vector<CCIndex> cb(edgeFace.begin(),edgeFace.end());
               if(indexAttr[cb[0]].measure>indexAttr[cb[1]].measure){
                  fpos=cb[0];
                  fneg=cb[1];
               }
               else if(indexAttr[cb[1]].measure>indexAttr[cb[0]].measure){
                   fpos=cb[1];
                   fneg=cb[0];
               } else
                   mdxInfo<<"WRONG ORIENTATION TRIANGULATION"<<endl;

          }

          auto eb=cs.edgeBounds(e);
          Point3d p1=indexAttr[eb.first].pos;
          Point3d p2=indexAttr[eb.second].pos;
          for (uint k = 0; k < initTriEdges.size(); k++){
          if((p1==initTriEdges[k].fBound && p2==initTriEdges[k].sBound) || (p2==initTriEdges[k].fBound && p1==initTriEdges[k].sBound)){
              eD.intercellularAuxin=initTriEdges[k].intercellularAuxin;
              eD.pin1[indexAttr[fpos].label]=initTriEdges[k].posPin1;
              eD.pin1Sensitivity[indexAttr[fpos].label]=initTriEdges[k].posSensitivity;
              eD.edge_export_prev[indexAttr[fpos].label]=initTriEdges[k].posPrevExport;
              eD.pin1[indexAttr[fneg].label]=initTriEdges[k].negPin1;
              eD.pin1Sensitivity[indexAttr[fneg].label]=initTriEdges[k].negSensitivity;
              eD.edge_export_prev[indexAttr[fneg].label]=initTriEdges[k].negPrevExport;
              break; //if found we can pass to the next edge;
          }
      }
      }

  //Link edges from both tisues
      for (CCIndex e1:csOut.edges()){
          auto eb1=csOut.edgeBounds(e1);
          if(edgeAttr[e1].wall){
              Point3d p11=indexAttr[eb1.first].pos;
              Point3d p12=indexAttr[eb1.second].pos;
          for (CCIndex e2:cs.edges()){
              auto eb2=cs.edgeBounds(e2);
              Point3d p21=indexAttr[eb2.first].pos;
              Point3d p22=indexAttr[eb2.second].pos;
                  if((p11==p21 && p12==p22) || (p12==p21 && p11==p22)){
                      edgeAttr[e1].sameEdge=e2;
                      edgeAttr[e2].sameEdge=e1;
                      break;
                  }
          }
          }
      }

      if(firstTriangulation)
          mesh->updateAll();

        return true;
  }







 REGISTER_PROCESS(General);

 REGISTER_PROCESS(Chemics);
 REGISTER_PROCESS(ChemicalSolver);

 REGISTER_PROCESS(Mechanics);
 REGISTER_PROCESS(MechanicalSolver);
 REGISTER_PROCESS(restUpdate);

 REGISTER_PROCESS(MySubdivide);
 REGISTER_PROCESS(MyDivide);
 REGISTER_PROCESS(MyTriangulate);

 REGISTER_PROCESS(SetCellType);
 REGISTER_PROCESS(SetAttr);




}
