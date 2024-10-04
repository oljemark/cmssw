#include <cstdio>
#include <vector>
#include <string>
#include <sstream>
#include <TXMLEngine.h>
#include <iostream>

const char* FINISHED_PATH =
    "data/369591.0-200/103,104,105,123,124,125/ev=1E4,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/"
    "finished";
const char* DATA_PATH =
    "data/369591.0-200/103,104,105,123,124,125/ev=1E4,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/";
// const char* RESULT_FILE = "results_cumulative_factored_Jan.xml";
const char* RESULT_FILE = "results_iteration_Jan.xml";

int getRomanPotNumber(unsigned long detectorId) {
  // Mapping of detector ids to their corresponding Roman Pot numbers
  std::map<unsigned long, int> detectorToRomanPot = {
      {1998585856, 104}, {1998618624, 104}, {1998651392, 104}, {1998684160, 104}, {1998716928, 104}, {1998749696, 104},
      {1998782464, 104}, {1998815232, 104}, {1998848000, 104}, {1998880768, 104}, {1999110144, 105}, {1999142912, 105},
      {1999175680, 105}, {1999208448, 105}, {1999241216, 105}, {1999273984, 105}, {1999306752, 105}, {1999339520, 105},
      {1999372288, 105}, {1999405056, 105}, {2006974464, 124}, {2007007232, 124}, {2007040000, 124}, {2007072768, 124},
      {2007105536, 124}, {2007138304, 124}, {2007171072, 124}, {2007203840, 124}, {2007236608, 124}, {2007269376, 124},
      {2007498752, 125}, {2007531520, 125}, {2007564288, 125}, {2007597056, 125}, {2007629824, 125}, {2007662592, 125},
      {2007695360, 125}, {2007728128, 125}, {2007760896, 125}, {2007793664, 125}, {2031616000, 103}, {2031681536, 103},
      {2031747072, 103}, {2031812608, 103}, {2031878144, 103}, {2031943680, 103}, {2040004608, 123}, {2040070144, 123},
      {2040135680, 123}, {2040201216, 123}, {2040266752, 123}, {2040332288, 123}};

  // Find the detector id in the map and return the corresponding Roman Pot number
  auto it = detectorToRomanPot.find(detectorId);
  if (it != detectorToRomanPot.end()) {
    return it->second;
  } else {
    // If the detector id is not found, return -1 or any other appropriate value to indicate an error
    return -1;
  }
}

bool file_exists(const char* name) {
  struct stat buffer;
  return (stat(name, &buffer) == 0);
}

void wait_for_finish(const char* finish_indicator) {
  while (!file_exists(finish_indicator)) {
  }
}

XMLNodePointer_t AccessXMLSection(TXMLEngine* xml,
                                  const char* inputFile,
                                  const char* sectionName,
                                  XMLDocPointer_t* xmldoc) {
  std::cout << "AccessXMLSection:inputFile=" << inputFile << std::endl;
  if (!file_exists(inputFile))
    return nullptr;
  // Parse an XML file
  *xmldoc = xml->ParseFile(inputFile);
  if (!*xmldoc) {
    std::cerr << "Error: could not parse the XML file." << std::endl;
    delete xml;
    return nullptr;
  }

  // Access the main node of the XML document
  XMLNodePointer_t mainNode = xml->DocGetRootElement(*xmldoc);

  // Navigate to the ConstantsSection node
  XMLNodePointer_t selectedSection = xml->GetChild(mainNode);
  while (selectedSection && std::string(xml->GetNodeName(selectedSection)) != sectionName) {
    selectedSection = xml->GetNext(selectedSection);
  }

  if (!selectedSection) {
    std::cerr << "Error: could not find " << sectionName << " node." << std::endl;
    xml->FreeDoc(*xmldoc);
    delete xml;
    return nullptr;
  }

  return selectedSection;
}

std::vector<int> GetSensors(const char* inputFile) {
  std::vector<int> sensors;

  TXMLEngine* xml = new TXMLEngine;
  XMLDocPointer_t xmldocInstance;
  XMLDocPointer_t* xmldoc = &xmldocInstance;

  auto dataSection = AccessXMLSection(xml, inputFile, "iov", xmldoc);

  if (dataSection == nullptr) {
    return sensors;
  }

  XMLNodePointer_t child = xml->GetChild(dataSection);

  while (child) {
    int sensorID = xml->GetIntAttr(child, "id");
    sensors.push_back(sensorID);
    child = xml->GetNext(child);
  }

  // Free memory
  xml->FreeDoc(*xmldoc);
  delete xml;

  return sensors;
}

void ModifyXMLValues(const char* inputFile, const char* outputFile, const float offset) {
  TXMLEngine* xml = new TXMLEngine;
  XMLDocPointer_t xmldocInstance;
  XMLDocPointer_t* xmldoc = &xmldocInstance;

  auto constantsSection = AccessXMLSection(xml, inputFile, "ConstantsSection", xmldoc);

  if (!constantsSection) {
    std::cerr << "Error: could not find 'ConstantsSection' node." << std::endl;
    xml->FreeDoc(xmldoc);
    delete xml;
    return;
  }

  // Traverse the ConstantsSection and modify all Constant values
  XMLNodePointer_t constant = xml->GetChild(constantsSection);
  while (constant) {
    if (std::string(xml->GetNodeName(constant)) == "Constant") {
      // Read the old value
      const char* raw_name = xml->GetAttr(constant, "name");
      std::string name(raw_name);

      if (name == "RP_220_Right_Det_Dist_3") {
        const char* original_raw = xml->GetAttr(constant, "value");
        float original;
        sscanf(original_raw, " %f*mm", &original);

        std::stringstream ss;
        ss << "+" << original + offset << "*mm";
        std::string newValue = ss.str();

        xml->FreeAttr(constant, "value");
        xml->NewAttr(constant, nullptr, "value", newValue.c_str());
      }
      // if (name == "RP_147_Right_Det_Dist_5"){
      //     const char* original_raw = xml->GetAttr(constant, "value");
      //     float original;
      //     sscanf(original_raw, " %f*mm", &original);

      //     std::stringstream ss;
      //     ss << "+" << original + offset << "*mm";
      //     std::string newValue = ss.str();

      //     xml->FreeAttr(constant, "value");
      //     xml->NewAttr(constant, nullptr, "value", newValue.c_str());
      // }
    }
    constant = xml->GetNext(constant);
  }

  // Save the modified XML file
  xml->SaveDoc(*xmldoc, outputFile);

  // Free memory
  xml->FreeDoc(*xmldoc);
  delete xml;
}

int GetShiftData(const char* inputFile,
                 std::map<int, float>& ShiftsX,
                 std::map<int, float>& ShiftsY,
                 std::map<int, float>& RotsZ) {
  TXMLEngine* xml = new TXMLEngine;
  XMLDocPointer_t xmldocInstance;
  XMLDocPointer_t* xmldoc = &xmldocInstance;
  std::cout << "Accesing shift data" << std::endl;

  auto dataSection = AccessXMLSection(xml, inputFile, "iov", xmldoc);
  std::cout << "Accesed section" << std::endl;

  if (dataSection == nullptr) {
    return 1;
  }

  XMLNodePointer_t child = xml->GetChild(dataSection);

  while (child) {
    int sensorID = xml->GetIntAttr(child, "id");

    // if (xml->GetNodeName(child) != std::string("rp")) {
    //     child = xml->GetNext(child);
    //     continue;
    // };
    // Read sh_x
    const char* rawShiftX = xml->GetAttr(child, "sh_x");
    float shiftX;
    sscanf(rawShiftX, "%f", &shiftX);
    ShiftsX[sensorID] = shiftX;

    // Read sh_y
    const char* rawShiftY = xml->GetAttr(child, "sh_y");
    float shiftY;
    sscanf(rawShiftY, "%f", &shiftY);
    ShiftsY[sensorID] = shiftY;

    // Read rot_z
    const char* rawRotZ = xml->GetAttr(child, "rot_z");
    float rotZ;
    sscanf(rawRotZ, "%f", &rotZ);
    RotsZ[sensorID] = rotZ;

    child = xml->GetNext(child);
  }

  // Free memory
  xml->FreeDoc(*xmldoc);
  delete xml;

  return 0;
}

// This is the main program
int script_plot_convergence_of_plane() {
  // Create a graph to store mean values
  TGraphErrors* graph = new TGraphErrors();

  const char* inputFile = "RP_Dist_Beam_Cent_original.xml";
  const char* outputFile = "RP_Dist_Beam_Cent.xml";
  int sensorID = 2040004608;

  std::stringstream ss;
  ss << DATA_PATH << "iteration1/" << RESULT_FILE;
  std::string result_file = ss.str();
  std::map<int, float> ShiftsX, ShiftsY, RotsZ;

  for (double offset = -0.03; offset < 0.05; offset = offset + 0.01) {
    std::cout << "Iteration " << offset << std::endl;

    // Change the values of a test files
    ModifyXMLValues(inputFile, outputFile, offset);

    remove(FINISHED_PATH);

    // Run tests
    system("./run");

    wait_for_finish(FINISHED_PATH);
    std::cout << "Finished waiting" << std::endl;

    if (GetShiftData(result_file.c_str(), ShiftsX, ShiftsY, RotsZ) != 0) {
      std::cout << "Skipping offset " << offset << ", couldn't access the data" << std::endl;
      continue;
    }

    // Converting offset in mm to offset in um
    graph->AddPoint(offset * 1000, ShiftsY[sensorID]);
    std::cout << "Calculated offset is: " << ShiftsY[sensorID] << std::endl;
  }

  // Fit a line to the plot
  TF1 f("f1", "[0]+x*[1]", 0, 0.5);
  f.SetLineColor(kRed);
  graph->Fit("f1");

  // Create a canvas and draw the graph
  TCanvas* canvas = new TCanvas("canvas", "", 1440, 650);
  graph->SetTitle("Calculated Y offset of the RP vs added offset from the beam");
  graph->GetXaxis()->SetTitle("RP 123 added offset from the beam [um]");
  graph->GetYaxis()->SetTitle("Calculated Y offset [um] of RP 123");
  graph->SetMarkerStyle(20);
  gStyle->SetOptFit(1111);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.3);

  graph->Draw("AP");

  // Save the plot as an image (optional)
  canvas->SaveAs("offset_factored_y_for_RP123.png");

  return 0;
}