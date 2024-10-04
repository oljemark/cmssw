#include <cstdio>
#include <vector>
#include <string>
#include <sstream>
#include <TXMLEngine.h>
#include <iostream>

const char* DATA_PATH =
    "data/369591.0-200/103,104,105,123,124,125-excl2007236608,2007269376/ev=-1,pl=3,units=2,ovlp=f,3rpsInO=t/"
    "s+sr,std,1rotzIt=0,eMuMvRot=t/";
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

XMLNodePointer_t AccessXMLSection(TXMLEngine* xml,
                                  const char* inputFile,
                                  const char* sectionName,
                                  XMLDocPointer_t* xmldoc) {
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

int GetShiftData(const char* inputFile,
                 std::map<int, float>& ShiftsX,
                 std::map<int, float>& ShiftsY,
                 std::map<int, float>& RotsZ) {
  TXMLEngine* xml = new TXMLEngine;
  XMLDocPointer_t xmldocInstance;
  XMLDocPointer_t* xmldoc = &xmldocInstance;

  auto dataSection = AccessXMLSection(xml, inputFile, "iov", xmldoc);

  if (dataSection == nullptr) {
    return 1;
  }

  XMLNodePointer_t child = xml->GetChild(dataSection);

  while (child) {
    int sensorID = xml->GetIntAttr(child, "id");

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
int script_plot_convergence_per_pot() {
  std::stringstream ss;
  ss << DATA_PATH << "iteration1/" << RESULT_FILE;
  std::string result_file = ss.str();

  std::vector<int> sensorsID = GetSensors(result_file.c_str());

  if (sensorsID.empty()) {
    std::cout << "No sensor data could be extracted, abandoning the script" << std::endl;
    return 1;
  }

  std::map<int, float> ShiftsX, ShiftsY, RotsZ;

  // Create maps for Roman Pots
  std::map<int, std::vector<int>> romanPots;
  for (const int& sensorID : sensorsID) {
    int romanPotNumber = getRomanPotNumber(sensorID);
    if (romanPotNumber != -1) {
      romanPots[romanPotNumber].push_back(sensorID);
    }
  }

  // Create a canvas and divide it into subplots
  TCanvas* canvas = new TCanvas("canvas", "Mean residuals of detector vs RP offset", 1000, 1000);
  int nRomanPots = romanPots.size();
  int nCols = static_cast<int>(std::sqrt(nRomanPots));
  int nRows = (nRomanPots + nCols - 1) / nCols;  // ensure all RomanPots have a subplot
  canvas->Divide(nCols, nRows);

  std::map<int, TMultiGraph*> multiGraphs;
  std::map<int, TGraphErrors*> graphs;

  // Creating graphs for each Roman Pot
  for (const auto& rp : romanPots) {
    int rpNumber = rp.first;
    TMultiGraph* mg = new TMultiGraph("mg", ("RP " + std::to_string(rpNumber)).c_str());
    multiGraphs[rpNumber] = mg;

    for (const int& sensorID : rp.second) {
      TGraphErrors* graph = new TGraphErrors();
      graph->SetMarkerStyle(8);
      graph->SetTitle(std::to_string(sensorID).c_str());
      graphs[sensorID] = graph;
      mg->Add(graph, "LP");
    }
  }

  // Main loop
  for (int iteration = 1; iteration < 11; iteration++) {
    std::cout << "Iteration " << iteration << std::endl;

    std::stringstream ss;
    ss << DATA_PATH << "iteration" << iteration << "/" << RESULT_FILE;
    result_file = ss.str();
    std::cout << result_file << std::endl;

    if (GetShiftData(result_file.c_str(), ShiftsX, ShiftsY, RotsZ) != 0) {
      std::cout << "Skipping iteration " << iteration << ", couldn't access the data" << std::endl;
      continue;
    }

    for (const int& sensorID : sensorsID) {
      auto graph = graphs[sensorID];
      graph->AddPoint(iteration, std::abs(ShiftsY[sensorID]));
    }
  }

  // Draw each subplot
  int padIndex = 1;
  for (const auto& rp : romanPots) {
    int rpNumber = rp.first;

    auto pad = canvas->cd(padIndex);
    pad->SetLogy();

    TMultiGraph* mg = multiGraphs[rpNumber];
    mg->GetXaxis()->SetTitle("Iteration");
    mg->GetYaxis()->SetTitle("abs(Shift Y) [um]");
    mg->Draw("A pmc plc");

    TLegend* legend = pad->BuildLegend();
    legend->SetNColumns(2);
    padIndex++;
  }

  // Save the plot as an image (optional)
  canvas->Print("shift_Y_for_all_per_RP.pdf");

  return 0;
}