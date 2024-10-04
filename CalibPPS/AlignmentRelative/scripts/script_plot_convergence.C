#include <cstdio>
#include <vector>
#include <string>
#include <sstream>
#include <TXMLEngine.h>
#include <iostream>

const char* DATA_PATH =
    "data/324579.0-25/103,104,105,123,124,125/ev=1E4,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/";
const char* RESULT_FILE = "results_iteration_Jan.xml";

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
int script_plot_convergence() {
  std::stringstream ss;
  ss << DATA_PATH << "iteration1/" << RESULT_FILE;
  std::string result_file = ss.str();

  std::vector<int> sensorsID = GetSensors(result_file.c_str());

  if (sensorsID.empty()) {
    std::cout << "No sensor data could be extracted, abandoning the script" << std::endl;
    return 1;
  }

  std::map<int, float> ShiftsX, ShiftsY, RotsZ;

  // Creating graphs
  TMultiGraph* mg = new TMultiGraph();
  // gStyle->SetPalette(73);
  std::map<int, TGraphErrors*> graphs;

  for (const int& sensorID : sensorsID) {
    TGraphErrors* graph = new TGraphErrors();
    graph->SetMarkerStyle(8);
    graph->SetTitle(std::to_string(sensorID).c_str());
    graphs[sensorID] = graph;
    mg->Add(graph, "LP");
  }

  // Main loop
  for (int iteration = 1; iteration < 6; iteration++) {
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
      graph->AddPoint(iteration, std::abs(ShiftsX[sensorID]));
    }
  }

  // Create a canvas and draw the graph
  TCanvas* canvas = new TCanvas("canvas", "Mean residuals of detector vs RP offset", 1000, 650);
  canvas->SetLogy();
  // graph->SetTitle("Mean residuals of detector vs RP offset");
  mg->GetXaxis()->SetTitle("Iteration");
  mg->GetYaxis()->SetTitle("abs(Shift X) [um]");
  mg->Draw("A pmc plc");

  // Save the plot as an image (optional)
  TLegend* legend = canvas->BuildLegend();
  legend->SetNColumns(4);
  // canvas->SaveAs("shift_X_for_all_legend.png");
  canvas->Print("shift_X_for_all_legend.pdf");

  return 0;
}