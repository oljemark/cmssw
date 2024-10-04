#include <cstdio>
#include <string>
#include <sstream>
#include <TXMLEngine.h>
#include <iostream>

const char* FINISHED_PATH =
    "data/324579.0-25/103,104,105,123,124,125/ev=1E4,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/"
    "finished";
const char* FILE_PATH =
    "data/324579.0-25/103,104,105,123,124,125/ev=1E4,pl=3,units=2,ovlp=f,3rpsInO=t/s+sr,std,1rotzIt=0,eMuMvRot=t/"
    "iteration5/diagnostics.root";

std::unique_ptr<TFile> OpenRootFile(const char* filePath) {
  auto file = std::unique_ptr<TFile>(TFile::Open(filePath));
  if (!file || file->IsZombie()) {
    std::cerr << "Error: could not open the file." << std::endl;
    return nullptr;
  }
  return file;
}

TDirectory* GetDirectory(TDirectory* parentDir, const char* dirName) {
  TDirectory* dir = parentDir->GetDirectory(dirName);
  if (!dir) {
    std::cerr << "Error: could not find '" << dirName << "' directory." << std::endl;
  }
  return dir;
}

TH1D* GetHistogram(TDirectory* dir, const char* histName) {
  TH1D* hist = nullptr;
  dir->GetObject(histName, hist);
  if (!hist) {
    std::cerr << "Error: could not find the histogram '" << histName << "'." << std::endl;
  }
  return hist;
}

int GetResidualsHistogramMean(const char* filePath, int detectorID, double& mean) {
  auto file = OpenRootFile(FILE_PATH);
  if (!file)
    return 0;

  // Navigate to the directories
  TDirectory* commonDir = GetDirectory(file.get(), "common");
  if (!commonDir)
    return 0;

  TDirectory* residualsDir = GetDirectory(commonDir, "residuals");
  if (!residualsDir)
    return 0;

  TDirectory* dirDetector = GetDirectory(residualsDir, to_string(detectorID).c_str());
  if (!dirDetector)
    return 0;

  // Histogram's name
  std::stringstream ss;
  ss << detectorID << ": total_selected;1";

  // Retrieve the histogram
  TH1D* hist = GetHistogram(dirDetector, ss.str().c_str());
  if (!hist)
    return 0;

  // // Get the mean value from the histogram
  mean = hist->GetMean();
  std::cout << "Mean value of the histogram: " << mean << std::endl;

  // Close the ROOT file
  file->Close();

  return 1;
}

void ModifyXMLValues(const char* inputFile, const char* outputFile, const float offset) {
  // Create XML engine
  TXMLEngine* xml = new TXMLEngine;

  // Parse an XML file
  XMLDocPointer_t xmldoc = xml->ParseFile(inputFile);
  if (!xmldoc) {
    std::cerr << "Error: could not parse the XML file." << std::endl;
    delete xml;
    return;
  }

  // Access the main node of the XML document
  XMLNodePointer_t mainNode = xml->DocGetRootElement(xmldoc);

  // Navigate to the ConstantsSection node
  XMLNodePointer_t constantsSection = xml->GetChild(mainNode);
  while (constantsSection && std::string(xml->GetNodeName(constantsSection)) != "ConstantsSection") {
    constantsSection = xml->GetNext(constantsSection);
  }
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
      const char* original_raw = xml->GetAttr(constant, "value");
      float original;
      sscanf(original_raw, " +%f*mm", &original);

      std::stringstream ss;
      ss << "+" << original + offset << "*mm";
      std::string newValue = ss.str();

      xml->FreeAttr(constant, "value");
      xml->NewAttr(constant, nullptr, "value", newValue.c_str());
    }
    constant = xml->GetNext(constant);
  }

  // Save the modified XML file
  xml->SaveDoc(xmldoc, outputFile);

  // Free memory
  xml->FreeDoc(xmldoc);
  delete xml;
}

bool file_exists(const char* name) {
  struct stat buffer;
  return (stat(name, &buffer) == 0);
}

void wait_for_finish(const char* finish_indicator) {
  while (!file_exists(finish_indicator)) {
  }
}

// This is the main program
int script_2() {
  // Create a graph to store mean values
  TGraphErrors* graph = new TGraphErrors();

  const char* inputFile = "RP_Dist_Beam_Cent_original.xml";
  const char* outputFile = "RP_Dist_Beam_Cent.xml";

  for (double i = 0; i < 0.1; i = i + 0.1) {
    std::cout << "Iteration " << i << std::endl;

    // Change the values of a test files
    ModifyXMLValues(inputFile, outputFile, i);

    remove(FINISHED_PATH);

    // Run tests
    system("./run");

    wait_for_finish(FINISHED_PATH);

    double mean;
    // Read data from test results
    auto result = GetResidualsHistogramMean(FILE_PATH, 1998585856, mean);
    if (!result)
      break;

    graph->AddPoint(i, mean);
  }

  // Create a canvas and draw the graph
  TCanvas* canvas = new TCanvas("canvas", "Mean residuals of detector vs RP offset", 1440, 650);
  graph->SetTitle("Mean residuals of detector vs RP offset");
  graph->GetXaxis()->SetTitle("RP offset [mm]");
  graph->GetYaxis()->SetTitle("Residuals of 1998585856 detector");
  graph->SetMarkerStyle(20);
  graph->Draw("AC*");

  // Save the plot as an image (optional)
  canvas->SaveAs("residuals_for_1998585856_v2.png");

  return 0;
}