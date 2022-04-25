

#ifndef myHistogram_h
#define myHistogram_h 1

#include <string>
#include <fstream>

// Written by Grant Berland
//
// Header-only histogramming class to record energy deposition
// per altitude bin


class myHistogram
{
public:

  myHistogram();
  myHistogram(double, double, int);
  
  ~myHistogram();

  // Tells program whether or not to write to histogram if other 
  // data collections methods are selected
  void InitializeHistogram(){ initializedFlag = 1; };
  
  // Method to fill histogram with data
  void AddCountToBin(unsigned int, double);


  void AddCountTo2DHistogram(unsigned int, double);

  // Writes vector to file name provided
  void WriteHistogramToFile(std::string);

  // Writes arbitrary data line directly to a filename specified
  void WriteDirectlyToFile(std::string, double*, size_t);


  void Write2DHistogram(std::string);
  
  // Always base 10
  void GenerateLogspaceBins(double, double, int, double[]);

private:
  // Need to run InitializeHistogram() method to write particle energy 
  // depostion to histogram
  int initializedFlag = 0;   
  
  // Array initialized to zeros (with fixed resolution)
  double histogramArray[1000] = {};

  double binEdges[101] = {};

  double twoDhistogramArray[1000][100] = {};
};

// Inline constructor and destructor methods
inline myHistogram::myHistogram()
  : histogramArray()
{}

inline myHistogram::myHistogram(double binStart, double binStop, int nBins) 
  : binEdges(),
    twoDhistogramArray()
{

  GenerateLogspaceBins(binStart, binStop, nBins, binEdges);

}


inline myHistogram::~myHistogram()
{}

inline void myHistogram::AddCountToBin(unsigned int binAddress, 
			               double       energyToAdd)
{
  histogramArray[binAddress] += energyToAdd;
}

inline void myHistogram::WriteHistogramToFile(std::string filename)
{
  //if(initializedFlag == 1)
  //{
    std::ofstream outputFile;
  
    // OPEN
    outputFile.open(filename, std::ios_base::app);

    for(unsigned int i=0; i<1000; i++)
    {
      outputFile << histogramArray[i] << "\n";
    }
  
    // CLOSE
    outputFile.close();
    //}
  
  /*
  else if(initializedFlag == 0)
  {
    std::cout << "Histogram not initialized" << std::endl;
  }
  */
}
  
inline void myHistogram::WriteDirectlyToFile(std::string filename, 
		                             double *data,
					     size_t  arraySize)
{
  std::ofstream outputFile;
  
  // OPEN
  outputFile.open(filename, std::ios_base::app);

  // Writes entire line in 'data[]' to file, comma delimited
  for(unsigned int i=0; i<arraySize-1; i++)
  {
    outputFile << data[i] << ',';
  }
  
  // Ends output with last entry and newline character
  outputFile << data[arraySize-1] << '\n';

  // CLOSE
  outputFile.close();
}



inline void myHistogram::AddCountTo2DHistogram(unsigned int address1, double value)
{

  if (value > 0.)
  {
    for(unsigned int i=0; i<101; i++)
    {   
      if(binEdges[i] > value)
      { 
	twoDhistogramArray[address1][i-1] += 1;
	break; 
      }
    }

   }

}

inline void myHistogram::Write2DHistogram(std::string filename)
{
    std::ofstream outputFile;
  
    // OPEN
    outputFile.open(filename, std::ios_base::app);

    for(unsigned int i=0; i<1000; i++)
    {
	for (unsigned int j=0; j<100-1; j++)
	{
	    outputFile << twoDhistogramArray[i][j] << ",";
	}
	outputFile << twoDhistogramArray[i][100-1] << "\n";
    }
  
    // CLOSE
    outputFile.close();
 
}

inline void myHistogram::GenerateLogspaceBins(G4double start, 
					      G4double end, 
					      G4int nBins, 
					      G4double binArray[])
{

  // step size 
  double c = (end - start) / (nBins - 1);
  
  // fill vector 
  for (int i = 0; i < nBins-1; ++i)
  {
    binArray[i] = std::pow(10., start + i * c);
  }

  // fix last entry to 10^b 
  binArray[nBins-1] = std::pow(10., end);

}

#endif
