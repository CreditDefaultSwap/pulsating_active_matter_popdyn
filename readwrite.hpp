#ifndef READWRITE_HPP
#define READWRITE_HPP

#include <fstream>
#include <string>
#include <iostream>

//wp: to read and write plain text files 
#include <sstream>
#include <vector>

/////////////
// CLASSES //
/////////////

class Read;
class Write;


/*  READ
 *  ----
 *  Read from binary input file.
 */

class Read {

  public:

    // CONSTRUCTORS

    Read(std::string filename) :
      inputFile(filename),
      inputStream(inputFile.c_str(), std::ios::in | std::ios::binary),
      fileSize (0) {

      // get file size
      if ( inputStream ) {
        inputStream.seekg(0, std::ios_base::end);
        fileSize = tellg();
        inputStream.seekg(0, std::ios_base::beg);
      }
    }

    // DESTRUCTORS

    ~Read() { close(); }

    // METHODS

    // READ WITHOUT OFFSET
    template<class InClass> void read(InClass* in) {
      // Read from input file.

      if ( inputStream ) {
        if ( ! inputStream.read((char*) in, sizeof(InClass)) ) {
          std::cerr
            << "[READ FAILURE]"
              << " CURSOR: " << tellg()
              << " | FILE: " << inputFile
            << std::endl;
        }
      }
    }
    template<class InClass> InClass read() {
      // Read from input file.

      InClass in (0);
      read<InClass>(&in);
      return in;
    }

    // READ WITH OFFSET
    template<class InClass> void read(InClass* in,
      long int offset, std::ios_base::seekdir way = std::ios_base::beg) {
      // Read from itput file at a given position given by offset relative to
      // way.

      if ( inputStream ) {
        inputStream.seekg(offset, way); // set position in input stream
        read<InClass>(in); // read
      }
    }

    template<class InClass> InClass read(
      long int offset, std::ios_base::seekdir way = std::ios_base::beg) {
      // Read from itput file at a given position given by offset relative to
      // way.

      InClass in (0);
      read<InClass>(&in, offset, way);
      return in;
    }

    std::string getInputFile() const { return inputFile; } // returns input file name

    long int getFileSize() { return fileSize; } // returns size of file
    long int tellg() { return inputStream.tellg(); } // returns position in input stream

    bool is_open() { return inputStream.is_open(); } // file stream opened
    void close() { inputStream.close(); } // close file stream
    void open() { // (re-)open file stream
      if ( ! is_open() ) {
        inputStream.open(inputFile.c_str(),
          std::ios::in | std::ios::binary);
      }
    }

  private:

    // ATTRIBUTES

    std::string const inputFile; // input file name
    std::ifstream inputStream; // input binary file stream
    long int fileSize; // size of file

};


/*  WRITE
 *  -----
 *  Write to binary output file.
 */

class Write {

  public:

    // CONSTRUCTORS

    Write(std::string filename) :
      outputFile(filename),
      outputStream(outputFile.c_str(), std::ios::out | std::ios::binary) {}
    Write() : Write("") {} 

    // DESTRUCTORS

    ~Write() { close(); }

    // METHODS 
    // WRITE WITHOUT OFFSET
	//wp: write  directly to plain text
    template<class OutClass> void write_txt(OutClass out) {
      // Write to output file in plain text 
        outputStream << out << "\n";            
    }

    // WRITE WITHOUT OFFSET
    template<class OutClass> void write(OutClass out) {
      // Write to output file.

      if ( outputStream ) {//wp: this is some condensed way of writing to file (using .write) + automatically catching execution flag 
        if ( ! outputStream.write((char*) &out, sizeof(OutClass)) ) {           
			std::cerr
            << "[WRITE FAILURE]"
              << " OUTPUT: " << out
              << " | CURSOR: " << tellp()
              << " | FILE: " << outputFile
            << std::endl;
        };
      }
    }

    // WRITE WITH OFFSET
    template<class OutClass> void write(OutClass out,
      long int offset, std::ios_base::seekdir way = std::ios_base::beg) {
      // Write to output file at a given position given by offset relative to
      // way.

      if ( outputStream ) {
        outputStream.seekp(offset, way); // set position in output stream
        write<OutClass>(out); // write
        outputStream.seekp(0, std::ios_base::end); // set position to end of file
      }
    }

    std::string getOutputFile() const { return outputFile; } // returns output file name
    void setOutputFile(std::string filename) { outputFile=filename; } //wp: re-sets name if necessary 

    long int tellp() { return outputStream.tellp(); } // returns position in output stream

    bool is_open() { return outputStream.is_open(); } // file stream opened
    void flush() { outputStream.flush(); } // flush file stream
    void close() { outputStream.close(); } // close file stream
    void open() { // (re-)open file stream
      if ( ! is_open() ) {
        outputStream.open(outputFile.c_str(),
          std::ios::in | std::ios::out | std::ios::ate | std::ios::binary);
      }
    }

  private:

    // ATTRIBUTES

    //std::string const outputFile; // output file name
    std::string outputFile; // output file name
    std::ofstream outputStream; // output binary file stream
      // WARNING: Content of file is erased.

};

//wp: read and write plain text files //wp: note that the function is fully defined in a header and creates linking problems
//wp: inlining fixes this but if this were to be called N times, then basically the compiler would be copy pasting this code N times wherever is called, making the executable bigger and other possible problems
//wp: Another solution is to keep the header here and then write the readwrite.cpp file, but then would need to modify the make file...
inline std::vector<std::vector<double>> read_file(std::string name_file){ 
    
	using Sequence = std::vector<double>; //wp: alias of a type
	std::vector<Sequence> sequences;

	std::ifstream fin(name_file); //wp: file in

    {
        std::string first_line;
        std::getline(fin, first_line);
        std::istringstream iss(first_line); // used to separate each element in the line
        
        double element;
        while (iss >> element)
        {
            sequences.push_back(Sequence()); // add empty sequence
            sequences.back().push_back(element); // insert first element
        }
    }
    
    // First line and all sequences are now created.
    // Now we just loop for the rest of the way.
    bool end = false;
    while (!end)
    {
        for (size_t i = 0; i < sequences.size(); i++)
        {
            double element;
            if (fin >> element)
            {
                sequences[i].push_back(element);
            }
            else
            {
				//printf("\tAll data read in\n");
                // end of data. could do extra error checking after this to make sure the columns are all equal in size
                end = true;
                break;
            } 
        }
    } 
	return sequences;
} //wp: end of read_file function

#endif
