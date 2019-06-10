/*SRIMscript.cpp
 *Takes a SRIM file and cleans it for use with anasen analysis. Mostly just makes everything have uniform units,
 *removes headers and unecessary info
 *
 *Maria A. -- a long time ago 
 *
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

int main(int argc, char **argv)
{
    string alouminum(argv[1]);//original filename
    
    vector<string> s2;//will store keV/MeV (unneccessary)
    vector<double> ionE,dE_dx_e, dE_dx_n;//will store values to save.
    string keV="keV";//this is just here to utilize in comparisons;
    
    ifstream in(alouminum.c_str());//open in file
    string input/*temporary*/,Straggling="Straggling"/*for comparison*/,commentline/*to save comment line*/;
    //This while loop reads the crap before the actual data.
    while (in>>input)//This saves in input automatically. if EOF will kill it (will never EOF)
    {
        if (input==Straggling)//Check if "Straggling is reached.
        {
            string junk;
            in>>input;//Second Straggling
            in>>commentline;//save this for later.
            getline(in, junk); //toss line of ---'s
            break;//kill loop since we can start reading data.
        }
        input.clear();
    }
    //Reads the rest until it reaches  a line that begins with '-'
    double ie;
    float dedxe;
    float dedxn;
    string units;
    string junk;
    while(in>>ie)
    {
        in>>units;
        if(units=="keV") ie = ie/1000.0;
        
        in>>dedxe;
        in>>dedxn;
        ionE.push_back(ie);
        dE_dx_e.push_back(dedxe);
        dE_dx_n.push_back(dedxn);
        
        in>>junk;
        in>>junk;
        in>>junk;
        in>>junk;
        in>>junk;
        in>>junk;
    }
     
    string outputfilename=alouminum.substr(0, alouminum.find("."));//output file name
    outputfilename += ".eloss";
    ofstream out(outputfilename.c_str());//open file.
    
    for(unsigned int i=0; i<ionE.size();i++)//save data.
    {
        out<<fixed<<showpoint;
        out<<setprecision(8);
        out<<setw(10)<<ionE[i]<<"\t"<<dE_dx_e[i]<<"\t"<<dE_dx_n[i]<<endl;
    }
    
	return 0;
}
