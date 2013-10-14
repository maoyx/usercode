#include <fstream>
#include <sstream>
#include <iostream>
#include "TGraphErrors.h"
#include "TCanvas.h"

using namespace std;

bool isData(string line);

int checkModes(string line, bool& doFocus, bool& doSys,bool& symStat,bool& symSys);

void checkTable(string filename = "figure1.txt", bool logCanvas = 0){
  bool rePrint = 0;
  ifstream in( filename.c_str() );
  string line;

  int iPlot = -1;
  int nPlot = -1;
  int mode = 0;

  TGraphErrors* gSys[20];
  TGraphErrors* gStat[20];

  int ipoint = 0;

  bool doSys = 0;
  bool symSys = 1;
  bool symStat = 1;
  bool doFocus = 0;

  bool start = 0;

  while ( getline( in, line ) ) {
    if ( !line.size() || line[0]=='#' ) continue;

    if ( line[0]=='x' ){

      mode = checkModes(line,doFocus,doSys,symStat,symSys);

      start = 1;
      iPlot++;
      nPlot++;
      gSys[iPlot] = new TGraphErrors();
      gStat[iPlot] = new TGraphErrors();

      ipoint = 0;
      continue;
    }

    if ( !isData(line) ){
      start = 0;
    }
   
    if(!start) continue;

    if(rePrint){ 
      cout<<"Reading point : "<<endl;
      cout<<line.data()<<endl;
    }
    istringstream ss(line);

    double xmin;
    double xmax;
    double y;
    double stat;
    double sys;

    double statup;
    double sysup;
    double statdown;
    double sysdown;

    double xfocus;

    if(mode == 0){
      ss >> xmin >> xmax >> y >> stat >> sys;
      xfocus = (xmax+xmin)/2.;
    }

    if(mode == 1){
      ss >> xmin >> xmax >> xfocus >> y >> stat >> sys;
    }

    gSys[iPlot]->SetPoint(ipoint,xfocus,y);
    gStat[iPlot]->SetPoint(ipoint,xfocus,y);

    gSys[iPlot]->SetPointError(ipoint,0,stat);
    gStat[iPlot]->SetPointError(ipoint,0,sys);


    ipoint++;

  }


  for(int i = 0; i < nPlot; ++i){

    TCanvas* c1 = new TCanvas(Form("c%d",i),Form("c%d",i),400,400);
  c1->SetLogy(logCanvas);

  gStat[i]->Draw("ap");
  if(doSys){
    gSys[i]->Draw("2");
  }
  
  }

}


bool isData(string line){

  bool data = 0;
  if(line.size() == 0) return 0;

  if(line[0] >= '0' && line[0] <= '9') data = 1;  
  if(line[0] == 'x') data = 1;

  return data;
}


int checkModes(string line, bool& doFocus, bool& doSys,bool& symStat,bool& symSys){
  
  int mode = 0;
  doFocus = line.find("xfocus", 0) != string::npos;
  doSys = line.find("sys", 0) != string::npos;
  symSys = !(line.find("+sys", 0) != string::npos);
  symStat = !(line.find("+stat", 0) != string::npos);

  if(!doFocus && !doSys && symStat && symSys) mode = 0;
  if( doFocus && !doSys && symStat && symSys) mode = 1;
  if( !doFocus && doSys && symStat && symSys) mode = 2;
  if( doFocus && doSys && symStat && symSys) mode = 3;
  if( !doFocus && doSys && symStat && !symSys) mode = 4;
  if( !doFocus && doSys && !symStat && !symSys) mode = 5;
  if( doFocus && doSys && symStat && !symSys) mode = 6;
  if( doFocus && doSys && !symStat && !symSys) mode = 7;

  return mode;
}

