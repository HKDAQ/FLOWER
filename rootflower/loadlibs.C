{
  TString libs0  = gSystem->GetDynamicPath();
  TString libswc = gSystem->Getenv("WCSIM_BUILD_DIR");
  libswc += "/lib";
  TString libsb  = gSystem->Getenv("BONSAIDIR");
  TString libseb = gSystem->Getenv("FLOWERDIR");
  TString libs   = libseb + ":" + libsb + ":" + libswc + ":" + libs0 + ":/usr/lib:/usr/local/lib:/opt/lib:/opt/local/lib";
  gSystem->SetDynamicPath(libs.Data());

  gSystem->Load("libGpad");
  gSystem->Load("libPhysics");
  gSystem->Load("libMatrix");
  gSystem->Load("libHist");
  gSystem->Load("libGraf");
  gSystem->Load("libTree");
  gSystem->Load("libRIO");
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");
  //gSystem->Load("libMinuit2");
  gSystem->Load("libMathMore"); 

  gSystem->Load("libWCSimRoot.so");
  
  gSystem->Load("libWCSimBonsai.so");

  gSystem->Load("libWCSimFLOWER.so");
}

