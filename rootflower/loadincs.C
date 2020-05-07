{
  TString wcsim_topdir   = gSystem->Getenv("WCSIMDIR");
  TString bonsai_topdir  = gSystem->Getenv("BONSAIDIR");
  TString flower_topdir = gSystem->Getenv("FLOWERDIR");

  TString mp = gROOT->GetMacroPath();
  TString ip;

  TString wcsim_inc   = wcsim_topdir   + "/include";
  TString bonsai_inc  = bonsai_topdir  + "/bonsai";
  TString flower_inc = flower_topdir;

  const char* p = wcsim_inc.Data();
  if (p) {
    mp += ":";
    mp += p;
    ip += " -I";
    ip += p;
  }
  else {
    cerr << "Could not find WCSim include path" << endl;
  }
  const char* b = bonsai_inc.Data();
  if (b) {
    mp += ":";
    mp += b;
    ip += " -I";
    ip += b;
  }
  else {
    cerr << "Could not find BONSAI include path" << endl;
  }
  const char* eb = flower_inc.Data();
  if (eb) {
    mp += ":";
    mp += eb;
    ip += " -I";
    ip += eb;
  }
  else {
    cerr << "Could not find FLOWER include path" << endl;
  }

  mp += ":/usr/local/include/";
  ip += "  -I/usr/local/include/";

  gROOT->SetMacroPath(mp.Data());
  gSystem->SetIncludePath(ip);

  // additions to .include must be done individually or CINT will
  // try to quote all the spaces as a single path

  TString dip = ".include ";
  dip += wcsim_inc.Data();
  gROOT->ProcessLine(dip.Data());

  dip = ".include ";
  dip += bonsai_inc.Data();
  gROOT->ProcessLine(dip.Data());

  dip = ".include ";
  dip += flower_inc.Data();
  gROOT->ProcessLine(dip.Data());

  dip = ".include /usr/local/include/";
  gROOT->ProcessLine(dip.Data());
}

