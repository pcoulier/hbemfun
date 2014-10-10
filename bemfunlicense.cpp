/* bemfunlicense.cpp
 *
 * Mattias Schevenels
 * July 2009
 */

#include <ctime>
#include "rsa.h"
#include "ripemd128.h"
#include "getmac.h"
#include "checklicense.h"
#include "isinf.h"
#include "mex.h"
#include "string.h"

#ifdef _WIN32
#include <windows.h>
#endif

#ifndef __GNUC__
#define strcasecmp _strcmpi
#endif

using namespace std;

/******************************************************************************/

struct LICENSE
{
  int status;
  time_t lastverified;
  string text;
  string number;
  string bemfunversion;
  int bemfunbuild;
  int expdate;
  bool trialmode;
  bool renewable;
  string message;
  string macaddress;
  string signature;
};

static LICENSE license = {-1,0,"","","",0,0,false,false,"","",""};

/******************************************************************************/

void strtrim(string& str)
{
  size_t startpos = str.find_first_not_of(" \t\f\v\n\r");
  size_t endpos = str.find_last_not_of(" \t\f\v\n\r");
  if((string::npos==startpos) || (string::npos==endpos))
    {
      str = "";
    }
    else
    {
      str = str.substr(startpos,endpos-startpos+1);
    }
}

/******************************************************************************/

void trimleadingzeros(string& str)
{
  size_t startpos = str.find_first_not_of("0 \t\f\v\n\r");
  size_t endpos = str.find_last_not_of(" \t\f\v\n\r");
  if((string::npos==startpos) || (string::npos==endpos))
    {
      str = "";
    }
    else
    {
      str = str.substr(startpos,endpos-startpos+1);
    }
}

/******************************************************************************/

void txtread(const string& file, string& str)
{
  str.clear();
  FILE* fid = fopen(file.c_str(),"r");
  if (fid!=NULL)
  {
    char c;
    while (true)
    {
      c = fgetc(fid);
      if (c==EOF) break;
      str.push_back(c);
    }
    fclose (fid);
  }
}

/******************************************************************************/

void strtodate(string str, int& date)
{
  string day;
  string month;
  string year;
  struct tm* timeinfo;
  time_t rawtime;
  unsigned int i;
  double gdate;
  double diff;

  i=str.find('/');
  day=str.substr(0,i);

  str.erase(0,i+1);
  i=str.find('/');
  month=str.substr(0,i);
  str.erase(0,i+1);
  year=str;

  timeinfo = gmtime(&rawtime);
  timeinfo->tm_year = atoi(year.c_str())-1900;
  timeinfo->tm_mon = atoi(month.c_str())-1;
  timeinfo->tm_mday = atoi(day.c_str());
  gdate=(double)mktime(timeinfo);

  time_t now=time(0);
  tm* t;
  t = localtime(&now);
  time_t ltime = mktime(t);
  t = gmtime(&now);
  time_t gtime = mktime(t);
  diff=difftime(gtime,ltime);

  date = (int)((gdate-diff+1)/3600/24+719529);
}


/******************************************************************************/

void datetostr(int date, string& str)
{
  time_t t;
  t=(date-719529)*3600*24;
  char c[40];
  strftime(c,40,"%d %B %Y",localtime(&t));
  str=c;
  if (str[0]=='0') str.erase(0,1);
}

/******************************************************************************/

void readlictxt(string& lictxt) // Find and read license file
{
  string licfile;
  unsigned int k;

  mxArray* lhs[1];
  mxArray* rhs[1];
  rhs[0] = mxCreateString(mexFunctionName());
  mexCallMATLAB(1,lhs,1,rhs,"which");
  char* c = new(nothrow) char[mxGetNumberOfElements(lhs[0])+1];
  if (c==0) throw("Out of memory.");
  mxGetString(lhs[0],c,mxGetNumberOfElements(lhs[0])+1);
  string licdir(c);
  k = licdir.find_last_of("/\\");
  mxDestroyArray(lhs[0]);
  mxDestroyArray(rhs[0]);
  delete [] c;
  lictxt.clear();

  while (lictxt.empty() && (k!=(unsigned int)string::npos))
  {
    licdir = licdir.substr(0,k);
    #ifdef _WIN32
    licfile = licdir + "\\license.dat";
    #else
    licfile = licdir + "/license.dat";
    #endif

    //mexPrintf("%s\n",licfile.c_str());
    txtread(licfile,lictxt);
    k = licdir.find_last_of("/\\");
    if (lictxt.find("BEMFUN LICENSE")==string::npos) lictxt.clear();
  }
}

/******************************************************************************/

void trimlictxt(string& lictxt) // Extract the license part involving BEMFUN
{
  unsigned int i;
  unsigned int j;
  unsigned int k;
  string line;
  string newlic;
  bool copy;

  // Extract relevant part from license
  newlic.clear();
  i = 0;
  j = lictxt.find('\n');
  copy = false;
  while (j!=(unsigned int)string::npos)
  {
    line=lictxt.substr(i,j-i);
    strtrim(line);
    if (!line.empty())
    {
      k = line.find(':');
      if (k==(unsigned int)string::npos) copy = (strcasecmp(line.c_str(),"BEMFUN LICENSE")==0);
      if (copy)
      {
        newlic.append(line);
        newlic.push_back('\n');
      }
    }
    i=j;
    j=lictxt.find('\n',j+1);
  }
  lictxt = newlic;
}

/******************************************************************************/

void striplictxt(string &lictxt) // Strip the signature from the license
{
  string newlic = "";
  unsigned int i;
  unsigned int j;
  unsigned int k;
  string line;
  string keyname;
  string keyvalue;

  i = 0;
  j = lictxt.find('\n');
  while (j!=(unsigned int)string::npos)
  {
    line=lictxt.substr(i,j-i);
    strtrim(line);
    k = line.find(':');
    if ((k!=(unsigned int)string::npos) && (k>0))
    {
      keyname = line.substr(0,k);
      strtrim(keyname);
      if (strcasecmp(keyname.c_str(),"Signature")!=0)
      {
        keyvalue = line.substr(k+1);
        strtrim(keyvalue);
        newlic.append(keyname);
        newlic.append(": ");
        newlic.append(keyvalue);
        newlic.push_back('\n');
      }
    }
    i=j;
    j=lictxt.find('\n',j+1);
  }
  lictxt = newlic;
}

/******************************************************************************/

void getkey(const string& keys, const string& keyname, string& keyvalue)
{
  unsigned int i;
  unsigned int j;
  unsigned int k;
  string line;
  string linename;

  keyvalue.clear();
  i = 0;
  j = keys.find('\n');
  while (j!=(unsigned int)string::npos)
  {
    line=keys.substr(i,j-i);
    strtrim(line);
    k = line.find(':');
    if ((k!=(unsigned int)string::npos) && (k>0))
    {
      linename = line.substr(0,k);
      strtrim(linename);
      if (strcasecmp(linename.c_str(),keyname.c_str())==0)
      {
        keyvalue = line.substr(k+1);
        strtrim(keyvalue);
        return;
      }
    }
    i=j;
    j=keys.find('\n',j+1);
  }
}

/******************************************************************************/

void readlicense(LICENSE& license)
{
  string lictxt;
  string bemfunbuild;
  string expdate;
  string renewable;
  string trialmode;

  license.status = -1;
  license.lastverified = 0;

  readlictxt(lictxt);
  trimlictxt(lictxt);
  getkey(lictxt,"License number",license.number);
  getkey(lictxt,"BEMFUN version",license.bemfunversion);
  getkey(lictxt,"BEMFUN build",bemfunbuild);
  license.bemfunbuild = atoi(bemfunbuild.c_str());
  getkey(lictxt,"Expiration date",expdate);
  strtodate(expdate,license.expdate);
  getkey(lictxt,"Trial mode",trialmode);
  license.trialmode = (strcasecmp(trialmode.c_str(),"yes") == 0);
  getkey(lictxt,"Renewable",renewable);
  license.renewable = (strcasecmp(renewable.c_str(),"yes") == 0);
  getkey(lictxt,"Message",license.message);
  if (!license.message.empty())
  {
    unsigned int i;
    i = license.message.find("\\n"); while (i != (unsigned int)string::npos) { license.message.replace(i,2,"\n"); i = license.message.find("\\n"); }
    i = license.message.find("\\r"); while (i != (unsigned int)string::npos) { license.message.replace(i,2,"\r"); i = license.message.find("\\r"); }
    i = license.message.find("\\t"); while (i != (unsigned int)string::npos) { license.message.replace(i,2,"\t"); i = license.message.find("\\t"); }
    i = license.message.find("\\b"); while (i != (unsigned int)string::npos) { license.message.replace(i,2,"\b"); i = license.message.find("\\b"); }
    i = license.message.find("\\f"); while (i != (unsigned int)string::npos) { license.message.replace(i,2,"\f"); i = license.message.find("\\f"); }
    i = license.message.find("\\\\");while (i != (unsigned int)string::npos) { license.message.replace(i,2,"\\"); i = license.message.find("\\\\");}
    i = license.message.find("%%");  while (i != (unsigned int)string::npos) { license.message.replace(i,2,"%");  i = license.message.find("%%");  }
  }
  getkey(lictxt,"MAC address",license.macaddress);
  getkey(lictxt,"Signature",license.signature);
  striplictxt(lictxt);
  license.text = lictxt;
}

/******************************************************************************/

int today()
{
  time_t now=time(0);
  tm* t;
  t = localtime(&now);
  time_t ltime = mktime(t);
  t = gmtime(&now);
  time_t gtime = mktime(t);

  int date=(int)((now-difftime(gtime,ltime))/3600/24+719529);

  #ifdef _WIN32

  HKEY hKey = 0;
  char subkey[27];
  subkey[0] = 'S';
  subkey[1] = 'o';
  subkey[2] = 'f';
  subkey[3] = 't';
  subkey[4] = 'w';
  subkey[5] = 'a';
  subkey[6] = 'r';
  subkey[7] = 'e';
  subkey[8] = '\\';
  subkey[9] = 'B';
  subkey[10] = 'W';
  subkey[11] = 'M';
  subkey[12] = '\\';
  subkey[13] = 'm';
  subkey[14] = 'e';
  subkey[15] = 'e';
  subkey[16] = 'u';
  subkey[17] = 'w';
  subkey[18] = 'e';
  subkey[19] = 'n';
  subkey[20] = 'v';
  subkey[21] = 'r';
  subkey[22] = 'i';
  subkey[23] = 'e';
  subkey[24] = 'n';
  subkey[25] = 'd';
  subkey[26] = 0;
  int regdate = 0;
  const int refdate = 733547;
  DWORD regdatesize = sizeof(regdate);

  RegCreateKeyEx(HKEY_CURRENT_USER,subkey,0,0,0,KEY_ALL_ACCESS,0,&hKey,0);
  if (RegQueryValueEx(hKey,0,0,0,(LPBYTE)&regdate,&regdatesize)==ERROR_SUCCESS)
  {
    if (regdate < 300) throw("BEMFUN license error.");
    regdate += refdate;
    if (regdate > date) date = regdate;
  }
  regdate = date - refdate;
  RegSetValueEx(hKey,0,0,REG_DWORD,(LPBYTE)&regdate,sizeof(regdate));
  RegCloseKey(hKey);

  #endif

  return date;
}

/******************************************************************************/

string url(const string& link, const string& disp)
{
  string browser;
  string result;

  #ifdef _WIN32
  result="<a href=\"matlab:bemfunDUMMY=system('cmd.exe /c rundll32 url.dll,FileProtocolHandler "+link+"');clear('bemfunDUMMY');\">"+disp+"</a>";
  #else
  if (system(NULL))
  {
    if (system("which xdg-open")==0) browser = "xdg-open";
    else if (system("which gnome-open")==0) browser = "gnome-open";
    else if (system("which exo-open")==0) browser = "exo-open";
    else if (system("which kfmclient")==0) browser = "kfmclient exec";
    else if (system("which firefox")==0) browser = "firefox";
    else if (system("which mozilla")==0) browser = "mozilla";
    else if (system("which netscape")==0) browser = "netscape";
    else if (system("which konqueror")==0) browser = "konqueror";
  }

  if (browser.empty()) browser = "none";
  if (strcasecmp(browser.c_str(),"none")==0)
  {
    result=disp;
  }
  else
  {
    result="<a href=\"matlab:bemfunLD_LIBRARY_PATH=getenv('LD_LIBRARY_PATH');setenv('LD_LIBRARY_PATH','');[bemfunDUMMY1,bemfunDUMMY2]=system('"+browser+" ''"+link+"''');setenv('LD_LIBRARY_PATH',bemfunLD_LIBRARY_PATH);clear('bemfunDUMMY1','bemfunDUMMY2','bemfunLD_LIBRARY_PATH');\">"+disp+"</a>";
  }
  #endif
  return result;
}

/******************************************************************************/

int bemfunlicense(const char* const action)
{
  // If the license status is known and the action is VerifyOnce, suppress output
  const bool quiet = ((license.status!=-1) && (strcasecmp(action,"VerifyOnce")==0));

  // Reset license status to unknown every 2 hours
  if (time(0)-license.lastverified>=7200) license.status = -1;

  // If requested and allowed, return immediately.
  if ((license.status==0) && (strcasecmp(action,"VerifyOnce")==0)) return license.status;

  // If requested, reset license status and return.
  if (strcasecmp(action,"Reset")==0)
  {
    license.status = -1; // Unknown
    mexUnlock();
    return license.status;
  }

  // Lock mex file to prevent static variable from being deleted
  if (!mexIsLocked()) mexLock();

  // If the license status is unknown, read the license file.
  if (license.status==-1) readlicense(license);

  // If the license status is unknown, determine it
  if (license.status==-1)
  {
    license.lastverified=time(0); // Verifying now
    if (license.text.empty()) license.status = 1; // License file not found
  }

  // If the license file is found, check license signature
  if (license.status==-1)
  {
    string rsaexp = "10001";
    string rsamod = "A7E460AFADB544E156180D853ECA0940EF015C17F968B4069C9638C4625708B9";
    string hash1;
    string hash2;

    rsadec(license.signature,rsaexp,rsamod,hash1);
    rmdstring(license.text,hash2);
    trimleadingzeros(hash1);
    trimleadingzeros(hash2);
    if (strcasecmp(hash1.c_str(),hash2.c_str())!=0) license.status = 4; // Invalid license file
  }

  // If the license file is valid, check the license term
  if (license.status==-1)
  {
    if (license.expdate<today()) license.status = 2; // Expired
  }

  // If the license term is still running, check the mac address
  if (license.status==-1)
  {
    if (!license.macaddress.empty())
    {
      license.status = 3; // Assume incorrect mac address
      unsigned int i;
      string macaddress;
      getmac(macaddress,true);
      while (!macaddress.empty())
      {
        i=macaddress.find(" ");
        if (i==(unsigned int)string::npos) i = macaddress.length();
        if (license.macaddress.find(macaddress.substr(0,i))!=string::npos)
        {
          license.status = -1 ; // Undo incorrect mac address assumption
          break;
        }
        macaddress.erase(0,i+1);
      }
    }
  }

  // If the mac address is OK, the license status is OK
  if (license.status==-1) license.status = 0;

  // Show message/error depending on license status
  if (license.status==0) // License OK, show message from license file
  {
    if (!quiet)
    {
      if (license.trialmode) // trial mode
      {
        string macaddress;
        getmac(macaddress,false);
        string expdate;
        datetostr(license.expdate,expdate);

        string msg = license.message;
        if (!msg.empty()) msg = msg+"------------------------------------------------------------------------\n";
        msg="------------------------------------------------------------------------\n"
            "BEMFUN: MATLAB toolbox for Boundary elements in Elastodynamics\n"
            "Version 2.1 Build 16 (July 2009)\n" // BEMFUNVERSION BEMFUNBUILD
            "------------------------------------------------------------------------\n"
            + msg +
            "You are using BEMFUN in trial mode.\n"
            "Your trial license expires on "+expdate+".\n"
            "If you want to continue using BEMFUN after this date, you have to obtain a\n"
            "permanent license. You can order a permanent license on the following\n"
            "web page: "+url("http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html?macaddress="+macaddress,"http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html")+".\n"
            "On this web page, you will be asked for a MAC address to identify the\n"
            "computer you want a license for. This computer's MAC address is the\n"
            "following: "+macaddress+".\n\n";

        mexPrintf(msg.c_str());


      }
      else // full mode
      {
        string msg = license.message;
        int daysLeft = license.expdate-today();
        if (daysLeft<=14)
        {
          char c[10];
          sprintf(c,"%i",daysLeft);
          if (!msg.empty()) msg=msg+"------------------------------------------------------------------------\n";
          if (daysLeft==0) msg=msg+"Your license file expires today.\n";
          else if (daysLeft==1) msg=msg+"Your license file will expire tomorrow.\n";
          else msg=msg+"Your license file will expire in "+c+" days.\n";

          if (license.renewable) msg=msg+"You can obtain a new license file, free of charge, from the following\nweb page: ";
          else msg=msg+"You can order a new license file on the following web page:\n";
          msg=msg+url("http://bwk.kuleuven.be/apps/bwm/bemfun/renew.html?licensenumber="+license.number,"http://bwk.kuleuven.be/apps/bwm/bemfun/renew.html")+".\n"
                  "On this web page, you will be asked for your license number.\n"
                  "Your license number is "+license.number+".\n";
        }
        if (!msg.empty())
        {
          msg="------------------------------------------------------------------------\n"
              "BEMFUN: MATLAB toolbox for Boundary elements in Elastodynamics\n"
              "Version 2.1 Build 16 (July 2009)\n" // BEMFUNVERSION BEMFUNBUILD
              "------------------------------------------------------------------------\n" + msg + "\n";
          mexPrintf(msg.c_str());
        }
      }
    }
  }
  else if (license.status==1) // License file not found
  {
    string macaddress;
    string msg;
    getmac(macaddress,false);
    msg="A license file is required to use BEMFUN.\n"
        "You can either buy a permanent license, or request a free trial license for one month.\n\n"
        "A permanent license can be ordered here: "+url("http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html?macaddress="+macaddress,"http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html")+".\n\n"
        "A trial license can be requested here: "+url("http://bwk.kuleuven.be/apps/bwm/bemfun/try.html?macaddress="+macaddress,"http://bwk.kuleuven.be/apps/bwm/bemfun/try.html")+".\n\n"
        "On both web pages, you will be asked for a MAC address to identify the computer you want a license for.\n"
        "This computer's MAC address is the following: "+macaddress+".";
    throw(msg.c_str());
  }
  else if (license.status==2) // License expired
  {
    if (license.trialmode)
    {
      string macaddress;
      string msg;
      getmac(macaddress,false);
      msg = "Your trial license for BEMFUN has expired.\n"
            "If you want to continue using BEMFUN, you have to obtain a permanent license.\n\n"
            "A permanent license can be ordered here: "+url("http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html?macaddress="+macaddress,"http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html")+".\n\n"
            "On this web page, you will be asked for a MAC address to identify the computer you want a license for.\n"
            "This computer's MAC address is the following: "+macaddress+".";
      throw(msg.c_str());
    }
    else // full mode
    {
      string msg;
      msg="Your BEMFUN license file has expired. ";
      if (license.renewable) msg=msg+"You can obtain a new license file, free of charge, from the following web page: ";
      else msg=msg+"You can order a new license file on the following web page: ";
      msg=msg+url("http://bwk.kuleuven.be/apps/bwm/bemfun/renew.html?licensenumber="+license.number,"http://bwk.kuleuven.be/apps/bwm/bemfun/renew.html")+". "
              "On this web page, you will be asked for your license number. "
              "Your license number is "+license.number+".";
      throw(msg.c_str());
    }
  }
  else if (license.status==3) // Incorrect mac address
  {
    string msg;
    string macaddress;
    getmac(macaddress,false);
    msg="You do not have a license to use BEMFUN on this computer.\n"
        "You can either buy a permanent license, or request a free trial license for one month.\n\n"
        "A permanent license can be ordered here: "+url("http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html?macaddress="+macaddress,"http://bwk.kuleuven.be/apps/bwm/bemfun/buy.html")+".\n\n"
        "A trial license can be requested here: "+url("http://bwk.kuleuven.be/apps/bwm/bemfun/try.html?macaddress="+macaddress,"http://bwk.kuleuven.be/apps/bwm/bemfun/try.html")+".\n\n"
        "On both web pages, you will be asked for a MAC address to identify the computer you want a license for. "
        "This computer's MAC address is the following: "+macaddress+".";
    throw(msg.c_str());
  }
  else if (license.status==4) // Invalid license file
  {
    throw("You do not have a valid license to use BEMFUN.");
  }
  else
  {
    throw("Unable to determine license status.");
  }

  return license.status;
}
