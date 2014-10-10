/* getmac.cpp
 *
 * Mattias Schevenels
 * May 2009
 */

// Windows Visual C++:  cl /EHsc getmac.cpp
// Windows GCC:         g++ -O3 -Wall -o getmac.exe getmac.cpp -liphlpapi
// Linux:               g++ -O3 -Wall -o getmac getmac.cpp

#ifdef _WIN32

#include <string>
#include <string.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <iphlpapi.h>
#include <stdio.h>
#include <stdlib.h>

#else

#include <string>
#include <string.h>
#include <stdio.h>
#include <sys/fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <net/if.h>

#endif

using namespace std;

void getmac(string& addr, const bool all)
{

#ifdef _WIN32 // WINDOWS
#ifdef __GNUC__ // WINDOWS GCC

  // declare and initialize variables
  DWORD dwSize = 0;
  DWORD dwRetVal = 0;
  DWORD maxSpeed = 0;
  MIB_IFTABLE *pIfTable;
  MIB_IFROW *pIfRow;
  char caddr[20];

  // initialize list of mac addresses
  addr = "";

  // allocate memory for our pointers
  pIfTable = (MIB_IFTABLE*)malloc(sizeof (MIB_IFTABLE));
  if (pIfTable==0) throw("Out of memory.");

  // make an initial call to GetIfTable to get the necessary size into dwSize
  dwSize = sizeof(MIB_IFTABLE);
  if (GetIfTable(pIfTable,&dwSize,false)==ERROR_INSUFFICIENT_BUFFER)
  {
    free(pIfTable);
    pIfTable = (MIB_IFTABLE *) malloc(dwSize);
    if (pIfTable==0) throw("Out of memory.");
  }

  // make a second call to GetIfTable to get the actual data we want.
  if ((dwRetVal = GetIfTable(pIfTable, &dwSize, FALSE)) == NO_ERROR)
  {
    for (int i=0; i<(int)pIfTable->dwNumEntries; i++)
    {
      pIfRow = (MIB_IFROW*)&pIfTable->table[i];
      if (pIfRow->dwPhysAddrLen>=6)
      {
        sprintf(caddr,"%2.2X%2.2X%2.2X%2.2X%2.2X%2.2X",(int)pIfRow->bPhysAddr[0],(int)pIfRow->bPhysAddr[1],(int)pIfRow->bPhysAddr[2],(int)pIfRow->bPhysAddr[3],(int)pIfRow->bPhysAddr[4],(int)pIfRow->bPhysAddr[5]);
        if (all)
        {
          if (addr.find(caddr)==string::npos)
          {
            if (addr.length()>0) addr.push_back(' ');
            addr.append(caddr);
          }
        }
        else
        {
          if ((pIfRow->dwType == IF_TYPE_ETHERNET_CSMACD) && (pIfRow->dwSpeed > maxSpeed))
          {
            addr=caddr;
            maxSpeed = pIfRow->dwSpeed;
          }
        }
      }
    }
  }
  if (pIfTable!=0)
  {
    free(pIfTable);
    pIfTable=0;
  }

#else // WINDOWS VISUAL C++
  #pragma comment(lib, "IPHLPAPI.lib")

  // Declare and initialize variables
  char caddr[20];
  PIP_ADAPTER_INFO pAdapterInfo;
  PIP_ADAPTER_INFO pAdapter = NULL;
  DWORD dwRetVal = 0;
  UINT i;

  // variables used to print DHCP time info
  struct tm newtime;
  char buffer[32];
  errno_t error;

  ULONG ulOutBufLen = sizeof (IP_ADAPTER_INFO);
  pAdapterInfo = (IP_ADAPTER_INFO *) malloc(sizeof (IP_ADAPTER_INFO));
  if (pAdapterInfo == NULL) throw("Unable to determine the computer's MAC address.");

  // Make an initial call to GetAdaptersInfo to get
  // the necessary size into the ulOutBufLen variable
  if (GetAdaptersInfo(pAdapterInfo, &ulOutBufLen) == ERROR_BUFFER_OVERFLOW)
  {
    free(pAdapterInfo);
    pAdapterInfo = (IP_ADAPTER_INFO *) malloc(ulOutBufLen);
    if (pAdapterInfo == NULL) throw("Unable to determine the computer's MAC address.");
  }

  if ((dwRetVal = GetAdaptersInfo(pAdapterInfo, &ulOutBufLen)) == NO_ERROR)
  {
    pAdapter = pAdapterInfo;
    while (pAdapter)
    {
      if (pAdapter->AddressLength>=6)
      {
        sprintf(caddr,"%2.2X%2.2X%2.2X%2.2X%2.2X%2.2X",(int)pAdapter->Address[0],(int)pAdapter->Address[1],(int)pAdapter->Address[2],(int)pAdapter->Address[3],(int)pAdapter->Address[4],(int)pAdapter->Address[5]);
        if (all)
        {
          if (addr.find(caddr)==string::npos)
          {
            if (addr.length()>0) addr.push_back(' ');
            addr.append(caddr);
          }
        }
        else
        {
          if (addr.empty() && (pAdapter->Type==MIB_IF_TYPE_ETHERNET))
          {
            addr=caddr;
          }
        }
      }
      pAdapter = pAdapter->Next;
    }
  }
  else throw("Unable to determine the computer's MAC address.");
  if (pAdapterInfo) free(pAdapterInfo);

#endif
#else // LINUX

  char devices[40960];
  int fid;
  int devlength;
  char* p1;
  char* p2;
  char* device;
  int sock;
  struct ifreq ifr;
  u_char uaddr[6];
  char caddr[20];

  // initialize list of mac addresses
  addr = "";

  // open socket
  sock = socket(AF_INET, SOCK_STREAM, 0);
  if (sock==-1) return;

  // load list of network devices
  if ((fid = open("/proc/net/dev",O_RDONLY))>=0)
  {
    devlength=read(fid, devices, sizeof(devices)-1);
    close(fid);
    devices[devlength]='\0';

    // loop over network devices
    p1=devices;
    while (p1!=0)
    {
      // extract line from /proc/net/dev
      p2=strchr(p1,'\n');
      if (p2)
      {
        *p2='\0';
        device=p1;
        p1=p2+1;
      }
      else
      {
        device=p1;
        p1=0;
      }

      // extract device name from line
      p2=strchr(device,':');
      if (p2)
      {
        *p2='\0';
        while ((strlen(device)>0) && (device[0]==' ')) device++;

        // get mac address
        strcpy(ifr.ifr_name,device);
        if ((ioctl(sock, SIOCGIFFLAGS, &ifr)==0) && (!(ifr.ifr_flags & IFF_LOOPBACK))
          && (ioctl(sock, SIOCGIFHWADDR, &ifr)==0) && (sizeof(ifr.ifr_hwaddr.sa_data)>=6))
        {
          memcpy(uaddr,ifr.ifr_hwaddr.sa_data,6);
          sprintf(caddr,"%2.2X%2.2X%2.2X%2.2X%2.2X%2.2X",uaddr[0],uaddr[1],uaddr[2],uaddr[3],uaddr[4],uaddr[5]);
          if (all)
          {
            if (addr.find(caddr)==string::npos)
            {
              if (addr.length()>0) addr.push_back(' ');
              addr.append(caddr);
            }
          }
          else
          {
            if (addr.empty() || (strcasecmp(device,"eth0")==0))
            {
              addr=caddr;
            }
          }
        }
      }
    }
  }

  // close socket
  close(sock);

#endif

  // check result
  if (addr.empty()) throw("Unable to determine the computer's MAC address. Please make sure that your Ethernet adapter is enabled.");
}



//int main()
//{
//  string addr;
//  getmac(addr,true);
//  printf("%s\n",addr.c_str());
//  return 0;
//}




