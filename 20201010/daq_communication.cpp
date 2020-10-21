#define _XOPEN_SOURCE_EXTENDED 1

#define MAX_LINE_LENGTH 1024
#define READING_BUF_SIZE 16384
#define MAX_PARAM_LENGTH 24
#define BUF_SIZE 40
#define RBCP_VER 0xFF
#define RBCP_ID 1
#define RBCP_PORT 0x1234
#define RBCP_CMD_WR 0x80
#define RBCP_CMD_RD 0xC0
#define TCP_PORT 0x18
#define UDP_BUF_SIZE 2048

#define RBCP_DISP_MODE_PROCESS 0
#define RBCP_DISP_MODE_WAVE 2
#define RBCP_DISP_MODE_STEADY 3

#define file_header "TRANSVERS"
#define data_header "wave"
#define file_footer "DATA processed with the 16-pu-monitor circuit"
#define data_footer "data"

#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <time.h>
#include <string>
#include <string.h>

using namespace std;

struct rbcp_header{
  unsigned char type;
  unsigned char command;
  unsigned char id;
  unsigned char length;
  unsigned int address;
};

int myGetArg(char* inBuf, int i, char* argBuf);
int myScanf(char* inBuf, char* argBuf1);
int rbcp_com(char* ipAddr, unsigned int port, struct rbcp_header* sendHeader, char* sendData, char* recvData, char dispMode);
int set_mode(char* ipAddr, unsigned int port, struct rbcp_header* sendHeader, char dispMode);
int DispatchCommand(char* pszVerb, char* ipAddr, unsigned int rbcpPort, struct rbcp_header* sndHeader, char dispMode);
int OnHelp();

int main() {
    cout << "╔═══════════════════════════════════════════════════════════╗" << endl;
    cout << "║    16Pick Up Beam Monitor DAQ connecting Program.         ║" << endl;
    cout << "║    This Program is built for collecting data for 16PU.    ║" << endl;
    cout << "║    If you have any connection error, please check out     ║" << endl;
    cout << "║    the ip address setting of script and firmware.         ║" << endl;
    cout << "║                                                           ║" << endl;
    cout << "║    Special Help : Toyama, Nakamura                        ║" << endl;
    cout << "║    Senior : Nakanishi, Uno, Tajima                        ║" << endl;
    cout << "║    Made by YOHAN LEE                                      ║" << endl;
    cout << "║                                                           ║" << endl;
    cout << "║                               Last update 2020/10/16      ║" << endl;
    cout << "╚═══════════════════════════════════════════════════════════╝" << endl;

    char* sitcpIpAddr = "10.72.108.42";
    unsigned int sitcpPort = 4660;
    struct rbcp_header sndHeader;
    char tempKeyBuf[MAX_LINE_LENGTH];
    char szVerb[MAX_PARAM_LENGTH];
    int rtnValue;

    sndHeader.type=RBCP_VER;
    sndHeader.id=0;
    
    while(1){
        printf("16PU-MONITOR@SiTCP-RBCP$ ");
        fgets(tempKeyBuf, MAX_LINE_LENGTH, stdin);
        if((rtnValue=myScanf(tempKeyBuf,szVerb))<0){
            printf("Erro myScanf(): %i\n",rtnValue);
            return -1;
        }
        sndHeader.id++;
        if(DispatchCommand(szVerb, sitcpIpAddr, sitcpPort, &sndHeader,1)<0) break;
    }
    return 0;
}

int myScanf(char* inBuf, char* argBuf1) {
    int i=0;

    argBuf1[0]='\0';

    if((i=myGetArg(inBuf, i, argBuf1))>0) {
        return i;
    } else {
        return i;
    }
}

int myGetArg(char* inBuf, int i, char* argBuf) {
    int j;

    while(i<MAX_LINE_LENGTH) {
        if(isblank(inBuf[i])) {
            i++;
        } else if(iscntrl(inBuf[i])) {
            return 0;
        } else {
            break;
        }
    }

    for(j=0;j<MAX_PARAM_LENGTH;j++) {
        if(i<MAX_LINE_LENGTH) {
            if(isblank(inBuf[i])) {
                argBuf[j]='\0';

            } else if(iscntrl(inBuf[i])) {
                argBuf[j]='\0';
                return 0;
            } else {
                argBuf[j]=inBuf[i];
                i++;
            }
        } else {
            return -1;
        }
    }

    return i;
}

int check_address(char *host_name) {
    string name = string(host_name);
    if (name.find("01") != std::string::npos) {
        return 13;
    } else if (name.find("02") != std::string::npos) {
        return 15;
    }
}

int rbcp_com(char* ipAddr, unsigned int port, struct rbcp_header* sendHeader, char* sendData, char* recvData, char dispMode){

    struct sockaddr_in sitcpAddr;
    int sock;

    struct timeval timeout;
    fd_set setSelect;
    
    int sndDataLen;
    int cmdPckLen;

    char sndBuf[1024];
    int address, i, j = 0;
    int rcvdBytes;
    char rcvdBuf[READING_BUF_SIZE];
    char host_name[64];
    int numReTrans =0;

    time_t rawtime;
    struct tm * timeinfo;
    char nFile [48];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    /* Create a Socket */
    if(dispMode==RBCP_DISP_MODE_STEADY) puts("\nRunning in steady mode...\n");

    sock = socket(AF_INET, SOCK_DGRAM, 0);
    address = gethostname(host_name, 64);

    cout << "You are connected in " << host_name << endl;

    sitcpAddr.sin_family      = AF_INET;
    sitcpAddr.sin_port        = htons(port);
    sitcpAddr.sin_addr.s_addr = inet_addr(ipAddr);

    sndDataLen = (int)sendHeader->length;

    /* Copy header data */
    memcpy(sndBuf,sendHeader, sizeof(struct rbcp_header));

    if(sendHeader->command==RBCP_CMD_WR){
        memcpy(sndBuf+sizeof(struct rbcp_header),sendData,sndDataLen);
        cmdPckLen=sndDataLen + sizeof(struct rbcp_header);
    }else{
        cmdPckLen=sizeof(struct rbcp_header);
    }

    /* send a packet*/

    sendto(sock, sndBuf, cmdPckLen, 0, (struct sockaddr *)&sitcpAddr, sizeof(sitcpAddr));
    if(dispMode==RBCP_DISP_MODE_STEADY) puts("The packet have been sent!\n");

    /* Receive packets*/
    
    if(dispMode==RBCP_DISP_MODE_STEADY) puts("\nWait to receive the ACK packet...");


    while(numReTrans<3){

    FD_ZERO(&setSelect);
    FD_SET(sock, &setSelect);

    timeout.tv_sec  = 1;
    timeout.tv_usec = 0;
 
    if(select(sock+1, &setSelect, NULL, NULL,&timeout)==0){
        /* time out */
        puts("\n***** Timeout ! *****");
        sendHeader->id++;
        memcpy(sndBuf,sendHeader, sizeof(struct rbcp_header));
        sendto(sock, sndBuf, cmdPckLen, 0, (struct sockaddr *)&sitcpAddr, sizeof(sitcpAddr));
        numReTrans++;
        FD_ZERO(&setSelect);
        FD_SET(sock, &setSelect);
    } else {
        /* receive packet */
        if(FD_ISSET(sock,&setSelect)){
            rcvdBytes=recvfrom(sock, rcvdBuf, READING_BUF_SIZE, 0, NULL, NULL);

            if(rcvdBytes<sizeof(struct rbcp_header)){
                puts("ERROR: ACK packet is too short");
                close(sock);
                return -1;
            }

            if((0x0f & rcvdBuf[1])!=0x8){
                puts("ERROR: Detected bus error");
                close(sock);
                return -1;
            }

            rcvdBuf[rcvdBytes]=0;
            int offset = 8;
            if((int)rcvdBuf[offset+0] == 0) puts("MODE : process");
            else if((int)rcvdBuf[offset+0] == 2) puts("MODE : wave");
            else if((int)rcvdBuf[offset+0] == 3) puts("MODE : steady");

            int data_length = 65536*((int)rcvdBuf[offset+33]) + 256*((int)rcvdBuf[offset+34]) + ((int)rcvdBuf[offset+35]);
            cout << "DATA # " << data_length << endl;

            if (rcvdBuf[offset+0] == 2 || rcvdBuf[offset+0] == 1) {
                if(rcvdBuf[offset+37] == 0) {cout << "sample memory is neither full nor empty." << endl;}
                else if(rcvdBuf[offset+37] == 1) {cout << "sample memory is full." << endl;}
                else if(rcvdBuf[offset+37] == 2) {cout << "sample memory is empty." << endl;}
            }

            if(RBCP_CMD_RD){
                memcpy(recvData,rcvdBuf+sizeof(struct rbcp_header),rcvdBytes-sizeof(struct rbcp_header));
            }

            if(dispMode==RBCP_DISP_MODE_STEADY){
                puts("\n***** A packet is received ! *****.");
                puts("Setting is successfully accepted");

            }else if(dispMode==RBCP_DISP_MODE_PROCESS){
                cout << "Start getting process data..." << endl;
                if (address == 13) strftime(nFile, 48, "process_%G_%m_%e_%H_%M_%S_address13.dat", timeinfo);
                else if (address == 15) strftime(nFile, 48, "process_%G_%m_%e_%H_%M_%S_address15.dat", timeinfo);
                while(1){

                }
                while(1){

                }
                cout << "Receiving data finished." << endl;
            }else if(dispMode==RBCP_DISP_MODE_WAVE){
                cout << "Start getting wave data..." << endl;
                if (address == 13) strftime(nFile, 20, "wave_%G_%m_%e_%H_%M_%S_address13.dat", timeinfo);
                else if (address == 15) strftime(nFile, 20, "wave_%G_%m_%e_%H_%M_%S_address15.dat", timeinfo);
                for (int ch = 0; ch < 15; ch++) {
                    while(1){

                    }
                    while(1){
                        
                    }
                }
                cout << "Receiving data finished." << endl;
            }
            numReTrans = 4;
            close(sock);
            return(rcvdBytes);
        }
    }
  }
  close(sock);

  return -3;
}

int set_mode(char* ipAddr, unsigned int port, struct rbcp_header* sendHeader, char dispMode) {
    struct sockaddr_in sitcpAddr;
    int sock;
    
    int sndDataLen;
    int cmdPckLen;

    char sndBuf[1024];
    int rcvdBytes;
    char rcvdBuf[2048];
    int sendData = 0xFF & dispMode;

    sock = socket(AF_INET, SOCK_DGRAM, 0);

    sitcpAddr.sin_family      = AF_INET;
    sitcpAddr.sin_port        = htons(port);
    sitcpAddr.sin_addr.s_addr = inet_addr(ipAddr);

    sndDataLen = (int)sendHeader->length;

    /* Copy header data */
    memcpy(sndBuf,sendHeader, sizeof(struct rbcp_header));

    if(sendHeader->command==RBCP_CMD_WR){
        memcpy(sndBuf+sizeof(struct rbcp_header), &sendData, sndDataLen);
        cmdPckLen=sndDataLen + sizeof(struct rbcp_header);
    }else{
        cmdPckLen=sizeof(struct rbcp_header);
    }

    /* send a packet*/
    sendto(sock, sndBuf, cmdPckLen, 0, (struct sockaddr *)&sitcpAddr, sizeof(sitcpAddr));
    rcvdBytes=recvfrom(sock, rcvdBuf, 2048, 0, NULL, NULL);

    return 0;
}

int DispatchCommand(char* pszVerb, char* ipAddr, unsigned int rbcpPort, struct rbcp_header* sndHeader, char dispMode){ 
    char recvData[UDP_BUF_SIZE];
    char *sendData;

    unsigned int tempInt;

    if(strcmp(pszVerb, "read") == 0){
        cout << "Select read out mode(0: process, 2: wave) : ";
        cin >> dispMode;
        sndHeader->command= RBCP_CMD_WR;
        sndHeader->length=1;
        sndHeader->address=htonl(0xFFFFFF10);
        return rbcp_com(ipAddr, rbcpPort, sndHeader, sendData, recvData, dispMode);
    }
    else if(strcmp(pszVerb, "trigger") == 0){
        sndHeader->command= RBCP_CMD_WR;
        sndHeader->length=1;
        sndHeader->address=htonl(37);
        sendData = "0x1";
        return rbcp_com(ipAddr, rbcpPort, sndHeader, sendData, recvData,3);
    }
    else if(strcmp(pszVerb, "reset") == 0){
        sndHeader->command= RBCP_CMD_WR;
        sndHeader->length=1;
        sndHeader->address=htonl(0xFFFFFF10);
        sendData = "0x0";
        return rbcp_com(ipAddr, rbcpPort, sndHeader, sendData, recvData,3);
    }
    else if(strcmp(pszVerb, "help") == 0){
        return OnHelp();
    }
    else if(strcmp(pszVerb, "quit") == 0){
        return -1;
    }
    puts("No such command!\n");
    return 0;
  
}

int OnHelp()
{
	puts("\nCommand list:");
	puts("   read\t\t:\tRead out wave/process binary data");
	puts("   trigger\t:\tSend test trigger singal");
	puts("   reset\t:\tReset SiTCP");
    puts("   help\t\t:\tShow available command");
	puts("   quit\t\t:\tquit from this program\n");

	return 0;
}