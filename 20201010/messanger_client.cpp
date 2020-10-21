#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <string.h>
using namespace std;

#pragma comment(lib, "WSock32.lib")

int main()
{
    cout << "********** start chat program **********" << endl;
    cout << "If you want to quit, type 'end'" << endl;

    //IPv4のソケットアドレス情報を設定する
    struct sockaddr_in sockadd;
    memset(&sockadd, 0, sizeof(sockadd));
    sockadd.sin_family = AF_INET;
    sockadd.sin_port = htons(4660);
    sockadd.sin_addr.s_addr = inet_addr("192.168.11.8");

    const auto sock = socket(PF_INET, SOCK_STREAM, 0);

    //指定したソケットへ接続する
    connect(sock, (struct sockaddr *)&sockadd, sizeof(sockadd));

    //メッセージを送受信
    string buffer;
    while (true)
    {
        cout << "You\t> ";
        cin >> buffer;
        send(sock, buffer.c_str(), buffer.size(), 0);
        if (buffer == "end")
        {
            break;
        }
        cout << "Waiting for reply...\n";
        char buf[1024];
        memset(buf, '\0', sizeof(buf));
        const auto rec = recv(sock, buf, sizeof(buf)-1, 0);
        if (rec == -1)
        {
            cout << "error\n";
            break;
        }
        buf[rec] = '\0';
        if (strcmp(buf, "end") == 0)
        {
            cout << "server teminated.\n";
            break;
        }

        cout << "Someone\t> " << buf << endl;
    }
    close(sock);
    cout << "End\n";

    return 0;
}