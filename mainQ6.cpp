//Counter for Date
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
const int SIZE = 300;
struct blocks {
    string start_code;
    string end_code;
    string name;
};
int main() {
    blocks BLOCK[SIZE];
    //Decimal numbers of start_code and end_code
    long start_cp[SIZE] = { 0 };
    long end_cp[SIZE] = { 0 };
    ifstream file("xxxxxx.txt", ifstream::in);
    if (!file) {
        cout << "Failed to open the file.\n";
        exit(1);
    }
    string line;
    int i = 0;
    while (getline(file, line)) {
        if (line[0] == '#' || line[0] == '\0' || line[0] == '\r') {
            continue;
        }
        istringstream record(line);
        vector<string> infos;
        string info1;
        string info2;
        string info3;
        while (getline(record, info1, '\r')) {
            istringstream record1(info1);
            while (getline(record1, info2, ';')) {
                infos.push_back(info2);
            }
            BLOCK[i].name = infos[1];
            BLOCK[i].name.erase(0, BLOCK[i].name.find_first_not_of(" "));
            BLOCK[i].name.erase(BLOCK[i].name.find_last_not_of(" ") + 1);
            istringstream record2(infos[0]);
            while (getline(record2, info3, '.')) {
                infos.push_back(info3);
            }
            BLOCK[i].start_code = infos[2].c_str();
            BLOCK[i].end_code = infos[4].c_str();
            start_cp[i] = strtol(BLOCK[i].start_code.c_str(), NULL, 16);
            end_cp[i] = strtol(BLOCK[i].end_code.c_str(), NULL, 16);
            i++;
        }
    }
    file.close();
	return 0;
}