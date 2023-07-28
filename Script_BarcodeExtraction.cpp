//Oregon Health & Science University ("OHSU") owns the copyrights in this Software.  OHSU hereby licenses this Software under the Attribution-NonCommercial-NoDerivs 2.0 Generic license as provided here:    https://creativecommons.org/licenses/by-nc-nd/2.0/
//At line107, please specify the path for pool directory. 

#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <list>
#include <chrono>
#include <algorithm>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

using namespace std;

#if defined(_WIN32) || defined(_WIN64)

#include <windows.h>

char buffer[MAX_PATH];
ssize_t count = GetModuleFileName(NULL, buffer, MAX_PATH);
const string current_loc = string(buffer);
#elif defined(__APPLE__) || defined(__MACH__)
#include<mach-o/dyld.h>
char buffer[PATH_MAX];
uint32_t size = sizeof(buffer);
ssize_t count = _NSGetExecutablePath(buffer, &size);
const string current_loc = string( buffer);
#else
#include <unistd.h>
char buffer[PATH_MAX];
ssize_t count = readlink("/proc/self/exe", buffer, PATH_MAX);
const string current_loc = string( buffer);
#endif

vector <string> split(const string str, const string delim);

string revComp(string seq);

string replaceStrChar(string str, const string replace, char ch);

void print2dVec(const vector <vector<string>> vec);

int openDir(boost::filesystem::path &file_dir);

int openFilenames(const string &filenames_path, string (&filenames)[4]);

int
openSBCFile(const string (&filenames)[4], vector <vector<string>> &Sample_name_SBC_ID_SBC,
            vector <string> &Sample_name_SBC);

int openPrimersFile(const string (&filenames)[4], vector <string> &primers);

int openVBCFile(const string (&filenames)[4], vector <vector<string>> &VBC);

int openSeqData(list <string> &seqData, ifstream &inFile, const int &chunk);

void createPCRProducts(vector <vector<string>> &PCR_Products_BC, vector <vector<string>> &PCR_Products_BC_S,
                       vector <vector<string>> &PCR_Products_BC_AS,
                       const vector <vector<string>> &Sample_name_SBC_ID_SBC,
                       const vector <string> &primers, const vector <vector<string>> &VBC, const int type);

void
createSBCSearchTag(vector <vector<string>> &SBC_Search_Tag_BC, const vector <vector<string>> &Sample_name_SBC_ID_SBC,
                   const vector <string> &primers, const int type);

void sortSeq(map < string, vector < vector < string >> > &SBC_sorted_seq_BC,
const vector <vector<string>> &SBC_Search_Tag_BC,
const list <string> &seqData
);

void createVBCSearchTag(map < string, vector < vector < string >> > &VBC_Search_Tag_BC,
const vector <vector<string>> &Sample_name_SBC_ID_SBC,
const vector <vector<string>> &VBC,
const vector <string> &primers,
const int type
);

void initializeCount(map <string, map<string, int>> &count_BC1, map <string, map<string, int>> &count_BC2,
                     const vector <vector<string>> &Sample_name_SBC_ID_SBC, const vector <vector<string>> &VBC);

int findMatches(map <string, map<string, int>> &count_BC1, map <string, map<string, int>> &count_BC2,
                const vector <vector<string>> &Sample_name_SBC_ID_SBC,
                const map <string, vector<vector < string>>

> &SBC_sorted_seq_BC1,
const map <string, vector<vector < string>>> &SBC_sorted_seq_BC2,
const map <string, vector<vector < string>>> &VBC_Search_Tag_BC1,
const map <string, vector<vector < string>>> &VBC_Search_Tag_BC2,
const vector <vector<string>> &VBC
);

void outputData(const vector <string> &Sample_name_SBC, const vector <vector<string>> &VBC,
                const map <string, map<string, int>> &count_BC1,
                const map <string, map<string, int>> &count_BC2,
                const int &totalCounts, const int &totalReadCounts);

void outputPool(const vector <string> &Sample_name_SBC, const vector <vector<string>> &VBC,
                const map <string, map<string, int>> &count_BC1,
                const map <string, map<string, int>> &count_BC2,
                const boost::filesystem::path &pool_dir);


boost::filesystem::path file_dir("");
boost::filesystem::path pool_dir("");
const bool printReadFile = false;
const bool printValues = false;
const char PATH_SEPERATOR = boost::filesystem::path::preferred_separator;
int max_threads = boost::thread::hardware_concurrency();
boost::mutex mtx;

int main() {
    printf("Detected %d cores/threads...\n", max_threads);
    auto start = std::chrono::high_resolution_clock::now();
    string filenames[4];
    string filenames_path =
            current_loc.substr(0, current_loc.find_last_of(PATH_SEPERATOR) + 1) + "List_of_filenames.txt";
    if (openFilenames(filenames_path, filenames) == 0)return 0;
    file_dir.make_preferred();
    if (file_dir.string().back() != PATH_SEPERATOR) {
        file_dir += PATH_SEPERATOR;
    }
    printf("Attempting to use %d cores/threads...\n", max_threads);
    printf("Running in %s\n", file_dir.string().c_str());
    vector <vector<string>> Sample_name_SBC_ID_SBC;
    vector <string> Sample_name_SBC;
    if (openSBCFile(filenames, Sample_name_SBC_ID_SBC, Sample_name_SBC) == 0)return 0;

    vector <string> primers;
    if (openPrimersFile(filenames, primers) == 0) return 0;

    vector <vector<string>> VBC;
    if (openVBCFile(filenames, VBC) == 0) return 0;

    vector <vector<string>> PCR_Products_BC1, PCR_Products_BC1_S, PCR_Products_BC1_AS;
    createPCRProducts(PCR_Products_BC1, PCR_Products_BC1_S, PCR_Products_BC1_AS, Sample_name_SBC_ID_SBC, primers, VBC,
                      1);

    vector <vector<string>> PCR_Products_BC2, PCR_Products_BC2_S, PCR_Products_BC2_AS;
    createPCRProducts(PCR_Products_BC2, PCR_Products_BC2_S, PCR_Products_BC2_AS, Sample_name_SBC_ID_SBC, primers, VBC,
                      2);

    vector <vector<string>> SBC_Search_Tag_BC1;
    createSBCSearchTag(SBC_Search_Tag_BC1, Sample_name_SBC_ID_SBC, primers, 1);

    vector <vector<string>> SBC_Search_Tag_BC2;
    createSBCSearchTag(SBC_Search_Tag_BC2, Sample_name_SBC_ID_SBC, primers, 2);

    map < string, vector < vector < string >> > VBC_Search_Tag_BC1;
    createVBCSearchTag(VBC_Search_Tag_BC1, Sample_name_SBC_ID_SBC, VBC, primers, 1);

    map < string, vector < vector < string >> > VBC_Search_Tag_BC2;
    createVBCSearchTag(VBC_Search_Tag_BC2, Sample_name_SBC_ID_SBC, VBC, primers, 2);

    map <string, map<string, int>> count_BC1, count_BC2;
    initializeCount(count_BC1, count_BC2, Sample_name_SBC_ID_SBC, VBC);
    //if (openSeqData(filenames, seqData) == 0)return 0;

    printf("Opening data file... ");
    ifstream inFile;
    try {
        //inFile.open(file_dir.string() + filenames[3]);
        boost::filesystem::path seq_data_file(filenames[3]);
        seq_data_file.make_preferred();
        inFile.open(seq_data_file.string());
        if (!inFile.is_open()) {
            cout << "Error opening file!";
            return 0;
        }
    }
    catch (const ifstream::failure &e) {
        cout << "Error opening file!";
        return 0;
    }
    printf("Finished.\n");
    int chunk = 1;
    int totalCounts = 0;
    int totalReadCounts = 0;
    list <string> seqData;
    while (!inFile.eof()) {
        seqData.clear();
        openSeqData(seqData, inFile, chunk);
        ++chunk;
        totalReadCounts += seqData.size();

        map < string, vector < vector < string >> > SBC_sorted_seq_BC1;
        sortSeq(SBC_sorted_seq_BC1, SBC_Search_Tag_BC1, seqData);

        map < string, vector < vector < string >> > SBC_sorted_seq_BC2;
        sortSeq(SBC_sorted_seq_BC2, SBC_Search_Tag_BC2, seqData);

        int subTotalCounts = findMatches(count_BC1, count_BC2, Sample_name_SBC_ID_SBC, SBC_sorted_seq_BC1,
                                         SBC_sorted_seq_BC2,
                                         VBC_Search_Tag_BC1, VBC_Search_Tag_BC2, VBC);
        if (subTotalCounts < 0) return 0;
        totalCounts += subTotalCounts;
    }
    outputData(Sample_name_SBC, VBC, count_BC1, count_BC2, totalCounts, totalReadCounts);

    pool_dir.make_preferred();
    if (pool_dir.string().back() != PATH_SEPERATOR) {
        pool_dir += PATH_SEPERATOR;
    }
    if (boost::filesystem::exists(pool_dir)) {
        outputPool(Sample_name_SBC, VBC, count_BC1, count_BC2, pool_dir);
    } else {
        printf("Cannot find pool directory: %s\n", pool_dir.string().c_str());
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    printf("Process finished in %lld microseconds\n", duration.count());
}

int openFilenames(const string &filenames_path, string (&filenames)[4]) {
    printf("Reading filenames... ");
    ifstream inFile;
    try {
        inFile.open(filenames_path);
        if (!inFile.is_open()) {
            cout << "Error opening file!";
            return 0;
        }
        int filenumber = -2;
        for (string filename; getline(inFile, filename);) {
            if (filename.substr(0, 1) == "#") {
                continue;
            }
            if (filenumber < -1) {
                filename = boost::regex_replace(filename, boost::regex("^ +| +$|\t|\n|\r"), "");
                max_threads = stoi(filename);
            }
            else if (filenumber < 0) {
                filename = boost::regex_replace(filename, boost::regex("^ +| +$|\t|\n|\r"), "");
                file_dir = filename;
            } else {
                filename = boost::regex_replace(filename, boost::regex("^ +| +$|\t|\n|\r"), "");
                filenames[filenumber] = filename;
            }
            ++filenumber;
        }
        inFile.close();
    }
    catch (const ifstream::failure &e) {
        cout << "Error opening file!";
        return 0;
    }
    printf("Finished.\n");
    return 1;
}

int
openSBCFile(const string (&filenames)[4], vector <vector<string>> &Sample_name_SBC_ID_SBC,
            vector <string> &Sample_name_SBC) {
    printf("Reading SBC file... ");
    ifstream inFile;
    try {
        inFile.open(file_dir.string() + filenames[0]);
        if (!inFile.is_open()) {
            cout << "Error opening file!";
            return 0;
        }
        for (string line; getline(inFile, line);) {
            if (line.substr(0, 1) == "#") {
                continue;
            }
            line = boost::regex_replace(line, boost::regex("^ +| +$|\n|\r"), "");
            vector <string> line1 = split(line, "\t");
            string SBC_rev_comp = revComp(line1[2]);
            vector <string> line2(line1.begin(), line1.end());
            line2.push_back(SBC_rev_comp);
            string line3 = line1[0] + "_" + line1[1];
            Sample_name_SBC_ID_SBC.push_back(line2);
            Sample_name_SBC.push_back(line3);
        }
        inFile.close();
    }
    catch (const ifstream::failure &e) {
        cout << "Error opening file!";
        return 0;
    }
    printf("Finished.\n");
    if (printReadFile)print2dVec(Sample_name_SBC_ID_SBC);
    return 1;
}

int openPrimersFile(const string (&filenames)[4], vector <string> &primers) {
    printf("Reading primers file... ");
    ifstream inFile;
    try {
        inFile.open(file_dir.string() + filenames[1]);
        if (!inFile.is_open()) {
            cout << "Error opening file!";
            return 0;
        }
        for (string line; getline(inFile, line);) {
            if (line.substr(0, 1) == "#") {
                continue;
            }
            line = boost::regex_replace(line, boost::regex("^ +| +$|\t|\n|\r"), "");
            primers.push_back(line);
        }
        inFile.close();
    }
    catch (const ifstream::failure &e) {
        cout << "Error opening file!";
        return 0;
    }
    primers[1] = revComp(primers[1]);
    primers[3] = revComp(primers[3]);
    printf("Finished.\n");
    if (printReadFile) {
        for (string primer: primers)printf("%s\n", primer.c_str());
    }
    return 1;
}

int openVBCFile(const string (&filenames)[4], vector <vector<string>> &VBC) {
    printf("Reading VBC file... ");
    ifstream inFile;
    try {
        inFile.open(file_dir.string() + filenames[2]);
        if (!inFile.is_open()) {
            cout << "Error opening file!";
            return 0;
        }
        for (string line; getline(inFile, line);) {
            if (line.substr(0, 1) == "#") {
                continue;
            }
            line = boost::regex_replace(line, boost::regex("^ +| +$|\n|\r"), "");
            vector <string> line1 = split(line, "\t");
            VBC.push_back(line1);
        }
        inFile.close();
    }
    catch (const ifstream::failure &e) {
        cout << "Error opening file!";
        return 0;
    }
    printf("Finished.\n");
    if (printReadFile)print2dVec(VBC);
    return 1;
}

int openSeqData(list <string> &seqData, ifstream &inFile, const int &chunk) {
    printf("Reading sequence data chunk: %d... ", chunk);
    vector <string> rawData;
    rawData.reserve(20000000);
    int lineNum = 0;
    for (string line; getline(inFile, line);) {
        rawData.push_back(line);
        ++lineNum;
        if (lineNum == 20000000)break;
    }
    printf("Finished\n");
    boost::regex expSpec("[#@\\+]*");
    boost::regex expNuc("[ACGTN]{55,}");
    int threadNum = 0;
    int chunkSize = ceil((float) rawData.size() / max_threads);
    auto parseData = [&]() {
        mtx.lock();
        int num = threadNum++;
        mtx.unlock();
        list <string> parsedThread;
        for (int i = num * (chunkSize); i < min((num + 1) * (chunkSize), (int) rawData.size()); ++i) {
            string line = rawData[i];
            boost::smatch matchNuc, matchSpec;
            if (boost::regex_match(line, matchSpec, expSpec) || !boost::regex_match(line, matchNuc, expNuc)) continue;
            line = boost::regex_replace(line, boost::regex("^ +| +$|\n|\r"), "");
            parsedThread.push_back(line);
        }
        mtx.lock();
        seqData.splice(seqData.end(), parsedThread);
        mtx.unlock();
    };
    printf("Parsing sequence data chunk:%d... ", chunk);
    boost::thread threads[max_threads];
    for (int j = 0; j < max_threads; j++) {
        threads[j] = boost::thread(parseData);
    }
    for (int j = 0; j < max_threads; j++) {
        threads[j].join();
    }
    printf("Finished.\n");
    return 1;
}

void createPCRProducts(vector <vector<string>> &PCR_Products_BC, vector <vector<string>> &PCR_Products_BC_S,
                       vector <vector<string>> &PCR_Products_BC_AS,
                       const vector <vector<string>> &Sample_name_SBC_ID_SBC,
                       const vector <string> &primers, const vector <vector<string>> &VBC, const int type) {
    printf("Creating PCR products... ");
    for (vector <string> sampleBarcode: Sample_name_SBC_ID_SBC) {
        string sampleName = sampleBarcode[0] + "_" + sampleBarcode[1];
        for (vector <string> virusBarcode:VBC) {
            string PCR_S;
            if (type == 1)
                PCR_S = sampleBarcode[2] + primers[0] + "GCTAGC" + virusBarcode[1] + primers[1] + sampleBarcode[3];
            else if (type == 2)
                PCR_S = sampleBarcode[2] + primers[2] + "GCTAGC" + virusBarcode[2] + primers[3] + sampleBarcode[3];
            string PCR_AS = revComp(PCR_S);
            vector <string> VBCID_PCR_S;
            VBCID_PCR_S.push_back(virusBarcode[0]);
            VBCID_PCR_S.push_back(PCR_S);
            vector <string> VBCID_PCR_AS;
            VBCID_PCR_AS.push_back(virusBarcode[0]);
            VBCID_PCR_AS.push_back(PCR_AS);
            vector <string> line_S, line_AS;
            line_S.push_back(sampleName);
            line_S.insert(line_S.end(), VBCID_PCR_S.begin(), VBCID_PCR_S.end());
            line_AS.push_back(sampleName);
            line_AS.insert(line_AS.end(), VBCID_PCR_AS.begin(), VBCID_PCR_AS.end());
            PCR_Products_BC_S.push_back(line_S);
            PCR_Products_BC_AS.push_back(line_AS);
        }
    }
    for (int i = 0; i < PCR_Products_BC_S.size(); ++i) {
        vector <string> line;
        line.push_back(PCR_Products_BC_S[i][0]);
        line.push_back(PCR_Products_BC_S[i][1]);
        line.push_back(PCR_Products_BC_S[i][2]);
        line.push_back(PCR_Products_BC_AS[i][2]);
        PCR_Products_BC.push_back(line);
    }
    if (printValues)print2dVec(PCR_Products_BC);
    printf("Finished.\n");
}

void
createSBCSearchTag(vector <vector<string>> &SBC_Search_Tag_BC, const vector <vector<string>> &Sample_name_SBC_ID_SBC,
                   const vector <string> &primers, const int type) {
    printf("Creating SBC BC%d search tag... ", type);
    for (vector <string> sampleBarcode: Sample_name_SBC_ID_SBC) {
        string sampleName = sampleBarcode[0] + "_" + sampleBarcode[1];
        string PCR_S;
        if (type == 1)PCR_S = sampleBarcode[2] + primers[0] + primers[1] + sampleBarcode[3];
        else if (type == 2) PCR_S = sampleBarcode[2] + primers[2] + primers[3] + sampleBarcode[3];
        string PCR_AS = revComp(PCR_S);
        string searchTagS = PCR_S.substr(1, 12);
        string searchTagAS = PCR_AS.substr(1, 12);
        vector <string> line;
        line.push_back(sampleName);
        line.push_back(searchTagS);
        line.push_back(searchTagAS);
        SBC_Search_Tag_BC.push_back(line);
    }
    printf("Finished.\n");
    if (printValues)print2dVec(SBC_Search_Tag_BC);
}

void sortSeq(map < string, vector < vector < string >> > &SBC_sorted_seq_BC,
const vector <vector<string>> &SBC_Search_Tag_BC,
const list <string> &seqData
) {
printf("Sorting sequence data... ");
int threadNum = 0;
auto thread_search = [&]() {
    mtx.lock();
    int num = threadNum++;
    mtx.unlock();
    int chunkSize = ceil((float) SBC_Search_Tag_BC.size() / max_threads);
    for (int i = num * (chunkSize);
         i < min((num + 1) * (chunkSize), (int) SBC_Search_Tag_BC.size()); ++i) {
        vector <string> SBC_sorted_seq_BC_S, SBC_sorted_seq_BC_AS;
        string search_tag_S = SBC_Search_Tag_BC[i][1];
        string search_tag_AS = SBC_Search_Tag_BC[i][2];
        for (string sequence:seqData) {
            if (boost::algorithm::contains(sequence, search_tag_S))SBC_sorted_seq_BC_S.push_back(sequence);
            else if (boost::algorithm::contains(sequence, search_tag_AS)) SBC_sorted_seq_BC_AS.push_back(sequence);
        }
        mtx.lock();
        SBC_sorted_seq_BC[SBC_Search_Tag_BC[i][0]].push_back(SBC_sorted_seq_BC_S);
        SBC_sorted_seq_BC[SBC_Search_Tag_BC[i][0]].push_back(SBC_sorted_seq_BC_AS);
        mtx.unlock();
    }
};
vector <boost::thread> threads;
for (
int j = 0;
j<max_threads;
j++) {
threads.
push_back(boost::thread(thread_search)
);
}
for (
int j = 0;
j<max_threads;
j++) {
threads[j].

join();

}
printf("Finished.\n");
}

void createVBCSearchTag(map < string, vector < vector < string >> > &VBC_Search_Tag_BC,
const vector <vector<string>> &Sample_name_SBC_ID_SBC,
const vector <vector<string>> &VBC,
const vector <string> &primers,
const int type
) {
printf("Creating SBC BC%d search tag... ", type);
for (
vector <string> sampleBarcode
: Sample_name_SBC_ID_SBC) {
string sampleName = sampleBarcode[0] + "_" + sampleBarcode[1];
for (
vector <string> virusBarcode
: VBC) {
vector <string> temp;
string PCR_S;
if (type == 1)
PCR_S = sampleBarcode[2] + primers[0] + "GCTAGC" + virusBarcode[1] + primers[1] + sampleBarcode[3];
else if (type == 2)
PCR_S = sampleBarcode[2] + primers[2] + virusBarcode[2] + "TGTACA" + primers[3] + sampleBarcode[3];
string PCR_AS = revComp(PCR_S);
PCR_S = PCR_S.substr(1, 50);
PCR_AS = PCR_AS.substr(1, 50);
temp.
push_back(virusBarcode[0]);
temp.
push_back(PCR_S);
temp.
push_back(PCR_AS);
VBC_Search_Tag_BC[sampleName].
push_back(temp);
}
}
printf("Finished.\n");
if (printValues) {
for (
auto const &x
: VBC_Search_Tag_BC)
print2dVec(x
.second);
}
}

void initializeCount(map <string, map<string, int>> &count_BC1, map <string, map<string, int>> &count_BC2,
                     const vector <vector<string>> &Sample_name_SBC_ID_SBC, const vector <vector<string>> &VBC) {
    for (int i = 0; i < Sample_name_SBC_ID_SBC.size(); ++i) {
        string sampleName = Sample_name_SBC_ID_SBC[i][0] + "_" + Sample_name_SBC_ID_SBC[i][1];
        for (int j = 0; j < VBC.size(); ++j) {
            string VBC_ID = VBC[j][0];
            count_BC1[sampleName][VBC_ID] = 0;
            count_BC2[sampleName][VBC_ID] = 0;
        }
    }
}

int findMatches(map <string, map<string, int>> &count_BC1, map <string, map<string, int>> &count_BC2,
                const vector <vector<string>> &Sample_name_SBC_ID_SBC,
                const map <string, vector<vector < string>>

> &SBC_sorted_seq_BC1,
const map <string, vector<vector < string>>> &SBC_sorted_seq_BC2,
const map <string, vector<vector < string>>> &VBC_Search_Tag_BC1,
const map <string, vector<vector < string>>> &VBC_Search_Tag_BC2,
const vector <vector<string>> &VBC
) {
int totalCounts = 0;
for (
int i = 0;
i<Sample_name_SBC_ID_SBC.

size();

++i) {
string sampleName = Sample_name_SBC_ID_SBC[i][0] + "_" + Sample_name_SBC_ID_SBC[i][1];
printf("Currently counting VBCs (BC1) of: %s... ", sampleName.

c_str()

);
vector <string> seq_data_S = SBC_sorted_seq_BC1.at(sampleName)[0];
vector <string> seq_data_AS = SBC_sorted_seq_BC1.at(sampleName)[1];
map<string, int> VBC_check;
for(
int j = 0;
j<VBC.

size();

j++){
VBC_check[VBC[j][0]] = 0;
}
int threadNum = 0;
int chunkSize = ceil((float) VBC.size() / max_threads);
auto thread_search_BC1 = [&]() {
    mtx.lock();
    int num = threadNum++;
    mtx.unlock();
    for (int j = num * chunkSize; j < min((int) VBC.size(), (num + 1) * chunkSize); ++j) {
        string VBC_ID = VBC_Search_Tag_BC1.at(sampleName)[j][0];
        string search_tag_S = VBC_Search_Tag_BC1.at(sampleName)[j][1];
        string search_tag_AS = VBC_Search_Tag_BC1.at(sampleName)[j][2];
        int count = 0;
        for (string sequence:seq_data_S) {
            if (sequence.find(search_tag_S) != string::npos) {
                ++count;
            }
        }
        for (string sequence:seq_data_AS) {
            if (sequence.find(search_tag_AS) != string::npos) {
                ++count;
            }
        }
        mtx.lock();
        VBC_check[VBC_ID] += 1;
        count_BC1[sampleName][VBC_ID] += count;
        totalCounts += count;
        mtx.unlock();
    }
};
boost::thread threads[max_threads];
for (
int j = 0;
j<max_threads;
j++) {
threads[j] =
boost::thread(thread_search_BC1);
}
for (
int j = 0;
j<max_threads;
j++) {
threads[j].

join();

}
for(
int j = 0;
j<VBC.

size();

j++){
if(VBC_check[VBC[j][0]]!=1){
printf("ABORT. Duplicate/missing data detected. Please consult following DEBUG:\n");
for(
int k = 0;
k<VBC.

size();

++k){
string VBC_ID = VBC[k][0];
printf("%s %d\n", VBC_ID.

c_str(), VBC_check[VBC_ID]

);
}
return -1;
}
}
printf("Finished.\n");
printf("Currently counting VBCs (BC2) of: %s... ", sampleName.

c_str()

);
seq_data_S = SBC_sorted_seq_BC2.at(sampleName)[0];
seq_data_AS = SBC_sorted_seq_BC2.at(sampleName)[1];
for(
int j = 0;
j<VBC.

size();

j++){
VBC_check[VBC[j][0]] = 0;
}
threadNum = 0;
auto thread_search_BC2 = [&]() {
    mtx.lock();
    int num = threadNum++;
    mtx.unlock();
    for (int j = num * (chunkSize); j < min((int) VBC.size(), (num + 1) * (chunkSize)); ++j) {
        string VBC_ID = VBC_Search_Tag_BC2.at(sampleName)[j][0];
        string search_tag_S = VBC_Search_Tag_BC2.at(sampleName)[j][1];
        string search_tag_AS = VBC_Search_Tag_BC2.at(sampleName)[j][2];
        int count = 0;
        for (string sequence:seq_data_S) {
            if (sequence.find(search_tag_S) != string::npos) {
                ++count;
            }
        }
        for (string sequence:seq_data_AS) {
            if (sequence.find(search_tag_AS) != string::npos) {
                ++count;
            }
        }
        mtx.lock();
        VBC_check[VBC_ID] += 1;
        count_BC2[sampleName][VBC_ID] += count;
        totalCounts += count;
        mtx.unlock();
    }
};
for (
int j = 0;
j<max_threads;
j++) {
threads[j] =
boost::thread(thread_search_BC2);
}
for (
int j = 0;
j<max_threads;
j++) {
threads[j].

join();

}
for(
int j = 0;
j<VBC.

size();

j++){
if(VBC_check[VBC[j][0]]!=1){
printf("ABORT. Duplicate/missing data detected. Please consult following DEBUG:\n");
for(
int k = 0;
k<VBC.

size();

++k){
string VBC_ID = VBC[k][0];
printf("%s %d\n", VBC_ID.

c_str(), VBC_check[VBC_ID]

);
}
return -1;
}
}
printf("Finished.\n");
}
return
totalCounts;
}

void outputData(const vector <string> &Sample_name_SBC, const vector <vector<string>> &VBC,
                const map <string, map<string, int>> &count_BC1,
                const map <string, map<string, int>> &count_BC2,
                const int &totalCounts, const int &totalReadCounts) {
    printf("Creating output files... ");
    ofstream outputFile;
    boost::filesystem::path dir(file_dir.string());
    if (!boost::filesystem::exists(file_dir.string() + "OUTPUT")) {
        boost::filesystem::create_directory(file_dir.string() + "OUTPUT");
    }
    if (!boost::filesystem::exists(file_dir.string() + "OUTPUT" + PATH_SEPERATOR + "Pool")) {
        boost::filesystem::create_directory(file_dir.string() + "OUTPUT" + PATH_SEPERATOR + "Pool");
    }
    for (string sampleName: Sample_name_SBC) {
        string filename_BC1 = sampleName + "_BC1.txt";
        string filename_BC2 = sampleName + "_BC2.txt";
        outputFile.open(file_dir.string() + "OUTPUT" + PATH_SEPERATOR + "Pool" + PATH_SEPERATOR + filename_BC1);
        outputFile << "\t" << sampleName + "_BC1" << endl;
        for (vector <string> x: VBC) {
            outputFile << x[0] << "\t" << count_BC1.at(sampleName).at(x[0]) << "\n";
        }
        outputFile.close();
        outputFile.open(file_dir.string() + "OUTPUT" + PATH_SEPERATOR + "Pool" + PATH_SEPERATOR + filename_BC2);
        outputFile << "\t" << sampleName + "_BC2" << endl;
        for (vector <string> x: VBC) {
            outputFile << x[0] << "\t" << count_BC2.at(sampleName).at(x[0]) << "\n";
        }
        outputFile.close();
    }
    outputFile.open(file_dir.string() + "OUTPUT" + PATH_SEPERATOR + "Summary.txt");
    outputFile << "Total read counts: " << totalReadCounts << endl;
    outputFile << "Total matched counts: " << totalCounts << endl;
    outputFile.close();
    printf("Finished.\n");
}

void outputPool(const vector <string> &Sample_name_SBC, const vector <vector<string>> &VBC,
                const map <string, map<string, int>> &count_BC1,
                const map <string, map<string, int>> &count_BC2,
                const boost::filesystem::path &pool_dir) {
    printf("Writing pool data... ");
    ofstream outputFile;
    for (string sampleName: Sample_name_SBC) {
        string filename_BC1 = sampleName + "_BC1.txt";
        string filename_BC2 = sampleName + "_BC2.txt";
        outputFile.open(pool_dir.string() + filename_BC1);
        outputFile << "\t" << sampleName + "_BC1" << endl;
        for (vector <string> x: VBC) {
            outputFile << x[0] << "\t" << count_BC1.at(sampleName).at(x[0]) << "\n";
        }
        outputFile.close();
        outputFile.open(pool_dir.string() + filename_BC2);
        outputFile << "\t" << sampleName + "_BC2" << endl;
        for (vector <string> x: VBC) {
            outputFile << x[0] << "\t" << count_BC2.at(sampleName).at(x[0]) << "\n";
        }
        outputFile.close();
    }
    printf("Finished.\n");
}

vector <string> split(const string str, const string delim) {
    vector <string> parts;
    size_t start, end = 0;
    while (end < str.size()) {
        start = end;
        while (start < str.size() && (delim.find(str[start]) != string::npos)) {
            ++start;  // skip initial whitespace
        }
        end = start;
        while (end < str.size() && (delim.find(str[end]) == string::npos)) {
            ++end; // skip to end of word
        }
        if (end - start != 0) {  // just ignore zero-length strings.
            parts.push_back(string(str, start, end - start));
        }
    }
    return parts;
}

char complement(char n) {
    switch (n) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
    }
    return ' ';
}

string revComp(string seq) {
    string seq_AS(seq);
    reverse(seq_AS.begin(), seq_AS.end());
    transform(seq_AS.begin(), seq_AS.end(), seq_AS.begin(), complement);
    return seq_AS;
}

string replaceStrChar(string str, const string replace, char ch) {
    size_t found = str.find_first_of(replace);
    while (found != string::npos) {
        str[found] = ch;
        found = str.find_first_of(replace, found + 1);
    }
    return str;
}

void print2dVec(const vector <vector<string>> vec) {
    for (vector <string> x: vec) {
        for (string str: x) {
            printf("%s ", str.c_str());
        }
        printf("\n");
    }
}
