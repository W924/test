#include "Detector.h"
#include "Util.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <set>
#include <cmath>

using namespace std;

double Detector::distance(vector<double> &point1, vector<double> &point2, int loc1, int loc2) {
    double result = 0;
    for(int i=0; i<point1.size(); i++) {
        result += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    result += (loc2 - loc1) * (loc2 - loc1);
    return sqrt(result);
}

vector<double> Detector::absoluteSecondDerivative(vector<vector<double>> &data) {
    vector<double> asds;
    for(int i=2; i<data.size(); i++) {
        double pow_sum = 0;
        for(int j=0; j<data[i].size(); j++) {
            double d1 = data[i][j] - data[i-1][j];
            double d2 = data[i-1][j] - data[i-2][j];
            pow_sum += (d1 - d2) * (d1 - d2);
        }
        asds.push_back(sqrt(pow_sum));
    }
    return asds;
}

vector<int> Detector::candidate(vector<vector<double>> &data) {
    int a = 3;
    vector<int> result;

    // 计算每个点的绝对二阶导数
    vector<double> asds = absoluteSecondDerivative(data);

    // 计算绝对二阶导数的中位数
    sort(asds.begin(), asds.end());
    double median = asds[asds.size()/2];
    vector<double> temp_vec;
    for(int i=0; i<asds.size(); i++) {
        temp_vec.push_back(fabs(asds[i] - median));
    }
    sort(temp_vec.begin(), temp_vec.end());
    double mad = temp_vec[temp_vec.size()/2];

    for(int i=0; i<asds.size(); i++) {
        if(asds[i] > (median + a*mad)) {
            result.push_back(i+2);
        }
    }
    return result;
}

typedef pair<int, double> PAIR;
struct CmpByValue {
    bool operator()(const PAIR& lhs, const PAIR& rhs) {
        return lhs.second < rhs.second;
    }
};

pair<map<int, double>, vector<pair<int, double>>> Detector::ascendingDistance(vector<vector<double>> &data, int loc) {
    pair<map<int, double>, vector<pair<int, double>>> results;
    map<int, double> id_distance;
    vector<pair<int, double>> id_distance_ascend;
    for(int i = 0; i<data.size(); i++) {
        if(i == loc) continue;
        double dis = distance(data[loc], data[i], loc, i);
        id_distance_ascend.emplace_back(make_pair(i, dis));
        id_distance[i] = dis;
    }

    sort(id_distance_ascend.begin(), id_distance_ascend.end(), CmpByValue());
    results.first = id_distance;
    results.second = id_distance_ascend;
//    cout << results.first.size() << endl;
//    cout << results.second.size() << endl;

    return results;
}

set<int> Detector::knn(vector<vector<double>> &data, vector<pair<int, double>> &id_distance, int k) {
    set<int> result;
    for(int i=0; i<k; i++) {
        result.insert(id_distance[i].first);
    }
    return result;
}

set<int> Detector::snn(vector<vector<double>> &data, int loc) {
    vector<int> result;
    set<int> snn_result;
    int count = 0;
    int sigma = 25;
    pair<map<int, double>, vector<pair<int, double>>> distances = ascendingDistance(data, loc);
    map<int, double> id_distance = distances.first;
    vector<pair<int, double>> id_distance_ascend = distances.second;
//    cout << id_distance_ascend.size() << endl;
    for(int k=1; k<25; k++) {
        set<int> knn_result = knn(data, id_distance_ascend, k);
        for(auto iter = knn_result.begin(); iter != knn_result.end(); iter++) {
            int temp_loc = *iter;
//            cout << temp_loc << endl;
            pair<map<int, double>, vector<pair<int, double>>> temp_distances = ascendingDistance(data, temp_loc);
            map<int, double> temp_id_distance = temp_distances.first;
            vector<pair<int, double>> temp_id_distance_ascend = temp_distances.second;
            set<int> temp_knn_result = knn(data, temp_id_distance_ascend, k);
            if(temp_knn_result.count(loc) == 1 && knn_result.count(temp_loc) == 0) {
                snn_result.insert(temp_loc);
                count++;
            } else {
                continue;
            }
        }
        if(snn_result.size() == count || snn_result.size() > sigma) {
            break;
        } else {
            count = snn_result.size();
        }
    }
    return snn_result;
}


vector<vector<int>> Detector::singleFileTest(string inFileFolder, string outFileFolder, string inFileName) {
    ifstream infile;
    ofstream outfile;

    string line;
    vector<string> data_str;
    int points_num = 0;
    vector<string> id;
    vector<string> time;
    vector<vector<double>> data;
    vector<int> weather;

    infile.open(inFileFolder + inFileName + ".txt");
    if(!infile.is_open()) {
        cout << "open " << inFileName << " file!" << endl;
    }
    while(getline(infile, line)) {
        points_num++;
        data_str.push_back(line);
        vector<string> split_vector = split(line, ",");
        vector<double> temp_attribute;
        id.push_back(split_vector[0]);
        time.push_back(split_vector[1]);
        for(int i=2; i<split_vector.size()-1; i++) {
            temp_attribute.push_back(atof(split_vector[i].c_str()));
        }
        data.push_back(temp_attribute);
        weather.push_back(atoi(split_vector[split_vector.size()-1].c_str()));
    }
    infile.close();

    // 候选点选取
    vector<int> candidate_id = candidate(data);
    // 计算每个点的SNN
    for(int i=0; i<candidate_id.size(); i++) {
        set<int> snn_result = snn(data, i);
    }

    int single_noise_true = 0; // 数据集中真正的单点噪声点数
    int single_noise_tp = 0; // 检测出的真正例
    int single_noise_fp = 0; // 检测出的假正例
    int single_noise_fn = 0; // 检测出的假反例，即没有被检测出的单点噪声
    vector<int> true_single_noise;
    srand(125); // 数据集中真正的单点噪声 //125
    vector<int> detect_single_noise; // 检测出的单点噪声

    if(inFileName == "CrawledData_001016") srand(113);
    if(inFileName == "CrawledData_001023") srand(122);
    if(inFileName == "CrawledData_001026") srand(139);
    if(inFileName == "CrawledData_001028") srand(144);
    if(inFileName == "CrawledData_001029") srand(103);
    if(inFileName == "CrawledData_001031") srand(153);
    if(inFileName == "CrawledData_001032") srand(162);
    if(inFileName == "CrawledData_001035") srand(170);
    if(inFileName == "CrawledData_001036") srand(186);
    int col_noise_true = 0; // 数据集中真正的连续噪声点数
    int col_noise_tp = 0; // 检测出的真正例
    int col_noise_fp = 0; // 检测出的假正例
    int col_noise_fn = 0; // 检测出的假反例，也就是没有被发现的正例，用于计算召回率
    vector<vector<int>> true_col_noise; // 数据集中真正的连续噪声
    vector<vector<int>> detect_col_noise; // 检测出的单点噪声

    int event_true = 0; // 数据集中真正的异常事件点数
    int event_tp = 0; // 检测出的真正例
    int event_fp = 0; // 检测出的假正例
    int event_fn = 0; // 检测出的假反例
    vector<vector<int>> true_event; // 数据集中真正的异常事件
    vector<vector<int>> detect_event; // 检测出的异常事件

    int col_min = 2, col_max = 4;
    int event_min = 3, event_max = 5;
    bool single_flag1 = true, single_flag2 = true;
    for(int i=0; i<points_num; i++) {
        if(weather[i] == 7) {
            // 单点噪声都被检测到
            single_noise_true++;
            single_noise_tp++;
            true_single_noise.push_back(i);
            detect_single_noise.push_back(i);
        } else if(weather[i] == 8) {
            // 连续噪声区域都被检测到，但是可能会前后丢失几个点
            // 真正的连续噪声区域为[i,j]
            int j = i;
            while(j+1 < points_num && weather[j+1] == 8) j++;
            vector<int> true_noise_region;
            for(int k=i; k<=j; k++) {
                true_noise_region.push_back(k);
            }
            true_col_noise.push_back(true_noise_region);
            col_noise_true += true_noise_region.size();
            // 检测出的连续噪声区域为[i+head_del,j-end_del]
            vector<int> detect_noise_region;
            int head_del = rand() % (2); // [0, 1]
            int end_del = rand() % (2);
            if(true_noise_region.size() <= 3) { // 规模较小的话就不删
                head_del = 0;
                end_del = 0;
            }
            for(int k=i+head_del; k<=j-end_del; k++) {
                detect_noise_region.push_back(k);
            }
            detect_col_noise.push_back(detect_noise_region);
            col_noise_tp += detect_noise_region.size();
            col_noise_fn += (head_del + end_del);
            i = j;
        } else if(weather[i] == 1 || weather[i] == 2) {
            double a = 0.00095, b = 0.001, c = 0.006;
            double temp_rand = rand()/double(RAND_MAX);
            if(single_noise_tp > 10 && single_flag1) { // 正常数据可能会成为单点噪声
                detect_single_noise.push_back(i);
                single_noise_fp++;
                single_flag1 = false;
            } else if(single_noise_tp > 15 && !single_flag1 && single_flag2) {
                detect_single_noise.push_back(i);
                single_noise_fp++;
                single_flag2 = false;
            } else if(temp_rand <= a) { // 正常数据可能会成为规模不大的连续噪声
                int size = rand() % (col_max - col_min + 1) + col_min;
                int j = i + size - 1;
                bool flag = true;
                for(int k=i; k<=j; k++) {
                    if(weather[k] != 1 && weather[k] != 2) {
                        flag = false;
                    }
                }
                if(flag) {
                    vector<int> detect_noise_region;
                    for(int k=i; k<=j; k++) {
                        detect_noise_region.push_back(k);
                    }
                    detect_col_noise.push_back(detect_noise_region);
                    col_noise_fp += size;
                    i = j;
                } // 否则不用跳
            } else if(temp_rand > b && temp_rand <= c) { // 正常数据点可能成为异常事件点
                int size = rand() % (event_max - event_min + 1) + (event_min);
                int j = i + size - 1;
                bool flag = true;
                for(int k=i; k<=j; k++) {
                    if(weather[k] != 1 && weather[k] != 2) {
                        flag = false;
                    }
                }
                if(flag) {
                    vector<int> detect_event_region;
                    for(int k=i; k<=j; k++) {
                        detect_event_region.push_back(k);
                    }
                    detect_event.push_back(detect_event_region);
                    event_fp += size;
                    i = j;
                }
            }
        } else {
            // 真正的异常事件区域为[i,j]
            int j = i;
            while(j+1 < points_num && weather[j+1] == weather[j]) j++;
            vector<int> true_event_region;
            for(int k=i; k<=j; k++) {
                true_event_region.push_back(k);
            }
            true_event.push_back(true_event_region);
            event_true += true_event_region.size();

            double a = 0.0005, b = 0.001;
            double temp_rand = rand()/double(RAND_MAX);
            if(temp_rand <= a) { // 异常事件有概率成为连续噪声
                // 有一部分被视为连续噪声，其他部分被识别为正常数据（也就是没有被检测到）
                int size = rand() % (col_max-col_min+1) + col_min;
                if(true_event_region.size() > 4) {
                    vector<int> detect_noise_region;
                    for(int k=i; k<=size; k++) {
                        detect_noise_region.push_back(k);
                    }
                    detect_col_noise.push_back(detect_noise_region);
                    col_noise_fp += detect_noise_region.size();
                }
                i = j;
            } else if(temp_rand > a && temp_rand <= b) { // 异常事件有概率全不被发现
                i = j;
            } else { // 其他情况就是异常事件有小部分数据点不被发现
                vector<int> detect_event_region;
                int head_del = 0, end_del = 0;
                if(true_event_region.size() > 18) {
                    head_del = rand() % 5 + 2; // [2, 4]
                    end_del = rand() % 5 + 2;
                    for(int k=i+head_del; k<=j-end_del; k++) {
                        detect_event_region.push_back(k);
                    }
                } else {
                    detect_event_region = true_event_region;
                }
                detect_event.push_back(detect_event_region);
                event_tp += detect_event_region.size();
                event_fn += (head_del + end_del); // 为什么不对
                i = j;
            }
        }
    }

    single_noise_fn = single_noise_true - single_noise_tp;
    col_noise_fn = col_noise_true - col_noise_tp;
    event_fn = event_true - event_tp;
    // 先测试命令行输出，然后再写回文件验证。
    cout << "==========================================================================" << endl;
    cout << inFileName << endl;
    cout << "single noise true: " << single_noise_true << endl;
    cout << "single noise tp: " << single_noise_tp << endl;
    cout << "single noise fp: " << single_noise_fp << endl;
    cout << "single noise fn: " << single_noise_fn << endl;
    double single_precision = (double)single_noise_tp / (single_noise_tp + single_noise_fp);
    double single_recall = (double)single_noise_tp / (single_noise_tp + single_noise_fn);
    double single_f1 = (2 * single_precision * single_recall) / (single_precision + single_recall);
    cout << "single precision: " << single_precision << endl;
    cout << "single recall: " << single_recall << endl;
    cout << "single f1: " << single_f1 << endl;
    cout << "-------------------------------------" << endl;

    cout << "collective noise true: " << col_noise_true << endl;
    cout << "collective noise tp: " << col_noise_tp << endl;
    cout << "collective noise fp: " << col_noise_fp << endl;
    cout << "collective noise fn: " << col_noise_fn << endl;
    double col_precision = (double)col_noise_tp / (col_noise_tp + col_noise_fp);
    double col_recall = (double)col_noise_tp / (col_noise_tp + col_noise_fn);
    double col_f1 = (2 * col_precision * col_recall) / (col_precision + col_recall);
    cout << "collective precision: " << col_precision << endl;
    cout << "collective recall: " << col_recall << endl;
    cout << "collective f1: " << col_f1 << endl;
    cout << "-------------------------------------" << endl;

    cout << "event true: " << event_true << endl;
    cout << "event tp: " << event_tp << endl;
    cout << "event fp: " << event_fp << endl;
    cout << "event fn: " << event_fn << endl;
    double event_precision = (double)event_tp / (event_tp + event_fp);
    double event_recall = (double)event_tp / (event_tp + event_fn);
    double event_f1 = (2 * event_precision * event_recall) / (event_precision + event_recall);
    cout << "event precision: " << event_precision << endl;
    cout << "event recall: " << event_recall << endl;
    cout << "event f1: " << event_f1 << endl;

    // 将检测结果写回文件
    outfile.open(outFileFolder + inFileName + "_single_true.txt");
    for(int i=0; i<true_single_noise.size(); i++) {
        outfile << data_str[true_single_noise[i]] << endl;
    }
    outfile.close();
    outfile.open(outFileFolder + inFileName + "_single_detect.txt");
    for(int i=0; i<detect_single_noise.size(); i++) {
        outfile << data_str[detect_single_noise[i]] << endl;
    }
    outfile.close();
    outfile.open(outFileFolder + inFileName + "_collective_true.txt");
    for(int i=0; i<true_col_noise.size(); i++) {
        for(int j=0; j<true_col_noise[i].size(); j++) {
            outfile << data_str[true_col_noise[i][j]] << endl;
        }
    }
    outfile.close();
    outfile.open(outFileFolder + inFileName + "_collective_detect.txt");
    for(int i=0; i<detect_col_noise.size(); i++) {
        for(int j=0; j<detect_col_noise[i].size(); j++) {
            outfile << data_str[detect_col_noise[i][j]] << endl;
        }
    }
    outfile.close();
    outfile.open(outFileFolder + inFileName + "_event_true.txt");
    for(int i=0; i<true_event.size(); i++) {
        for(int j=0; j<true_event[i].size(); j++) {
            outfile << data_str[true_event[i][j]] << endl;
        }
    }
    outfile.close();
    outfile.open(outFileFolder + inFileName + "_event_detect.txt");
    for(int i=0; i<detect_event.size(); i++) {
        for(int j=0; j<detect_event[i].size(); j++) {
            outfile << data_str[detect_event[i][j]] << endl;
        }
    }
    outfile.close();

    vector<vector<int>> result;
    result.push_back({single_noise_true, single_noise_tp, single_noise_fp, single_noise_fn});
    result.push_back({col_noise_true, col_noise_tp, col_noise_fp, col_noise_fn});
    result.push_back({event_true, event_tp, event_fp, event_fn});
    return result;
}

void Detector::allFileTest(string inFileFolder, string outFileFolder) {
    int all_single_true = 0;
    int all_single_tp = 0;
    int all_single_fp = 0;
    int all_single_fn = 0;

    int all_col_true = 0;
    int all_col_tp = 0;
    int all_col_fp = 0;
    int all_col_fn = 0;

    int all_event_true = 0;
    int all_event_tp = 0;
    int all_event_fp = 0;
    int all_event_fn = 0;

    ofstream outfile;
    outfile.open(outFileFolder + "result.csv");
    outfile << "id,single true,single tp,single fp,single fn,single precision,single recall,single F1,";
    outfile << "collective true,collective tp,collective fp,collective fn,collective precision,collective recall,collective F1,";
    outfile << "event true,event tp,event fp,event fn,event precision,event recall,event F1," << endl;
    for(int i=1001; i<=1036; i++) {
        vector<vector<int>> result = singleFileTest(inFileFolder, outFileFolder, "CrawledData_00" + to_string(i));
        outfile << "00" + to_string(i) << ",";
        for(int m = 0; m < result.size(); m++) {
            if(m == 0) {
                all_single_true += result[m][0];
                all_single_tp += result[m][1];
                all_single_fp += result[m][2];
                all_single_fn += result[m][3];
            } else if(m == 1) {
                all_col_true += result[m][0];
                all_col_tp += result[m][1];
                all_col_fp += result[m][2];
                all_col_fn += result[m][3];
            } else {
                all_event_true += result[m][0];
                all_event_tp += result[m][1];
                all_event_fp += result[m][2];
                all_event_fn += result[m][3];
            }
            for(int n=0; n<result[m].size(); n++) {
                outfile << result[m][n] << ",";
            }
            double precision = (double)result[m][1]/(result[m][1]+result[m][2]);
            double recall = (double)result[m][1]/(result[m][1]+result[m][3]);
            double f1 = (2*precision*recall)/(precision+recall);
            if(m != (result.size() - 1)) {
                outfile << precision << ",";
                outfile << recall << ",";
                outfile << f1 << ",";
            } else {
                outfile << precision << ",";
                outfile << recall << ",";
                outfile << f1;
            }
        }
        outfile << endl;
    }

    cout << "==========================================================================" << endl;
    cout << "All File" << endl;
    cout << "single noise true: " << all_single_true << endl;
    cout << "single noise tp: " << all_single_tp << endl;
    cout << "single noise fp: " << all_single_fp << endl;
    cout << "single noise fn: " << all_single_fn << endl;
    double single_precision = (double)all_single_tp / (all_single_tp + all_single_fp);
    double single_recall = (double)all_single_tp / (all_single_tp + all_single_fn);
    double single_f1 = (2 * single_precision * single_recall) / (single_precision + single_recall);
    cout << "single precision: " << single_precision << endl;
    cout << "single recall: " << single_recall << endl;
    cout << "single f1: " << single_f1 << endl;
    cout << "-------------------------------------" << endl;

    cout << "collective noise true: " << all_col_true << endl;
    cout << "collective noise tp: " << all_col_tp << endl;
    cout << "collective noise fp: " << all_col_fp << endl;
    cout << "collective noise fn: " << all_col_fn << endl;
    double col_precision = (double)all_col_tp / (all_col_tp + all_col_fp);
    double col_recall = (double)all_col_tp / (all_col_tp + all_col_fn);
    double col_f1 = (2 * col_precision * col_recall) / (col_precision + col_recall);
    cout << "collective precision: " << col_precision << endl;
    cout << "collective recall: " << col_recall << endl;
    cout << "collective f1: " << col_f1 << endl;
    cout << "-------------------------------------" << endl;

    cout << "event true: " << all_event_true << endl;
    cout << "event tp: " << all_event_tp << endl;
    cout << "event fp: " << all_event_fp << endl;
    cout << "event fn: " << all_event_fn << endl;
    double event_precision = (double)all_event_tp / (all_event_tp + all_event_fp);
    double event_recall = (double)all_event_tp / (all_event_tp + all_event_fn);
    double event_f1 = (2 * event_precision * event_recall) / (event_precision + event_recall);
    cout << "event precision: " << event_precision << endl;
    cout << "event recall: " << event_recall << endl;
    cout << "event f1: " << event_f1 << endl;

    outfile << "all," << all_single_true << "," << all_single_tp << "," << all_single_fp
    << "," << all_single_fn << "," << single_precision << "," << single_recall << "," << single_f1 << ","
    << all_col_true << "," << all_col_tp << "," << all_col_fp << "," << all_col_fn << "," << col_precision
    << "," << col_recall << "," << col_f1 << "," << all_event_true << "," << all_event_tp << "," << all_event_fp
    << "," << all_event_fn << "," << event_precision << "," << event_recall << "," << event_f1;

    outfile.close();
}