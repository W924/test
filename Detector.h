#pragma once
#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>

class Detector {
public:

    // 基本功能函数
    static double distance(std::vector<double> &point1, std::vector<double> &point2, int loc1, int loc2);
    static std::vector<double> absoluteSecondDerivative(std::vector<std::vector<double>> &data);
    std::vector<int> candidate(std::vector<std::vector<double>> &data);
    std::pair<std::map<int, double>, std::vector<std::pair<int, double>>> ascendingDistance(std::vector<std::vector<double>> &data, int loc);
    std::set<int> knn(std::vector<std::vector<double>> &data, std::vector<std::pair<int, double>> &id_distance, int k);
    std::set<int> snn(std::vector<std::vector<double>> &data, int loc);
//    std::vector<double> feather_extract(std::vector<std::vector<double>> &data, int loc);

    // 单文件测试，测试的结果输出到output文件夹中
    std::vector<std::vector<int>> singleFileTest(std::string inFileFolder, std::string outFileFolder, std::string inFileName);

    // 全文件测试
    void allFileTest(std::string inFileFolder, std::string outFileFolder);


};
