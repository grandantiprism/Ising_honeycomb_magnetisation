#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

// シミュレーションパラメータの設定
const int L = 16;            // 格子サイズ (L x L)
const int N = L * L;        // 全スピン数
const int MCS = 1000000;      // 各温度でのモンテカルロステップ数
const int THERM = (L < 64) ? 1000 : L * 20;     // 熱平衡化のための捨てステップ数

const double beta_min = 0.6;
const double beta_max = 0.7;
const int num_beta = 20;

struct Ising2D {
    int L;
    vector<int> spins;
    mt19937 gen;
    uniform_int_distribution<int> dist_site;
    uniform_real_distribution<double> dist_prob;

    Ising2D(int l, int seed) : L(l), spins(l * l, 1), gen(seed),
                               dist_site(0, l * l - 1), dist_prob(0.0, 1.0) {}

    // 2次元インデックスを1次元に変換 (周期境界条件)
    int get_idx(int x, int y) {
        return ((x + L) % L) * L + ((y + L) % L);
    }

    // 全磁化の計算
    double calc_magnetization() {
        double m = 0;
        for (int s : spins) m += s;
        return m / (L * L);
    }

    // Wolffアルゴリズムによる1ステップの更新
    void wolff_step(double beta) {
        double p_add = 1.0 - exp(-2.0 * beta);
        int root = dist_site(gen);
        int old_spin = spins[root];
        
        vector<int> cluster_stack;
        cluster_stack.push_back(root);
        spins[root] *= -1;

        while (!cluster_stack.empty()) {
            int curr = cluster_stack.back();
            cluster_stack.pop_back();

            int x = curr / L;
            int y = curr % L;

            // 隣接3サイトの決定（レンガ格子構造）
            int neighbors[3];
            neighbors[0] = get_idx(x, y + 1); // 上
            neighbors[1] = get_idx(x, y - 1); // 下
            
            if ((x + y) % 2 == 0) {
                neighbors[2] = get_idx(x + 1, y); // 右
            } else {
                neighbors[2] = get_idx(x - 1, y); // 左
            }

            for (int next : neighbors) {
                if (spins[next] == old_spin) {
                    if (dist_prob(gen) < p_add) {
                        spins[next] *= -1;
                        cluster_stack.push_back(next);
                    }
                }
            }
        }
    }
};

int main() {
    string dir_name = "output/" + to_string(L) + "x" + to_string(L);
    if (!fs::exists(dir_name)) {
        fs::create_directories(dir_name);
    }
    
    auto start = chrono::high_resolution_clock::now();
    
    // 逆温度の設定
    double beta_step = (beta_max - beta_min) / num_beta;

    for (int i = 0; i < num_beta + 1; ++i) {
        double beta = beta_min + i * beta_step;
        Ising2D model(L, 12345); // 固定シード

        // ファイル出力の準備
        stringstream ss;
        ss << dir_name << "/" << L << "x" << L << "_" << fixed << setprecision(5) << beta << ".txt";
        ofstream ofs(ss.str(), ios::app);
        
        ofs << fixed << setprecision(20);

        cout << "Simulating beta = " << beta << "..." << endl;

        // 熱平衡化
        for (int t = 0; t < THERM; ++t) {
            model.wolff_step(beta);
        }

        // 本測定
        for (int t = 0; t < MCS; ++t) {
            model.wolff_step(beta);
            double m = model.calc_magnetization();
            ofs << m << "\n"; // 各ステップの磁化を出力
        }
        
        ofs.close();
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    // ログファイルへの追記
    ofstream log_file("log.txt", ios::app);
    if (log_file) {
        // 現在時刻の取得（オプション）
        auto now = chrono::system_clock::to_time_t(chrono::system_clock::now());
        
        log_file << "--- Simulation Log ---" << endl;
        log_file << "Date: " << ctime(&now); // 実行日時
        log_file << "L: " << L << ", MCS: " << MCS << ", THERM: " << THERM << endl;
        log_file << "Num_beta: " << num_beta << " (" << beta_min << " to " << beta_max << ")" << endl;
        log_file << "Elapsed time: " << fixed << setprecision(2) << elapsed.count() << " seconds" << endl;
        log_file << "-----------------------" << endl << endl;
        log_file.close();
    }

    cout << "Simulation completed in " << elapsed.count() << " seconds." << endl;
    return 0;
}
