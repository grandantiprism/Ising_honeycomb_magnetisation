import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re

def get_binder_data(folder_path):
    """特定のフォルダ内の全ファイルを解析してbetaとU4のリストを返す"""
    betas = []
    u4_values = []
    
    # ファイル名形式: 32x32_0.441.txt から beta を抽出
    file_list = glob.glob(os.path.join(folder_path, "*.txt"))
    
    def extract_beta(filename):
        match = re.search(r"_(\d+\.\d+)\.txt", filename)
        return float(match.group(1)) if match else None

    # betaの昇順にソート
    file_list = [f for f in file_list if extract_beta(f) is not None]
    file_list.sort(key=extract_beta)

    for file_path in file_list:
        beta = extract_beta(file_path)
        try:
            m_data = np.loadtxt(file_path)
            # 磁化が空でないかチェック
            if m_data.size == 0: continue
            
            m2_avg = np.mean(m_data**2)
            m4_avg = np.mean(m_data**4)
            # u4 = 1.0 - (m4_avg / (3.0 * (m2_avg**2)))
            u4 = m4_avg / (m2_avg**2)
            
            betas.append(beta)
            u4_values.append(u4)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    return np.array(betas), np.array(u4_values)

def plot_all_sizes(root_dir):
    plt.figure(figsize=(10, 7))
    
    # 1. 使用するマーカーのリストを定義
    # 'o':円, 's':正方形, '^':三角形, 'v':逆三角形, 'D':菱形, 'p':五角形, '*':星型
    markers = ['o', 's', '^', 'v', 'D', 'p', '*', 'h']
    
    folders = [d for d in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, d)) and 'x' in d]
    folders.sort(key=lambda x: int(x.split('x')[0]))

    if not folders:
        print(f"No size folders found in {root_dir}")
        return

    # 2. enumerateを使ってインデックス i を取得し、マーカーを切り替える
    for i, folder in enumerate(folders):
        path = os.path.join(root_dir, folder)
        
        if folder.split("x")[0] not in [str(l) for l in l_list]:
            continue

        print(f"Processing folder: {folder}...")
        betas, u4 = get_binder_data(path)

        if len(betas) > 0:
            # markers[i % len(markers)] でリストを使い回す
            plt.plot(betas, u4, 
                     marker=markers[i % len(markers)], 
                     linestyle='None',  # 線をなしにする別の書き方
                     label=f'L = {folder.split("x")[0]}', 
                     markersize=6, 
                     alpha=0.8)

    # 理論的な臨界点
    beta_c = 0.5 * np.log(2.0 + np.sqrt(3))
    plt.axvline(x=beta_c, color='black', linestyle='--', alpha=0.6, label=f'Theory βc ≈ {beta_c:.4f}')

    plt.title('Binder Cumulant Crossing', fontsize=14)
    plt.xlabel(r'Inverse Temperature $\beta$', fontsize=12)
    plt.ylabel(r'Binder Ratio $U_4$', fontsize=12)
    plt.grid(True, which='both', linestyle=':', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("binder_plot.png", dpi=300)
    plt.show()

# --- 実行 ---
# C++コードが生成したディレクトリがカレントディレクトリにある場合は "."
# "output" というフォルダにまとめている場合は "output" を指定してください

l_list = [16, 24, 32, 48, 64, 96, 128]

root_directory = "." 
plot_all_sizes(root_directory)