import os
import numpy as np
def read_lattice_vector(filename='final.config'):
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        exit()
    with open(filename, 'r') as file:
        lines = file.readlines()
    read=False
    lattice_vector=[]
    for line in lines:
        if 'lattice' in line.lower():
            read=True
            num=0
            continue
        if read and num<3:
            parts = line.split()
            lattice_vector.append([float(parts[0]), float(parts[1]), float(parts[2])])
            num+=1
            if num==3:
                break
    return lattice_vector

def get_metric_tensor(lattice_vectors):
    return np.dot(lattice_vectors, lattice_vectors.T)

# 计算两点之间的最小距离（考虑周期性）
def min_distance(frac1, frac2, lattice_vectors):
    min_d = float('inf')
    G = get_metric_tensor(np.array(lattice_vectors))
    
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                # 计算分数坐标差
                delta_frac = np.array([
                    frac2[0] + dx - frac1[0],
                    frac2[1] + dy - frac1[1],
                    frac2[2] + dz - frac1[2]
                ])

                delta_frac = delta_frac - np.round(delta_frac)               
                distance_sq = np.dot(delta_frac, np.dot(G, delta_frac))
                if distance_sq < min_d:
                    min_d = distance_sq

    return np.sqrt(min_d)

def read_vatom_v2(defect_center, lattice_vectors):
    with open('OUT.VATOM', 'r') as file:
        lines = file.readlines()
        atomic_list = []
        atomic_coordinates_frac = []
        vatom = []
        
        # 读取原子信息（保持原始分数坐标）
        for i, line in enumerate(lines[1:], 1):
            parts = line.split()
            atomic_list.append(i)
            atomic_coordinates_frac.append([float(p) for p in parts[1:4]])
            vatom.append(parts[4:])

    distances = []
    for frac_coord in atomic_coordinates_frac:
        # 使用修改后的距离计算函数
        min_d = min_distance(np.array(defect_center), np.array(frac_coord), lattice_vectors)
        distances.append(min_d)
        
    # 找到最大距离对应的原子
    max_index = distances.index(max(distances))

    print(f'{atomic_list[max_index]+1:4} '
        f'{" ".join(f"{coord:.5f}".rjust(8) for coord in atomic_coordinates_frac[max_index])} '
        f'{" ".join(str(coord).rjust(8) for coord in vatom[max_index])}')


def read_vatom_v0(defect_center,lattice_vectors):
    lattice_vectors = np.array(lattice_vectors)
    with open('OUT.VATOM', 'r') as file:
        lines = file.readlines()
        distances= []
        atomic_list=[]
        atomic_coordinates=[]
        vatom=[]
        i=1
        for line in lines[1:]:
            parts = line.split()
            atomic_list.append(i)
            atomic_coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])
            vatom.append(parts[4:])

            far_point = (np.array(defect_center) + 0.5) % 1.0
            delta_frac_adjusted = np.array(atomic_coordinates[i-1]) - far_point
            delta_frac_adjusted -= np.round(delta_frac_adjusted)
            A = np.array(lattice_vectors).T
            delta_cart = np.dot(A, delta_frac_adjusted)
            distance=np.linalg.norm(delta_cart)
            distances.append(distance)
            i+=1

    number = distances.index(min(distances))
    print(f'{atomic_list[number]+1:4} ' 
          f'{" ".join(f"{coord:.5f}".rjust(8) for coord in atomic_coordinates[number])} '
          f'{" ".join(str(coord).rjust(8) for coord in vatom[number])}')

def read_vatom_v1(defect_center,lattice_vectors):
    with open('OUT.VATOM', 'r') as file:
        lines = file.readlines()
        distances= []
        atomic_list=[]
        atomic_coordinates=[]
        vatom=[]
        i=1
        for line in lines[1:]:
            parts = line.split()
            atomic_list.append(i)
            atomic_coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])
            vatom.append(parts[4:])
            i+=1
   
        for j in range(len(atomic_coordinates)):
            # 计算分数坐标差并调整到最短距离
            delta_frac = np.array(atomic_coordinates[j]) - np.array(defect_center)
            delta_frac -= np.round(delta_frac)
            # 转换回真实坐标并计算距离
            real_delta = delta_frac @ lattice_vectors
            distance = np.linalg.norm(real_delta)
            distances.append(distance) 
    number = distances.index(max(distances))   
    print(f'{atomic_list[number]+1:4} ' 
          f'{" ".join(f"{coord:.5f}".rjust(8) for coord in atomic_coordinates[number])} '
          f'{" ".join(str(coord).rjust(8) for coord in vatom[number])}')

defect_center=[]
print("the farest atom from the defect coordinates, pot:")
input_coords = input().strip()
defect_center = [float(coord) for coord in input_coords.split()]
lattice_vectors = read_lattice_vector(filename='final.config')
#任选其一
read_vatom_v0(defect_center, lattice_vectors)
#read_vatom_v1(defect_center,lattice_vectors)
#read_vatom_v2(defect_center,lattice_vectors)
