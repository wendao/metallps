def find_nearby_residues(pdb_file, threshold=6):
    # 读取文件
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    # 解析金属离子坐标（最后一行）
    last_line = lines[-1].strip()
    if not last_line.startswith('HETATM'):
        raise ValueError("最后一行的记录不是HETATM，无法获取金属离子坐标。")
    try:
        x = float(last_line[30:38].strip())
        y = float(last_line[38:46].strip())
        z = float(last_line[46:54].strip())
    except:
        raise ValueError("无法解析金属离子的坐标。")
    
    residues = set()
    
    # 遍历所有行处理ATOM记录
    for line in lines:
        line = line.strip()
        if line.startswith('ATOM'):
            chain = line[21].strip()
            if chain not in ('A', 'B'):
                continue
            residue_num = line[22:26].strip()  # 残基编号，保留插入码
            # 解析原子坐标
            try:
                x_atom = float(line[30:38].strip())
                y_atom = float(line[38:46].strip())
                z_atom = float(line[46:54].strip())
            except:
                continue  # 坐标解析失败则跳过
            
            # 计算距离
            distance = ((x_atom - x)**2 + (y_atom - y)**2 + (z_atom - z)**2) ** 0.5
            if distance < threshold:
                residue_id = f"{chain}{residue_num}"
                residues.add(residue_id)
    
    # 排序函数：按链和残基编号的自然顺序排序
    def sort_key(res):
        chain_part = res[0]
        num_part = res[1:]
        # 分离数字和字母部分（例如处理'100A'）
        digits = ''
        chars = ''
        for c in num_part:
            if c.isdigit():
                digits += c
            else:
                chars += c
        return (chain_part, int(digits) if digits else 0, chars)
    
    # 转换为列表并排序
    sorted_residues = sorted(residues, key=sort_key)
    
    return sorted_residues

# 示例用法
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("请提供PDB文件路径。")
        sys.exit(1)
    pdb_path = sys.argv[1]
    threshold = 6  # 可调整阈值
    nearby = find_nearby_residues(pdb_path, threshold)
    print(" ".join(nearby))
