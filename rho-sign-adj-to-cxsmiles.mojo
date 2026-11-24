from std.collections import List

# 1. 隣接行列から共有結合リストと水素結合ペア(元インデックス)を抽出
def extract_bonds(
    atoms: List[String], adj_matrix: List[List[Int]]
) -> Tuple[List[List[Int]], List[Tuple[Int, Int]]]:
    var n = len(atoms)
    var cov_adj = List[List[Int]]()
    for _ in range(n):
        cov_adj.append(List[Int]())

    var raw_hbonds = List[Tuple[Int, Int]]()

    for i in range(n):
        for j in range(i + 1, n):
            var val = adj_matrix[i][j]
            if val == -1 or val == 1:
                cov_adj[i].append(j)
                cov_adj[j].append(i)
            if val == 1:
                # 元のH -> 重原子の向きに整えるロジックをそのまま活用
                if atoms[i] == "H" and atoms[j] != "H":
                    raw_hbonds.append((i, j))
                elif atoms[j] == "H" and atoms[i] != "H":
                    raw_hbonds.append((j, i))
                else:
                    raw_hbonds.append((i, j))

    return (cov_adj^, raw_hbonds^)


# 2. 走査計画: DFSで書き込み順(old_to_new)・親・環閉鎖(非木辺)を決める
#    元のコードは木辺しか出力せず、環状クラスターでは辺が1本欠落していた。
#    ここで非木辺を検出して両端の原子に環番号を割り当てておく。
def plan_dfs(
    cov_adj: List[List[Int]],
    u: Int,
    par: Int,
    mut visited: List[Bool],
    mut old_to_new: List[Int],
    mut parent: List[Int],
    mut depth: List[Int],
    mut digits_at: List[List[Int]],
    mut ring_u: List[Int],
    mut ring_v: List[Int],
    mut counter: Int,
):
    visited[u] = True
    old_to_new[u] = counter
    counter += 1
    parent[u] = par

    for v in cov_adj[u]:
        if not visited[v]:
            depth[v] = depth[u] + 1
            plan_dfs(
                cov_adj.copy(),
                v,
                u,
                visited,
                old_to_new,
                parent,
                depth,
                digits_at,
                ring_u,
                ring_v,
                counter,
            )
        elif v != par and depth[v] < depth[u]:
            # 後方辺 (u, v) = 環閉鎖。環番号は ring_u/ring_v の格納順 (1始まり)
            ring_u.append(u)
            ring_v.append(v)
            var rn = len(ring_u)
            digits_at[u].append(rn)
            digits_at[v].append(rn)


def plan_traversal(
    cov_adj: List[List[Int]],
) -> Tuple[List[Int], List[Int], List[List[Int]], List[Int], List[Int]]:
    var n = len(cov_adj)
    var visited = List[Bool]()
    var old_to_new = List[Int]()
    var parent = List[Int]()
    var depth = List[Int]()
    var digits_at = List[List[Int]]()
    for _ in range(n):
        visited.append(False)
        old_to_new.append(-1)
        parent.append(-1)
        depth.append(0)
        digits_at.append(List[Int]())

    var ring_u = List[Int]()
    var ring_v = List[Int]()
    var counter = 0

    for i in range(n):
        if not visited[i]:
            plan_dfs(
                cov_adj.copy(),
                i,
                -1,
                visited,
                old_to_new,
                parent,
                depth,
                digits_at,
                ring_u,
                ring_v,
                counter,
            )

    return (old_to_new^, parent^, digits_at^, ring_u^, ring_v^)


# 3. SMILES文字列の構築 (木辺 + 環閉鎖数字)
def ring_digit_str(rn: Int) -> String:
    if rn < 10:
        return String(rn)
    return "%" + String(rn)


def write_dfs(
    atoms: List[String],
    children: List[List[Int]],
    digits_at: List[List[Int]],
    u: Int,
    mut out_str: String,
):
    var sym = atoms[u]
    out_str += "[" + sym + "]" if sym == "H" else sym
    for rn in digits_at[u]:
        out_str += ring_digit_str(rn)

    # 修正1: var ではなく ref を使って参照を受け取る（コピーが発生しない）
    ref ch = children[u]
    for k in range(len(ch)):
        var v = ch[k]
        var is_branch = k < len(ch) - 1
        if is_branch:
            out_str += "("
        # 修正2: 引数は読み取り専用で渡せるため .copy() を削除して高速化
        write_dfs(atoms, children, digits_at, v, out_str)
        if is_branch:
            out_str += ")"



def build_smiles(
    atoms: List[String],
    old_to_new: List[Int],
    parent: List[Int],
    digits_at: List[List[Int]],
) -> String:
    var n = len(atoms)

    var new_to_old = List[Int]()
    for _ in range(n):
        new_to_old.append(-1)
    for old in range(n):
        new_to_old[old_to_new[old]] = old

    # 子リスト (書き込み順)
    var children = List[List[Int]]()
    for _ in range(n):
        children.append(List[Int]())
    for ni in range(n):
        var old = new_to_old[ni]
        if parent[old] != -1:
            children[parent[old]].append(old)

    var smiles_parts = List[String]()
    for i in range(n):
        if parent[i] == -1:
            var mol_str = String("")
            write_dfs(atoms.copy(), children.copy(), digits_at.copy(), i, mol_str)
            smiles_parts.append(mol_str)

    var full_smiles = String("")
    for idx in range(len(smiles_parts)):
        if idx > 0:
            full_smiles += "."
        full_smiles += smiles_parts[idx]

    return full_smiles


# 4. 結合インデックスの算出
#    CXSMILESの |H:| 拡張は「原子インデックス.結合インデックス」形式で、
#    結合インデックスはSMILES文字列に現れる順 (木辺は子の出現順、
#    環閉鎖は閉じる数字の出現順) を指す。RDKitもこの順で解釈する。
def compute_bond_indices(
    old_to_new: List[Int],
    parent: List[Int],
    digits_at: List[List[Int]],
    ring_u: List[Int],
    ring_v: List[Int],
) -> Tuple[List[Int], List[Int]]:
    var n = len(old_to_new)

    var new_to_old = List[Int]()
    for _ in range(n):
        new_to_old.append(-1)
    for old in range(n):
        new_to_old[old_to_new[old]] = old

    var tree_bond_idx = List[Int]()
    for _ in range(n):
        tree_bond_idx.append(-1)
    var ring_bond_idx = List[Int]()
    for _ in range(len(ring_u)):
        ring_bond_idx.append(-1)

    var counter = 0
    for ni in range(n):
        var old = new_to_old[ni]
        if parent[old] != -1:
            tree_bond_idx[old] = counter
            counter += 1
        for rn in digits_at[old]:
            var k = rn - 1
            var other = ring_v[k] if ring_u[k] == old else ring_u[k]
            if old_to_new[other] < ni:
                ring_bond_idx[k] = counter
                counter += 1

    return (tree_bond_idx^, ring_bond_idx^)


# 5. CXSMILESの水素結合拡張ブロックを生成
#    各ペアは (水素の新インデックス).(その水素結合の結合インデックス)
def format_hbonds_extension(
    raw_hbonds: List[Tuple[Int, Int]],
    old_to_new: List[Int],
    parent: List[Int],
    tree_bond_idx: List[Int],
    ring_u: List[Int],
    ring_v: List[Int],
    ring_bond_idx: List[Int],
) -> String:
    if len(raw_hbonds) == 0:
        return String("")

    var hb_str = String("")
    for k in range(len(raw_hbonds)):
        if k > 0:
            hb_str += ","
        var pair = raw_hbonds[k]
        var h = pair[0]
        var heavy = pair[1]
        var bond_idx = -1
        if parent[h] == heavy:
            bond_idx = tree_bond_idx[h]
        elif parent[heavy] == h:
            bond_idx = tree_bond_idx[heavy]
        else:
            for e in range(len(ring_u)):
                if (ring_u[e] == h and ring_v[e] == heavy) or (
                    ring_u[e] == heavy and ring_v[e] == h
                ):
                    bond_idx = ring_bond_idx[e]
                    break
        hb_str += String(old_to_new[h]) + "." + String(bond_idx)

    return " |H:" + hb_str + "|"


# エントリポイント（元のインターフェースを維持）
def adj_matrix_to_cxsmiles(
    atoms: List[String], adj_matrix: List[List[Int]]
) -> String:
    var bonds = extract_bonds(atoms, adj_matrix)
    var cov_adj = bonds[0].copy()
    var raw_hbonds = bonds[1].copy()

    var plan = plan_traversal(cov_adj)
    var old_to_new = plan[0].copy()
    var parent = plan[1].copy()
    var digits_at = plan[2].copy()
    var ring_u = plan[3].copy()
    var ring_v = plan[4].copy()

    var smiles_part = build_smiles(atoms, old_to_new, parent, digits_at)

    var bidx = compute_bond_indices(old_to_new, parent, digits_at, ring_u, ring_v)
    var tree_bond_idx = bidx[0].copy()
    var ring_bond_idx = bidx[1].copy()

    var ext_part = format_hbonds_extension(
        raw_hbonds,
        old_to_new,
        parent,
        tree_bond_idx,
        ring_u,
        ring_v,
        ring_bond_idx,
    )
    return smiles_part + ext_part



def main():
    # ==========================================
    # 使用例：水の3量体クラスター (全9原子)
    # ==========================================
    var atoms: List[String] = ["O", "H", "H", "O", "H", "H", "O", "H", "H"]
    var n = len(atoms)

    # 9x9 のゼロ行列を初期化
    var adj = List[List[Int]]()
    for _ in range(n):
        var row = List[Int]()
        for _ in range(n):
            row.append(0)
        adj.append(row^)

    # 共有結合 (-1) の設定
    adj[0][1] = -1; adj[1][0] = -1
    adj[0][2] = -1; adj[2][0] = -1
    adj[3][4] = -1; adj[4][3] = -1
    adj[3][5] = -1; adj[5][3] = -1
    adj[6][7] = -1; adj[7][6] = -1
    adj[6][8] = -1; adj[8][6] = -1

    # 水素結合 (1) の設定 (環状クラスター)
    adj[1][3] = 1; adj[3][1] = 1
    adj[4][6] = 1; adj[6][4] = 1
    adj[7][0] = 1; adj[0][7] = 1

    # CXSMILESの生成と出力
    var cxsmiles = adj_matrix_to_cxsmiles(atoms, adj)
    print(cxsmiles)
