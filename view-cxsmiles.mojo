from std.python import Python

def main() raises:
    var Chem = Python.import_module("rdkit.Chem")
    var Draw = Python.import_module("rdkit.Chem.Draw")
    var rdDepictor = Python.import_module("rdkit.Chem.rdDepictor")
    rdDepictor.SetPreferCoordGen(True)
    var cxsmiles = "O1([H]O([H]O([H]1)[H])[H])[H] |H:5.5,1.1,3.3|"

    # 水素原子を削除しないようにパーサー設定を変更
    var params = Chem.SmilesParserParams()
    params.removeHs = False

    var mol = Chem.MolFromSmiles(cxsmiles, params)

    if mol:
        var img = Draw.MolToImage(mol)
        img.save("output.png")
        print("output.png に保存しました")
    else:
        print("CXSMILESのパースに失敗しました")
