model2rrn = {"16srna": "rrnL", "12srna": "rrnS"}

rrn2symbol = {"rrnS": "MT-RNR1", "rrnL": "MT-RNR2"}

rrn2product = {
    "12srna": "12S rRNA",
    "rrn12": "12S rRNA",
    "16srna": "16S rRNA",
    "rrn16": "16S rRNA",
    "rrnS": "12S rRNA",
    "rrnL": "16S rRNA",
}
for k,v in model2rrn.items():
    rrn2symbol[k] = rrn2symbol[v]


trn2symbol = {
    "trnA": "MT-TA",
    "trnC": "MT-TC",
    "trnD": "MT-TD",
    "trnE": "MT-TE",
    "trnF": "MT-TF",
    "trnG": "MT-TG",
    "trnH": "MT-TH",
    "trnI": "MT-TI",
    "trnK": "MT-TK",
    "trnL2": "MT-TL2",
    "trnL1": "MT-TL1",
    "trnM": "MT-TM",
    "trnN": "MT-TN",
    "trnP": "MT-TP",
    "trnQ": "MT-TQ",
    "trnR": "MT-TR",
    "trnS1": "MT-TS1",
    "trnS2": "MT-TS2",
    "trnT": "MT-TT",
    "trnV": "MT-TV",
    "trnW": "MT-TW",
    "trnY": "MT-TY",
}

trn2product = {
    "trnA": "tRNA-Ala",
    "trnC": "tRNA-Cys",
    "trnD": "tRNA-Asp",
    "trnE": "tRNA-Glu",
    "trnF": "tRNA-Phe",
    "trnG": "tRNA-Gly",
    "trnH": "tRNA-His",
    "trnI": "tRNA-Ile",
    "trnK": "tRNA-Lys",
    "trnL2": "tRNA-Leu",
    "trnL1": "tRNA-Leu",
    "trnM": "tRNA-Met",
    "trnN": "tRNA-Asn",
    "trnP": "tRNA-Pro",
    "trnQ": "tRNA-Gln",
    "trnR": "tRNA-Arg",
    "trnS1": "tRNA-Ser",
    "trnS2": "tRNA-Ser",
    "trnT": "tRNA-Thr",
    "trnV": "tRNA-Val",
    "trnW": "tRNA-Trp",
    "trnY": "tRNA-Tyr",
}

cds2symbol = {
    "ATP6": "MT-ATP6",
    "ATP8": "MT-ATP8",
    "COX1": "MT-CO1",
    "COX2": "MT-CO2",
    "COX3": "MT-CO3",
    "CYTB": "MT-CYTB",
    "ND1": "MT-ND1",
    "ND2": "MT-ND2",
    "ND3": "MT-ND3",
    "ND4": "MT-ND4",
    "ND4L": "MT-ND4L",
    "ND5": "MT-ND5",
    "ND6": "MT-ND6",
}

cds2product = {
    "ND1": "NADH dehydrogenase subunit 1",
    "nad1": "NADH dehydrogenase subunit 1",
    "ND2": "NADH dehydrogenase subunit 2",
    "nad2": "NADH dehydrogenase subunit 2",
    "ND3": "NADH dehydrogenase subunit 3",
    "nad3": "NADH dehydrogenase subunit 3",
    "ND4": "NADH dehydrogenase subunit 4",
    "nad4": "NADH dehydrogenase subunit 4",
    "ND4L": "NADH dehydrogenase subunit 4L",
    "nad4L": "NADH dehydrogenase subunit 4L",
    "ND5": "NADH dehydrogenase subunit 5",
    "nad5": "NADH dehydrogenase subunit 5",
    "ND6": "NADH dehydrogenase subunit 6",
    "nad6": "NADH dehydrogenase subunit 6",
    "COX1": "cytochrome c oxidase subunit I",
    "cox1": "cytochrome c oxidase subunit I",
    "COX2": "cytochrome c oxidase subunit II",
    "cox2": "cytochrome c oxidase subunit II",
    "COX3": "cytochrome c oxidase subunit III",
    "cox3": "cytochrome c oxidase subunit III",
    "ATP6": "ATP synthase F0 subunit 6",
    "atp6": "ATP synthase F0 subunit 6",
    "ATP8": "ATP synthase F0 subunit 8",
    "atp8": "ATP synthase F0 subunit 8",
    "CYTB": "cytochrome b",
    "cob": "cytochrome b",
}
def main():
    print("""
from organelle_svg.og_colors import GeneColor, mregex
from organelle_svg import og_colors as og""")
    print("EmmaColors = [")
    for k,v in rrn2symbol.items():
        # print(f"{v} -> {rrn2product[k]}")
        print('GeneColor(')
        print('  type=mregex("rRNA"),')
        print(f'  pattern=mregex("^{v}$"),')
        print('  color_tuple=og.VIOLET,')
        print(f'  fullname="{rrn2product[k]}",')
        print('),')

    for k,v in cds2symbol.items():
        print('GeneColor(')
        print('  type=mregex("CDS"),')
        print(f'  pattern=mregex("^{v}$"),')
        print('  color_tuple=og.VIOLET,')
        print(f'  fullname="{cds2product[k]}",')
        print('),')

    for k,v in trn2symbol.items():
        print('GeneColor(')
        print('  type=mregex("tRNA"),')
        print(f'  pattern=mregex("^{v}$"),')
        print('  color_tuple=og.VIOLET,')
        print(f'  fullname="{trn2product[k]}",')
        print('),')
    print("]")

main()
