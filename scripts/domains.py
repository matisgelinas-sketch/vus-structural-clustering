"""
Curated domain/region annotations. TP53/BRCA1/PTEN are hand-curated
(predate the automated pipeline); every other gene is generated
automatically from UniProt features by scripts/phase2b/generate_domains.py
-- do not hand-edit those blocks, they get overwritten on the next run.
Ranges are UniProt canonical numbering (1-based, inclusive). Order
matters: first matching range wins (most specific/structurally
meaningful first).
"""

TP53_DOMAINS = [
    (102, 292, "DNA-binding domain"),
    (325, 356, "Oligomerization (tetramerization) domain"),
    (1, 44, "Transactivation domain (acidic)"),
    (50, 96, "Disordered (TAD2/proline-rich linker)"),
    (282, 325, "Disordered (DBD-tetramerization linker)"),
    (351, 393, "Disordered (C-terminal regulatory)"),
]

BRCA1_DOMAINS = [
    (24, 65, "RING-type zinc finger"),
    (1642, 1736, "BRCT domain 1"),
    (1756, 1855, "BRCT domain 2"),
    (230, 270, "Disordered"),
    (306, 338, "Disordered"),
    (534, 570, "Disordered"),
    (654, 709, "Disordered"),
    (1181, 1216, "Disordered"),
    (1322, 1387, "Disordered"),
    (1440, 1505, "Disordered"),
    (1565, 1596, "Disordered"),
]

PTEN_DOMAINS = [
    (14, 185, "Phosphatase (tensin-type) domain"),
    (190, 350, "C2 (tensin-type) domain"),
    (352, 403, "Disordered (C-terminal tail)"),
]

ABCA4_DOMAINS = [
    (929, 1160, "ABC transporter 1"),
    (1938, 2170, "ABC transporter 2"),
    (891, 911, "Disordered"),
    (1284, 1345, "Disordered"),
    (2244, 2249, "Essential for ATP binding and ATPase activity"),
]

ATP7A_DOMAINS = [
    (8, 74, "HMA 1"),
    (85, 151, "HMA 2"),
    (171, 237, "HMA 3"),
    (277, 343, "HMA 4"),
    (377, 443, "HMA 5"),
    (488, 554, "HMA 6"),
    (564, 630, "HMA 7"),
    (1486, 1500, "PDZD11-binding"),
]

COL1A1_DOMAINS = [
    (38, 96, "VWFC"),
    (1229, 1464, "Fibrillar collagen NC1"),
    (98, 1214, "Disordered"),
    (162, 178, "Nonhelical region (N-terminal)"),
    (179, 1192, "Triple-helical region"),
    (1193, 1218, "Nonhelical region (C-terminal)"),
]

COL2A1_DOMAINS = [
    (32, 90, "VWFC"),
    (1253, 1487, "Fibrillar collagen NC1"),
    (97, 1237, "Disordered"),
    (201, 1214, "Triple-helical region"),
    (1215, 1241, "Nonhelical region (C-terminal)"),
]

COL5A1_DOMAINS = [
    (72, 244, "Laminin G-like"),
    (1609, 1837, "Fibrillar collagen NC1"),
    (231, 443, "Nonhelical region"),
    (242, 269, "Disordered"),
    (281, 457, "Disordered"),
    (444, 558, "Interrupted collagenous region"),
    (470, 520, "Disordered"),
    (526, 545, "Disordered"),
    (559, 1574, "Disordered"),
    (559, 1570, "Triple-helical region"),
    (1571, 1605, "Nonhelical region"),
]

ENG_DOMAINS = [
    (363, 533, "ZP"),
    (26, 337, "Required for interaction with GDF2"),
    (26, 46, "OR1, N-terminal part"),
    (47, 199, "OR2"),
    (200, 330, "OR1, C-terminal part"),
    (270, 282, "Essential for interaction with GDF2"),
    (626, 658, "Disordered"),
]

HNF1A_DOMAINS = [
    (1, 32, "HNF-p1"),
    (87, 182, "POU-specific atypical"),
    (199, 279, "Homeobox; HNF1-type"),
    (1, 31, "Dimerization"),
    (40, 81, "Disordered"),
    (130, 132, "Interaction with DNA"),
    (143, 149, "Interaction with DNA"),
    (155, 158, "Interaction with DNA"),
    (183, 205, "Disordered"),
    (203, 206, "Interaction with DNA"),
    (263, 265, "Interaction with DNA"),
    (270, 273, "Interaction with DNA"),
    (283, 358, "Disordered"),
    (545, 567, "Disordered"),
]

LDLR_DOMAINS = [
    (25, 65, "LDL-receptor class A 1"),
    (66, 106, "LDL-receptor class A 2"),
    (107, 145, "LDL-receptor class A 3"),
    (146, 186, "LDL-receptor class A 4"),
    (195, 233, "LDL-receptor class A 5"),
    (234, 272, "LDL-receptor class A 6"),
    (274, 313, "LDL-receptor class A 7"),
    (314, 353, "EGF-like 1"),
    (354, 393, "EGF-like 2; calcium-binding"),
    (397, 438, "LDL-receptor class B 1"),
    (439, 485, "LDL-receptor class B 2"),
    (486, 528, "LDL-receptor class B 3"),
    (529, 572, "LDL-receptor class B 4"),
    (573, 615, "LDL-receptor class B 5"),
    (616, 658, "LDL-receptor class B 6"),
    (663, 712, "EGF-like 3"),
    (146, 233, "Binding to Getah virus E1-E2 spike glycoproteins"),
    (721, 768, "Clustered O-linked oligosaccharides"),
    (734, 755, "Disordered"),
    (811, 860, "Required for MYLIP-triggered down-regulation of LDLR"),
]

MLH1_DOMAINS = [
    (355, 378, "Disordered"),
    (400, 491, "Disordered"),
    (410, 650, "Interaction with EXO1"),
]

MSH2_DOMAINS = [
    (601, 671, "Interaction with EXO1"),
]

MSH6_DOMAINS = [
    (92, 154, "PWWP"),
    (1, 84, "Disordered"),
    (195, 362, "Disordered"),
]

MYBPC3_DOMAINS = [
    (153, 256, "Ig-like C2-type 1"),
    (362, 452, "Ig-like C2-type 2"),
    (453, 543, "Ig-like C2-type 3"),
    (544, 633, "Ig-like C2-type 4"),
    (645, 771, "Ig-like C2-type 5"),
    (774, 870, "Fibronectin type-III 1"),
    (872, 967, "Fibronectin type-III 2"),
    (971, 1065, "Ig-like C2-type 6"),
    (1068, 1163, "Fibronectin type-III 3"),
    (1181, 1274, "Ig-like C2-type 7"),
    (107, 153, "Disordered"),
]

MYO7A_DOMAINS = [
    (65, 741, "Myosin motor"),
    (745, 765, "IQ 1"),
    (768, 788, "IQ 2"),
    (791, 811, "IQ 3"),
    (814, 834, "IQ 4"),
    (837, 857, "IQ 5"),
    (1017, 1253, "MyTH4 1"),
    (1258, 1602, "FERM 1"),
    (1603, 1672, "SH3"),
    (1747, 1896, "MyTH4 2"),
    (1902, 2205, "FERM 2"),
    (632, 639, "Actin-binding"),
    (858, 935, "SAH"),
]

PMS2_DOMAINS = [
    (391, 552, "Disordered"),
]

RB1_DOMAINS = [
    (1, 42, "Disordered"),
    (373, 771, "Pocket; binds T and E1A"),
    (373, 579, "Domain A"),
    (580, 639, "Spacer"),
    (610, 632, "Disordered"),
    (640, 771, "Domain B"),
    (763, 928, "Interaction with LIMD1"),
    (771, 928, "Domain C; mediates interaction with E4F1"),
    (860, 928, "Disordered"),
]

RET_DOMAINS = [
    (168, 272, "Cadherin"),
    (724, 1016, "Protein kinase"),
    (29, 153, "Cadherin-like region 1 (CLD1)"),
    (265, 379, "Cadherin-like region 3 (CLD3)"),
    (405, 506, "Cadherin-like region 4 (CLD4)"),
]

SCN1A_DOMAINS = [
    (110, 454, "I"),
    (750, 1022, "II"),
    (1200, 1514, "III"),
    (1523, 1821, "IV"),
    (1915, 1944, "IQ"),
    (28, 60, "Disordered"),
    (455, 529, "Disordered"),
    (584, 627, "Disordered"),
    (1129, 1163, "Disordered"),
    (1561, 1571, "S1-S2 loop of repeat IV"),
    (1619, 1636, "S3b-S4 loop of repeat IV"),
    (1986, 2009, "Disordered"),
]

SCN2A_DOMAINS = [
    (111, 456, "I"),
    (741, 1013, "II"),
    (1190, 1504, "III"),
    (1513, 1811, "IV"),
    (1905, 1934, "IQ"),
    (28, 61, "Disordered"),
    (494, 529, "Disordered"),
    (590, 610, "Disordered"),
    (917, 918, "Binds SCN2B"),
    (1120, 1165, "Disordered"),
    (1935, 2005, "Disordered"),
]

SCN8A_DOMAINS = [
    (114, 442, "I"),
    (735, 1007, "II"),
    (1180, 1495, "III"),
    (1504, 1801, "IV"),
    (1895, 1924, "IQ"),
    (1, 20, "Disordered"),
    (28, 62, "Disordered"),
    (446, 530, "Disordered"),
    (568, 602, "Disordered"),
    (1107, 1148, "Disordered"),
    (1922, 1980, "Disordered"),
]

TSC2_DOMAINS = [
    (1531, 1758, "Rap-GAP"),
    (1, 400, "Required for interaction with TSC1"),
    (655, 676, "Disordered"),
    (930, 964, "Disordered"),
    (1083, 1132, "Disordered"),
    (1150, 1174, "Disordered"),
    (1331, 1352, "Disordered"),
    (1364, 1488, "Disordered"),
    (1765, 1793, "Disordered"),
]

VHL_DOMAINS = [
    (14, 18, "1"),
    (19, 23, "2"),
    (24, 28, "3"),
    (29, 33, "4"),
    (34, 38, "5"),
    (39, 43, "6"),
    (44, 48, "7"),
    (49, 53, "8"),
    (1, 65, "Disordered"),
    (14, 53, "8 X 5 AA tandem repeats of G-[PAVG]-E-E-[DAYSLE]"),
    (100, 155, "Involved in binding to CCT complex"),
    (157, 166, "Interaction with Elongin BC complex"),
]

DOMAIN_MAPS = {
    "TP53": TP53_DOMAINS,
    "BRCA1": BRCA1_DOMAINS,
    "PTEN": PTEN_DOMAINS,
    "ABCA4": ABCA4_DOMAINS,
    "ATP7A": ATP7A_DOMAINS,
    "COL1A1": COL1A1_DOMAINS,
    "COL2A1": COL2A1_DOMAINS,
    "COL5A1": COL5A1_DOMAINS,
    "ENG": ENG_DOMAINS,
    "HNF1A": HNF1A_DOMAINS,
    "LDLR": LDLR_DOMAINS,
    "MLH1": MLH1_DOMAINS,
    "MSH2": MSH2_DOMAINS,
    "MSH6": MSH6_DOMAINS,
    "MYBPC3": MYBPC3_DOMAINS,
    "MYO7A": MYO7A_DOMAINS,
    "PMS2": PMS2_DOMAINS,
    "RB1": RB1_DOMAINS,
    "RET": RET_DOMAINS,
    "SCN1A": SCN1A_DOMAINS,
    "SCN2A": SCN2A_DOMAINS,
    "SCN8A": SCN8A_DOMAINS,
    "TSC2": TSC2_DOMAINS,
    "VHL": VHL_DOMAINS,
}


def annotate_domain(gene: str, position: int) -> str:
    for start, end, label in DOMAIN_MAPS[gene]:
        if start <= position <= end:
            return label
    return "Unannotated / linker"
