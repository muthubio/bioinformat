#!/usr/bin/env python3
import os
import argparse
import pandas as pd

def extract_snp_positions(vcf_file, output_txt):
    if not os.path.exists(vcf_file):
        print(f"File not found: {vcf_file}")
        return None

    positions = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):  # Skip header lines
                continue
            columns = line.strip().split("\t")
            if len(columns) > 1:
                positions.append(columns[1])  # Extract position (2nd column)

    if not positions:
        print(f"No SNPs found in {vcf_file}")
        return None

    with open(output_txt, "w") as out_f:
        out_f.write("\n".join(positions))
    return output_txt

def str_conv(file):
    if not os.path.exists(file):
        print(f"File not found: {file}")
        return set()
    df = pd.read_csv(file, sep='\t', header=None)
    if df.empty:
        print(f"File is empty: {file}")
        return set()
    return set(df[0].astype(str).tolist())

def process_vcf(vcf_file, lineage_references):
    sample_name = os.path.splitext(os.path.basename(vcf_file))[0]
    snp_file = f"{sample_name}_snp.txt"

    
    extracted_snp_file = extract_snp_positions(vcf_file, snp_file)
    if not extracted_snp_file:
        print(f"No valid SNP file generated for {sample_name}. Skipping.")
        return None

    sample_snps = str_conv(extracted_snp_file)
    try:
        sample_snps = set(map(int, sample_snps))
    except ValueError:
        print(f"Error converting SNP positions to integers for {sample_name}. Skipping.")
        return None

    lineage_probabilities = {}
    for lineage, snp_set in lineage_references.items():
        matching_count = len(snp_set.intersection(sample_snps))
        probability = (matching_count / len(snp_set)) * 100
        lineage_probabilities[lineage] = f"{probability:.2f}%"

    max_probability = max([float(prob.strip('%')) for prob in lineage_probabilities.values()])
    if max_probability == 0:
        predicted_lineage = "unknown_lineage"
    else:
        predicted_lineage = max(lineage_probabilities, key=lambda k: float(lineage_probabilities[k].strip('%')))

    lineage_probabilities["Predicted_lineage"] = predicted_lineage
    lineage_probabilities["Sample_Name"] = sample_name
    return lineage_probabilities

def main():
    parser = argparse.ArgumentParser(
        description="Predict lineage from one or more VCF files using SNP matching probabilities."
    )
    parser.add_argument("--vcf_files", type=str, nargs="+", required=True,
                        help="Paths to one or more VCF files")
    parser.add_argument("--output_dir", type=str, default=".",
                        help="Directory to save output result files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    lineage_references = {
        'lineage1': {615938, 646531, 4398732, 272678, 3830566, 4081987, 3233605, 244550, 4081996, 1560912, 812502, 3876953, 344288, 811492, 3647591, 2508395, 763886, 495473, 2897528, 1119739, 4022652},
        'lineage2': {1834177, 4243460, 811753, 518987, 282892, 4236237, 4308395, 2298194, 497491, 4290135, 4050811, 4158493, 4254431},
        'lineage3': {2840999, 633562, 3582694, 24007},
        'lineage4': {1170404, 1341624, 3381641, 3482432},
        'lineage5': {1911301, 352646, 1505806, 2304017, 4339610, 1097633, 2744225, 3840932, 801959, 3882025, 1216822, 4399422, 4387392, 344258, 3603523, 2485956, 2516804, 910282, 4086604, 1648089, 9566, 2463455, 1648224, 1145442, 1578212, 1555432, 1256176, 1649265, 1799921, 1352566},
        'lineage6': {1372002, 3862148, 3213255, 4086697, 4236891, 334445, 2427828, 507989, 2847737, 1069146, 1867707, 3155164, 1886077},
        'lineage7': {3009027, 1125894, 1365895, 349081, 1561245, 4070302, 22303, 2305184, 1867937, 1513250, 4113054, 2303524, 3667883, 8876, 3603631, 2406193, 2086202, 2871227, 1297084, 4308411, 2414272, 2759363, 2833478, 3324231, 3473482, 2385615, 4262608, 2380244, 3182040, 3860696, 2913373, 1799774, 3225566, 1463776, 3635935, 2035938, 2414434, 497126, 3360233, 1137518, 1237743, 73456, 2383599, 2999030, 811642, 3470461, 2842238},
        'lineage8': {2022784, 1704707, 1665285, 1914246, 1941256, 4410891, 3645324, 4074512, 343314, 742270, 4227348, 3786517, 1391766, 4034711, 1505182, 3830815, 3840800, 2862497, 4210592, 4254755, 1594788, 508454, 3569319, 459561, 2415402, 1505455, 2940079, 17333, 1922231, 23098, 1446846, 2295359, 2466760, 765772, 2086736, 442577, 2939600, 3469397, 2450263, 3333720, 3151706, 1391967, 1805152, 3074400, 2440802, 2306915, 1858921, 595051, 1553399, 1227518},
        'lineage9': {3094577, 3225892, 3370805, 3414553},
        'bovis': {4038403, 354054, 652935, 2716806, 45193, 3371401, 3876953, 1647117, 2444952, 2643109, 4034981, 62768, 3022532, 1800521, 1462220, 3623372, 2782927, 437465, 2844761, 648667, 3723227, 4229470, 602207, 1137632, 3199071, 2767842, 3323617, 1344741, 1671658, 1813494, 3667193, 905082, 1492605},
        'caprae': {1673766, 422182, 649585, 2292409, 4398204},
        'dassie': {2271988, 557574, 1218256, 1226531, 1390170, 1554550, 1663500, 1751197, 2825199, 2846051, 2862458, 2871821, 2909119, 3062450, 3155107, 3392478, 3635392, 4085665, 4242313},
        'microtii': {3688579, 5671, 2874797, 1364877, 247252, 1004570},
        'mungi': {4268672, 2924039, 284041, 2022409, 24974, 437264, 3625360, 1256212, 270869, 64790, 2747554, 2999330, 2018219, 2437426, 2686388, 333878, 1216954, 570683, 4409917, 3180094, 73027, 1567814, 649036, 3326925, 45904, 4224851, 2937949, 2780127, 352481, 4274530, 3829350, 2993651, 3372917, 2427519},
        'orygis': {247300, 3337480, 1236745, 53899, 44812, 3061003, 1555342, 3598221, 3667352, 268953, 3039261, 748320, 3148454, 1344807, 3058989, 4039853, 2486576, 2298548, 3562936, 2296634, 2309820, 1579197, 3147196, 3155260, 4067010, 3241414, 3337675, 1449038, 459600, 587601, 1490258, 498771, 2301395, 3770449, 1918811, 1125468, 6109, 2429276, 8930, 1216738, 500710, 1652839, 634348, 1810542, 4069621, 1383800, 1834363, 4409725},
        'pinepedii': {3337219, 1036936, 2301448, 736140, 2352292, 4161188, 458156, 1937070, 1648564, 2301749, 28344, 1650257, 1650258, 1924695, 3090780, 3145957, 1391466, 2942959, 918007, 348280, 458617, 2908666, 1577723, 3469951},
        'surricate': {2914435, 1883911, 3343629, 3396110, 4049678, 3060380, 21664, 2428451, 4396201, 1831339, 2762930, 1366453, 2410438, 344135, 3233351, 2716491, 765004, 495566, 741841, 3233747, 2746981, 2872678, 1650665, 2069610, 763375, 3603697, 286579, 351101}
    }
    all_results = []
    for vcf_file in args.vcf_files:
        print(f"Processing {vcf_file} ...")
        result = process_vcf(vcf_file, lineage_references)
        if result:
            all_results.append(result)
    if all_results:
        df = pd.DataFrame(all_results).set_index("Sample_Name")
        combined_file = os.path.join(args.output_dir, "combined_lineage_results.txt")
        df.to_csv(combined_file, sep="\t")
        print(f"\nCombined results saved to {combined_file}")

if __name__ == "__main__":
    main()

